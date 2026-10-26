# -*- coding: utf-8 -*-
"""
BM-01 外部门锚 v1：用 EPA SWMM（独立第三方 1D 动力波求解器）复建 BM-01 算例
================================================================================
三模型第四轮的"一票否决语"：
  "判定码的翻转会不会只是你自家求解器在闸门边界条件下的伪影？"
排除它唯一零现场成本的办法：用**另一个独立开发、独立验收的 1D 圣维南求解器**重跑同一算例。
本脚本用 EPA SWMM 5.2（经 pyswmm / swmm-toolkit 调用，**非本仓库任何代码**）。

复建原则（可比较性）
--------------------
1. 几何：棱柱矩形 B=50 m，L 可调，底坡由"正常水深反解"（与 bm01_sve1d_compare.normal_slope 同式）。
2. 池体离散：SWMM 每段 conduit 为一个计算元 ⇒ 用 ncell 段串联提供空间分辨率。
3. 闸门：SWMM ORIFICE（RECT_CLOSED），令其流量律与仓库闸门律**逐点相等**：
     仓库  q_gate = kg·a·√(h/h_t) = (kg/√h_t)·a·√h
     SWMM  q_orif = Cd·(a·W·H)·√(2g·h) = Cd·W·H·√(2g)·a·√h
   ⇒ Cd·W·H = kg/(√h_t·√(2g)) = 30.72（BM-01）；取 Cd=0.6, W=50 ⇒ H=1.024 m。
4. 反馈量：
     mode="pool" → 池均等效水位 = Σ(网段体积)/(B·L)（= 仓库"蓄量等效"估计量；棱柱池下 = 全池空间平均）
     mode="gate" → 闸前节点水位 depth（= 仓库 h[-1] 单点）
5. 控制器：与 bm01_sve1d_compare.Ctrl 逐行等价，在 Python 中实现，每路由步写 orifice.target_setting。

口径：纯数值外部对照，证据等级 external-1d-swmm；不改变 BM-01 判定码，不构成工程结论。
"""

from __future__ import annotations

import json
import math
import os
import sys
import time

# pyswmm 1.0.1 依赖已被 Py3.12+ 移除的 distutils；setuptools 提供其兼容 shim。
# 必须先 import setuptools（在 PYTHONPATH 方式下 distutils-precedence.pth 不会被处理），
# 否则 `import pyswmm` 会因 ModuleNotFoundError: distutils 失败。
try:
    import setuptools  # noqa: F401
except Exception:
    pass

G = 9.81

P = {
    "pool_len_m": 6000.0,
    "pool_width_m": 50.0,
    "q0": 180.0,
    "h_target": 7.0,
    "a0": 0.5,
    "kg": 360.0,
    "kp": 1.6,
    "ki": 0.0016,
    "ctrl_period_s": 300.0,
    "deadband_m": 0.01,
    "meas_resolution_m": 0.001,
    "rate_limit_per_s": 1.0 / 3000.0,
    "band_normal_m": 0.15,
    "step_dq": 45.0,
    "n_manning": 0.014,
}


def normal_slope(B, ht, Q, n):
    A = B * ht
    R = A / (B + 2.0 * ht)
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


class Ctrl:
    """与 bm01_sve1d_compare.Ctrl 逐行等价。"""

    def __init__(self):
        self.a = P["a0"]
        self.a_cmd = P["a0"]
        self.integral = 0.0
        self.last_ctrl_t = -1e18

    def control(self, t, h):
        if t - self.last_ctrl_t < P["ctrl_period_s"]:
            return
        self.last_ctrl_t = t
        meas = round(h / P["meas_resolution_m"]) * P["meas_resolution_m"]
        e = meas - P["h_target"]
        if abs(e) >= P["deadband_m"]:
            at_hi = self.a >= 1.0 - 1e-12 and e > 0
            at_lo = self.a <= 1e-12 and e < 0
            if not (at_hi or at_lo):
                self.integral += e * P["ctrl_period_s"]
        u = P["kp"] * e + P["ki"] * self.integral
        self.a_cmd = min(max(P["a0"] + u, 0.0), 1.0)

    def step_gate(self, dt):
        step = P["rate_limit_per_s"] * dt
        delta = self.a_cmd - self.a
        if abs(delta) < step:
            self.a = self.a_cmd
        else:
            self.a += step if delta > 0 else -step
        self.a = min(max(self.a, 0.0), 1.0)


def _hhmmss(sec):
    sec = float(sec)
    h = int(sec // 3600)
    m = int((sec - 3600 * h) // 60)
    s = sec - 3600 * h - 60 * m
    return "%02d:%02d:%06.3f" % (h, m, s)


def build_inp(path, ncell, dx, s0, t_spinup_s, t_step_s, t_end_s,
              routing_step_s=5.0, xs_depth=30.0, z0=100.0):
    B = P["pool_width_m"]
    kd = P["kg"] / (math.sqrt(P["h_target"]) * math.sqrt(2.0 * G))
    Cd = 0.6
    W = B
    H = kd / (Cd * W)
    L = []

    def add(*a):
        L.append(" ".join(str(x) for x in a) if len(a) > 1 else a[0])

    add("[TITLE]")
    add("BM-01 external anchor via EPA SWMM (independent 1D dynamic-wave solver)")
    add("")
    add("[OPTIONS]")
    add("FLOW_UNITS           CMS")
    add("FLOW_ROUTING         DYNWAVE")
    add("START_DATE           01/01/2000")
    add("START_TIME           00:00:00")
    add("REPORT_START_DATE    01/01/2000")
    add("REPORT_START_TIME    00:00:00")
    add("END_DATE             01/01/2000")
    add("END_TIME             %s" % _hhmmss(t_end_s))
    add("REPORT_STEP          %s" % _hhmmss(routing_step_s))
    add("WET_STEP             %s" % _hhmmss(routing_step_s))
    add("DRY_STEP             %s" % _hhmmss(routing_step_s))
    add("ROUTING_STEP         %g" % routing_step_s)
    add("RULE_STEP            %s" % _hhmmss(routing_step_s))
    add("INERTIAL_DAMPING     PARTIAL")
    add("NORMAL_FLOW_LIMITED  BOTH")
    add("FORCE_MAIN_EQUATION  H-W")
    add("VARIABLE_STEP        0.0")
    add("LENGTHENING_STEP     0")
    add("MIN_SURFAREA         12.557")
    add("MAX_TRIALS           50")
    add("HEAD_TOLERANCE       0.0015")
    add("SYS_FLOW_TOL         5")
    add("LAT_FLOW_TOL         5")
    add("MINIMUM_STEP         0.5")
    add("THREADS              1")
    add("")
    add("[EVAPORATION]")
    add("CONSTANT         0.0")
    add("DRY_ONLY         NO")
    add("")
    add("[JUNCTIONS]")
    add(";;Name           Elevation  MaxDepth  InitDepth  SurDepth  Aponded")
    for i in range(ncell + 1):
        z = z0 - s0 * (i * dx)
        add("J%-6d  %10.4f  %8.2f  %8.2f  %8.2f  %8.2f"
            % (i, z, xs_depth, P["h_target"], 0.0, 0.0))
    add("")
    add("[OUTFALLS]")
    add(";;Name           Elevation  Type       Stage Data       Gated    RouteTo")
    add("OUT              %10.4f  FREE                        NO" % (z0 - s0 * (ncell * dx) - 0.05))
    add("")
    add("[CONDUITS]")
    add(";;Name           From     To       Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow")
    for i in range(ncell):
        add("C%-6d  J%-6d  J%-6d  %9.3f  %9.4f  %9.2f  %9.2f  %9.2f  %9.2f"
            % (i, i, i + 1, dx, P["n_manning"], 0.0, 0.0, P["q0"], 0.0))
    add("")
    add("[XSECTIONS]")
    add(";;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels")
    for i in range(ncell):
        add("C%-6d  RECT_OPEN    %12.4f    %8.2f    0          0          1" % (i, xs_depth, B))
    add("")
    add("[ORIFICES]")
    add(";;Name           From     To       Type       Offset     Qcoeff     Gated      CloseTime")
    add("GATE             J%-6d  OUT      SIDE       %8.3f   %8.3f   NO         0.0" % (ncell, 0.0, Cd))
    add("")
    add("[XSECTIONS]")
    add(";;Link           Shape         Geom1            Geom2      Geom3      Geom4      Barrels")
    add("GATE             RECT_CLOSED   %12.4f    %8.2f    0          0" % (H, W))
    add("")
    add("[TIMESERIES]")
    add(";;Name           Date       Time       Value")
    add("TSIN             01/01/2000 00:00:00   %.3f" % P["q0"])
    add("TSIN             01/01/2000 %s   %.3f" % (_hhmmss(t_step_s), P["q0"]))
    add("TSIN             01/01/2000 %s   %.3f" % (_hhmmss(t_step_s + 1.0), P["q0"] + P["step_dq"]))
    add("TSIN             01/01/2000 %s   %.3f" % (_hhmmss(t_end_s), P["q0"] + P["step_dq"]))
    add("")
    add("[INFLOWS]")
    add(";;Node           Parameter        Time Series      Type       Mfactor  Sfactor  Baseline  Pattern")
    add("J0               FLOW             TSIN             FLOW       1.0      1.0")
    add("")
    add("[REPORT]")
    add("INPUT      NO")
    add("CONTROLS   NO")
    add("SUBCATCHMENTS ALL")
    add("NODES ALL")
    add("LINKS ALL")
    add("")
    add("[TAGS]")
    add("")

    with open(path, "w", encoding="ascii") as fh:
        fh.write("\n".join(L))
    return {"dx": dx, "ncell": ncell, "Cd": Cd, "W": W, "H": H, "kd_target": kd,
            "s0": s0, "t_spinup_s": t_spinup_s, "t_step_s": t_step_s, "t_end_s": t_end_s,
            "routing_step_s": routing_step_s}


def run(path, ncell, mode="pool", routing_step_s=5.0, t_spinup_s=7200.0,
        t_step_s=7260.0, verbose=True):
    from pyswmm import Simulation, Nodes, Links

    B, L_pool, ht = P["pool_width_m"], P["pool_len_m"], P["h_target"]
    cond_ids = ["C%d" % i for i in range(ncell)]
    gate_node = "J%d" % ncell
    ctrl = Ctrl()
    agg = {"peak_pool_mm": 0.0, "peak_gate_mm": 0.0, "base_pool_mm": None,
           "base_gate_mm": None, "a_final": None, "max_pool_mm": -1e9, "min_pool_mm": 1e9}
    csv_rows = []
    t0 = time.time()
    n_step = 0
    with Simulation(path, reportfile=path + ".rpt", outputfile=path + ".out") as sim:
        nodes = Nodes(sim)
        links = Links(sim)
        gate = links["GATE"]
        gate.target_setting = P["a0"]
        t_start = sim.start_time
        last_tk = -1.0
        for step in sim:
            # 真实仿真时钟（SWMM 可能自动缩步，不能假设恒定路由步）
            tk = (sim.current_time - t_start).total_seconds()
            dt_eff = tk - last_tk if last_tk >= 0 else routing_step_s
            last_tk = tk
            vol = 0.0
            for cid in cond_ids:
                vol += links[cid].volume
            pool_mean = vol / (B * L_pool)
            pool_mm = (pool_mean - ht) * 1000.0
            gate_mm = (nodes[gate_node].depth - ht) * 1000.0
            meas = pool_mean if mode == "pool" else nodes[gate_node].depth
            ctrl.control(tk, meas)
            ctrl.step_gate(dt_eff)
            gate.target_setting = ctrl.a
            if tk >= t_spinup_s:
                if agg["base_pool_mm"] is None:
                    agg["base_pool_mm"] = pool_mm
                    agg["base_gate_mm"] = gate_mm
                dp = pool_mm - agg["base_pool_mm"]
                dg = gate_mm - agg["base_gate_mm"]
                if abs(dp) > agg["peak_pool_mm"]:
                    agg["peak_pool_mm"] = abs(dp)
                if abs(dg) > agg["peak_gate_mm"]:
                    agg["peak_gate_mm"] = abs(dg)
                agg["max_pool_mm"] = max(agg["max_pool_mm"], pool_mm)
                agg["min_pool_mm"] = min(agg["min_pool_mm"], pool_mm)
                csv_rows.append((tk, pool_mm, gate_mm, ctrl.a, gate.flow))
            n_step += 1
            if verbose and n_step % 1000 == 0:
                print("    step %6d  t=%.0f s  a=%.3f  pool=%+.1f mm  gate=%+.1f mm"
                      % (n_step, tk, ctrl.a, pool_mm, gate_mm))
    agg["a_final"] = ctrl.a
    agg["wall_s"] = time.time() - t0
    agg["n_step"] = n_step
    agg["ncell"] = ncell
    # 落盘 CSV
    csvp = path.replace(".inp", "_%s.csv" % mode)
    with open(csvp, "w", encoding="utf-8") as fh:
        fh.write("t_s,pool_mm,gate_mm,a,q_out\n")
        for r in csv_rows:
            fh.write("%.1f,%.4f,%.4f,%.6f,%.4f\n" % r)
    return {"agg": agg, "csv": csvp}


def main():
    import argparse
    import shutil
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncell", type=int, default=20)
    ap.add_argument("--mode", default="pool", choices=["pool", "gate"])
    ap.add_argument("--spinup", type=float, default=7200.0)
    ap.add_argument("--horizon", type=float, default=14400.0)
    ap.add_argument("--routing", type=float, default=5.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    B, ht, q0 = P["pool_width_m"], P["h_target"], P["q0"]
    s0 = normal_slope(B, ht, q0, P["n_manning"])
    dx = P["pool_len_m"] / args.ncell
    here = os.path.dirname(os.path.abspath(__file__))
    # SWMM 的 C 引擎无法打开含非 ASCII 的路径 ⇒ 工作文件放纯 ASCII 临时目录
    workdir = os.path.join(os.environ.get("LOCALAPPDATA", os.path.expanduser("~")),
                           "Temp", "bm01_ext_swmm")
    os.makedirs(workdir, exist_ok=True)
    tag = "%s_n%d_r%g" % (args.mode, args.ncell, args.routing)
    inp = os.path.join(workdir, "bm01_ext_%s.inp" % tag)
    t_step_s = args.spinup + 60.0
    t_end_s = t_step_s + args.horizon
    meta = build_inp(inp, args.ncell, dx, s0, t_spinup_s=args.spinup, t_step_s=t_step_s,
                     t_end_s=t_end_s, routing_step_s=args.routing)

    print("SWMM 外部对照：ncell=%d  dx=%.1f m  S0=%.4e  Cd=%.3f W=%.1f H=%.4f  mode=%s"
          % (args.ncell, dx, s0, meta["Cd"], meta["W"], meta["H"], args.mode))
    c = math.sqrt(G * ht)
    print("波速 c=%.2f m/s  L/c=%.1f min  路由步=%.1f s  总时长=%.0f s"
          % (c, P["pool_len_m"] / c / 60.0, args.routing, t_end_s))
    print("工作目录（ASCII）：%s" % workdir)
    res = run(inp, args.ncell, mode=args.mode, routing_step_s=args.routing,
              t_spinup_s=args.spinup, t_step_s=t_step_s)
    a = res["agg"]
    print("\n结果：")
    print("  池均基线 %s mm → 峰值(相对基线) %.2f mm   区间[%.1f, %.1f]"
          % ("%.2f" % a["base_pool_mm"], a["peak_pool_mm"], a["min_pool_mm"], a["max_pool_mm"]))
    print("  闸前基线 %s mm → 峰值(相对基线) %.2f mm"
          % ("%.2f" % a["base_gate_mm"], a["peak_gate_mm"]))
    print("  末开度 a=%.3f  步数=%d  墙钟=%.1f s" % (a["a_final"], a["n_step"], a["wall_s"]))
    # 结果拷回案例目录
    dst_csv = os.path.join(here, "raw", "bm01_ext_swmm_%s.csv" % tag)
    try:
        shutil.copyfile(res["csv"], dst_csv)
        print("  CSV → %s" % dst_csv)
    except Exception as e:
        print("  CSV 拷贝失败：%s（原文件 %s）" % (e, res["csv"]))
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump({"meta": meta, "params": P, "agg": a, "csv": dst_csv}, fh, ensure_ascii=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
