# -*- coding: utf-8 -*-
"""
BM-01 一维圣维南（SVE）模型 —— 自包含单文件版
================================================
把原项目里三个文件（bm01_sve1d_compare.py 的 v1 求解器 + bm01_sve1d_v2.py 的
断面开关与修正）合并成一个文件，只依赖 Python 标准库（无 numpy / matplotlib）。

解什么：一维圣维南（浅水）方程，守恒形式 + Rusanov（局部 Lax-Friedrichs）通量 + 显式欧拉：
    ∂A/∂t + ∂Q/∂x = 0
    ∂Q/∂t + ∂(Q²/A + gA²/2B)/∂x = gA(S0 − Sf)

算例：一个 6 km × 50 m 的渠池，目标水位 7 m，来水 180 m³/s，t=900 s 时上游
阶跃 +45 m³/s（约 +25%）。下游闸门由一个 PI 控制器按"闸前水位"调节开度。
看池均/闸前水位有没有超过 ±150 mm 正常调节条带。

用法：
    python bm01_sve1d_standalone.py                 # 跑矩形基准 + 梯形开关（240/480 格）
    python bm01_sve1d_standalone.py --selftest      # 自检：单点传感器下复现 v1 记录值
    python bm01_sve1d_standalone.py --grids 120,240 # 更快，只看两级网格
    python bm01_sve1d_standalone.py --only S2c      # 只跑某一个开关

证据等级：synthetic-1d-sve（数值降阶，量级示例参数，非现场结论）。
"""
from __future__ import annotations

import argparse
import json
import math
import time

# ======================================================================
# 全局常量
# ======================================================================
G = 9.81
MANNING_N = 0.014          # 混凝土衬砌量级示例
CFL = 0.40
# 候选时间步：都能整除 14400 s（总时长）、900 s（阶跃时刻）、300 s（控制周期）
DT_CANDIDATES = [5, 4, 3.6, 3, 2.5, 2.25, 2, 1.8, 1.5, 1.25, 1.125, 1, 0.9, 0.75,
                 0.6, 0.5, 0.45, 0.4, 0.36, 0.3, 0.25, 0.225, 0.2]


# ======================================================================
# PI 控制器（与主案例 bm01_pool_step.py 逐行等价）
# ======================================================================
class Ctrl:
    def __init__(self, p: dict):
        self.p = p
        self.a = p["a0"]
        self.a_cmd = p["a0"]
        self.integral = 0.0
        self.last_ctrl_t = -1e18

    def control(self, t: float, h: float) -> None:
        p = self.p
        if t - self.last_ctrl_t < p["ctrl_period_s"]:
            return
        self.last_ctrl_t = t
        meas = round(h / p["meas_resolution_m"]) * p["meas_resolution_m"]
        e = meas - p["h_target"]
        if abs(e) >= p["deadband_m"]:
            at_hi = self.a >= 1.0 - 1e-12 and e > 0
            at_lo = self.a <= 1e-12 and e < 0
            if not (at_hi or at_lo):
                self.integral += e * p["ctrl_period_s"]
        u = p["kp"] * e + p["ki"] * self.integral
        self.a_cmd = min(max(p["a0"] + u, 0.0), 1.0)

    def step_gate(self, dt: float) -> None:
        step = self.p["rate_limit_per_s"] * dt
        delta = self.a_cmd - self.a
        if abs(delta) < step:
            self.a = self.a_cmd
        else:
            self.a += step if delta > 0 else -step
        self.a = min(max(self.a, 0.0), 1.0)


# ======================================================================
# 断面几何：矩形（m=0）/ 梯形（m>0）
# ======================================================================
class Section:
    """A(h)、B(h)、P(h)、静水压力积分 I1(h)=∫0^h (h−η)B(η)dη。"""

    def __init__(self, b: float, m: float = 0.0):
        self.b = b        # 底宽
        self.m = m        # 边坡（水平:垂直），m=0 → 矩形

    def A(self, h):
        return h * (self.b + self.m * h)

    def B(self, h):
        return self.b + 2.0 * self.m * h

    def P(self, h):
        return self.b + 2.0 * h * math.sqrt(1.0 + self.m * self.m)

    def I1(self, h):
        return 0.5 * self.b * h * h + (self.m / 3.0) * h * h * h

    def is_rect(self):
        return abs(self.m) < 1e-12


def normal_slope(p: dict, n: float = MANNING_N) -> float:
    """矩形断面：反解渠底坡度，让 (h_target, q0) 恰为均匀流（正常水深）。"""
    B, ht = p["pool_width_m"], p["h_target"]
    A = B * ht
    R = A / (B + 2.0 * ht)
    return (n * p["q0"] / (A * R ** (2.0 / 3.0))) ** 2


def section_normal_slope(sec: "Section", h: float, Q: float, n: float) -> float:
    """按**实际断面**反解纵坡，使 (h, Q) 恰为均匀流（梯形断面下必须用本式）。"""
    A = sec.A(h)
    P = sec.P(h)
    R = A / P if P > 0 else h
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


def pick_dt(dx: float, wmax: float) -> float:
    lo = CFL * dx / max(wmax, 1e-6)
    for c in DT_CANDIDATES:
        if c <= lo:
            return c
    raise ValueError(f"CFL 太紧：dx={dx} 需 dt<{lo:.3f} s")


# ======================================================================
# 核心求解器（矩形与梯形通用；m_slope=0 + ctrl_on="gate" 时与 v1 逐行等价）
# ======================================================================
def sve_1d(p: dict, ncell: int, step_dq: float, *,
           ctrl_on: str = "gate", m_slope: float = 0.0,
           n: float = MANNING_N, cfl: float = CFL,
           s0: float = None, sensor_len_m: float = 100.0,
           early_window_s: float = 7200.0) -> dict:
    """
    ctrl_on:
        "gate"     = 单点传感器（最末单元 h[-1]）—— 回路增益 ∝ 1/dx，
                     网格越细越容易在 300 s 采样下失稳（960 格出现增幅极限环）。
        "gate_avg" = 闸前固定物理长度 sensor_len_m 内平均（回路增益与 dx 无关）。
        "mean"     = 全池平均（对照组）。
        "off"      = 不开控制器。
    """
    B0, L, ht = p["pool_width_m"], p["pool_len_m"], p["h_target"]
    sec = Section(B0, m_slope)
    dx = L / ncell
    S0 = normal_slope(p, n) if s0 is None else s0
    h = [ht] * ncell
    Q = [p["q0"]] * ncell
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))

    wmax = p["q0"] / (B0 * ht) + math.sqrt(G * ht)
    lo = cfl * dx / wmax
    dt = pick_dt(dx, wmax) if abs(cfl - CFL) < 1e-15 else None
    if dt is None:
        dt = next(c for c in DT_CANDIDATES if c <= lo)
    nt = int(round(p["duration_s"] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = Ctrl(p)

    out = {"t": [], "dev_mean": [], "dev_gate": [], "dev_head": [], "gate": [], "q_out": []}
    peak_mean = peak_gate = 0.0
    early_peak = 0.0
    tail_mean = 0.0
    t_tail = p["duration_s"] - p["settle_window_s"]
    over_band_dt = 0
    vol_net = 0.0
    c1 = dt / dx
    sqrt = math.sqrt
    hmin = 1e-3
    friction_correct = not sec.is_rect()

    for k in range(nt):
        tk = k * dt
        q_up = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)
        q_gate = p["kg"] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)

        if ctrl_on != "off":
            if ctrl_on == "gate":
                meas = h[-1]
            elif ctrl_on == "gate_avg":
                meas = sum(h[ncell - nsens:]) / nsens
            else:
                meas = sum(h) / ncell
            ctrl.control(tk, meas)
            ctrl.step_gate(dt)

        # ---- Rusanov（局部 Lax-Friedrichs）数值通量 ----
        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif
        F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = hx[j], hx[j + 1]
            if h1 < hmin: h1 = hmin
            if h2 < hmin: h2 = hmin
            A1, A2 = sec.A(h1), sec.A(h2)
            if A1 < 1e-6: A1 = 1e-6
            if A2 < 1e-6: A2 = 1e-6
            Bc1, Bc2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / Bc1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / Bc2)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            P1 = 0.5 * G * A1 * A1 / Bc1
            P2 = 0.5 * G * A2 * A2 / Bc2
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2 + P1 + P2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate

        # ---- 单元推进（显式欧拉）----
        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            A = sec.A(hi)
            if A < 1e-6: A = 1e-6
            Qi = Q[i]
            if friction_correct:
                R = A / sec.P(hi) if sec.P(hi) > 0 else hi
            else:
                R = A / (B0 + 2.0 * hi)
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        vol_net += (q_up - q_gate) * dt

        dm = sum(h) / ncell - ht
        dg = h[-1] - ht
        if abs(dm) > peak_mean: peak_mean = abs(dm)
        if abs(dg) > peak_gate: peak_gate = abs(dg)
        if tk <= early_window_s and abs(dm) > early_peak: early_peak = abs(dm)
        if tk >= t_tail and abs(dm) > tail_mean: tail_mean = abs(dm)
        if abs(dm) > p["band_normal_m"]: over_band_dt += 1
        if k % rec_every == 0:
            out["t"].append(tk); out["dev_mean"].append(dm); out["dev_gate"].append(dg)
            out["dev_head"].append(h[0] - ht); out["gate"].append(ctrl.a)
            out["q_out"].append(q_gate)

    # 质量守恒：净入流 vs 池内蓄量增量（精确式，梯形下不要用 B·Δh）
    dv = 0.0
    for x in h:
        dv += (x - ht) * (sec.b + sec.m * (x + ht)) * dx
    return {"ncell": ncell, "dx_m": dx, "dt_s": dt, "steps": nt, "S0": S0,
            "ctrl_on": ctrl_on, "m_slope": m_slope, "sensor_len_m": sensor_len_m,
            "nsens_cells": nsens, "early_window_s": early_window_s,
            "series": out,
            "peak_mean_mm": peak_mean * 1000, "peak_gate_mm": peak_gate * 1000,
            "early_peak_mm": early_peak * 1000,
            "final_mean_mm": tail_mean * 1000,
            "over_band_s": over_band_dt * dt,
            "vol_net_m3": vol_net, "vol_storage_m3": dv,
            "mass_balance_rel_err": abs(dv - vol_net) / max(abs(vol_net), 1.0)}


# ======================================================================
# 默认算例参数（与 raw/bm01_main.json 一致；如需换参数直接改这里或传 --params）
# ======================================================================
DEFAULT_PARAMS = {
    "pool_len_m": 6000.0, "pool_width_m": 50.0,
    "q0": 180.0, "h_target": 7.0, "a0": 0.5, "kg": 360.0,
    "kp": 1.6, "ki": 0.0016,
    "ctrl_period_s": 300.0, "deadband_m": 0.01, "meas_resolution_m": 0.001,
    "rate_limit_per_s": 0.0003333333333333333,
    "dt": 5.0, "duration_s": 14400.0, "t_step_s": 900.0,
    "band_normal_m": 0.15, "band_alarm_m": 0.3,
    "settle_tol_m": 0.02, "settle_window_s": 1800.0,
    "area_m2": 300000.0,
    "gate_up_capacity_m3s": 180.0, "gate_down_capacity_m3s": 180.0,
}
DEFAULT_STEP_DQ = 45.0

SWITCHES = [
    ("S2a", "梯形断面 边坡 m=1.0", dict(m_slope=1.0)),
    ("S2b", "梯形断面 边坡 m=1.5", dict(m_slope=1.5)),
    ("S2c", "梯形断面 边坡 m=2.5", dict(m_slope=2.5)),
]
SENSOR_LEN_M = 100.0


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 一维圣维南（自包含单文件）")
    ap.add_argument("--grids", default="240,480")
    ap.add_argument("--only", default=None)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--params", default=None, help="可选：参数 JSON 文件路径，默认用内置参数")
    args = ap.parse_args()

    p = DEFAULT_PARAMS
    step_dq = DEFAULT_STEP_DQ
    if args.params:
        with open(args.params, encoding="utf-8") as fh:
            raw = json.load(fh)
        p = raw.get("params", raw)
        step_dq = float(raw.get("scenario", {}).get("step_dq", DEFAULT_STEP_DQ))
    grids = [int(x) for x in args.grids.split(",")]
    t0 = time.time()

    c = math.sqrt(G * p["h_target"])
    print("BM-01 一维圣维南（自包含单文件）")
    print(f"池长 {p['pool_len_m']:.0f} m × 底宽 {p['pool_width_m']:.0f} m，"
          f"h_target={p['h_target']:.2f} m，Q0={p['q0']:.0f} m³/s，"
          f"阶跃 +{step_dq:.0f} m³/s，条带 ±{p['band_normal_m']*1000:.0f} mm")
    print(f"波速 c={c:.2f} m/s，L/c={p['pool_len_m']/c/60:.1f} min，"
          f"2L/c={2*p['pool_len_m']/c/60:.1f} min\n")

    if args.selftest:
        v1ref = {"120": 142.88143604489179, "240": 156.6079038724011}
        print("[自检] 与 v1 记录值比对（矩形 + 单点传感器）")
        ok = True
        for nc in (120, 240):
            r = sve_1d(p, nc, step_dq, ctrl_on="gate", m_slope=0.0)
            d = abs(r["peak_mean_mm"] - v1ref[str(nc)])
            good = d < 0.5
            ok = ok and good
            print(f"   {nc:>4} 格  本文件 {r['peak_mean_mm']:.2f} mm vs "
                  f"v1 {v1ref[str(nc)]:.2f} mm　差 {d:.3f} mm → {'一致 ✓' if good else '不一致 ✗'}")
        print(f"   结论：{'与 v1 逐行等价 ✓' if ok else '存在差异 ✗'}")
        return 0 if ok else 1

    print(f"[基准] 矩形 m=0；传感器 = 闸前 {SENSOR_LEN_M:.0f} m 平均（网格无关）")
    base = {}
    for nc in grids:
        base[nc] = sve_1d(p, nc, step_dq, ctrl_on="gate_avg", m_slope=0.0,
                          sensor_len_m=SENSOR_LEN_M)
        r = base[nc]
        print(f"   {nc:>4} 格  dt={r['dt_s']:>6} s｜池均峰值 {r['peak_mean_mm']:7.1f} mm｜"
              f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰值 {r['peak_gate_mm']:7.1f} mm｜"
              f"末段 {r['final_mean_mm']:7.1f} mm｜守恒 {r['mass_balance_rel_err']:.1e}")
    base_conv = base[grids[-1]]["early_peak_mm"]
    spread = abs(base[grids[-1]]["early_peak_mm"] - base[grids[-2]]["early_peak_mm"])
    print(f"   收敛峰值（前2h窗，{grids[-1]} 格）= {base_conv:.1f} mm"
          f"（与 {grids[-2]} 格差 {spread:.1f} mm，{spread/base_conv*100:.1f}%）")

    only = args.only.split(",") if args.only else None
    for tag, name, kw in SWITCHES:
        if only and tag not in only:
            continue
        print(f"\n[{tag}] {name}")
        sec = Section(p["pool_width_m"], kw.get("m_slope", 0.0))
        s0_sec = section_normal_slope(sec, p["h_target"], p["q0"], MANNING_N)
        s0_leg = normal_slope(p)
        print(f"   纵坡：断面感知 {s0_sec:.6e}（1/{1/s0_sec:,.0f}）｜"
              f"矩形反解 {s0_leg:.6e}（1/{1/s0_leg:,.0f}）　矩形式偏大 {s0_leg/s0_sec:.2f} 倍")
        for nc in grids:
            r = sve_1d(p, nc, step_dq, ctrl_on="gate_avg", s0=s0_sec,
                       sensor_len_m=SENSOR_LEN_M, **kw)
            print(f"   {nc:>4} 格  dt={r['dt_s']:>6} s｜池均峰值 {r['peak_mean_mm']:7.1f} mm｜"
                  f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰值 {r['peak_gate_mm']:7.1f} mm｜"
                  f"末段 {r['final_mean_mm']:6.1f} mm｜越带 {r['over_band_s']/60:5.2f} min｜"
                  f"守恒 {r['mass_balance_rel_err']:.1e}")
        conv = sve_1d(p, grids[-1], step_dq, ctrl_on="gate_avg", s0=s0_sec,
                      sensor_len_m=SENSOR_LEN_M, **kw)["early_peak_mm"]
        print(f"   收敛峰值（前2h窗，{grids[-1]} 格）= {conv:.1f} mm　"
              f"对基准 {base_conv:.1f} mm 的增量 {conv-base_conv:+.1f} mm"
              f"（{(conv-base_conv)/base_conv*100:+.1f}%）")

    print(f"\n总用时 {time.time()-t0:.1f} s")
    print("证据等级：synthetic-1d-sve（数值降阶，非现场结论）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
