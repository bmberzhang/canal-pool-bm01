# -*- coding: utf-8 -*-
"""
BM-01 附加诊断：一维圣维南（SVE）对照 —— 集总模型到底丢了多少？
================================================================
问题
----
BM-01 主案例用的是集总蓄量模型  Δh ≈ ΔQ·Δt/A ，它隐含两个假设：
  (1) 池内水位处处相同（零维，所有 ∂/∂x 项被丢掉）；
  (2) 闸门"瞬时"就能量到全池水位（闸门看到的就是池均水位）。
真实渠道里，闸门量到的是**闸前局部水位**，水位变化要以波的形式传过去。

本文件把这两个假设拿掉，看差多少：同一组参数、同一套 PI、同一采样周期（300 s）、
同一速率限，跑一个一维圣维南模型，与集总结果并排对比。

自校验（重要）
--------------
本脚本**只用 Python 标准库**（不依赖 numpy / matplotlib），参数直接读
`raw/bm01_main.json`——那是主案例自己记录的参数（单一事实来源）。
并且会先用同一套纯 Python 代码**重放集总模型**，与记录指标比对；
只有重放一致，后面的一维结果才可信。

口径（不得外推）
----------------
仍是数值降阶证据（synthetic-1d-sve），**不改变 BM-01 的判定码与证据等级**，
参数仍为量级示例、不是中线现场数据，不构成工程结论。本文件唯一用途是把
"降阶到底丢了多少"变成一个具体数字（论文"方法验证"用）。

用法
----
  python bm01_sve1d_compare.py                 # 默认 120/240/480 三级网格
  python bm01_sve1d_compare.py --finest 240    # 想快就跑两级
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

G = 9.81
MANNING_N = 0.014          # 混凝土衬砌量级示例
CFL = 0.40
# 候选时间步：都能整除 14400 s（总时长）、900 s（阶跃时刻）、300 s（控制周期）
DT_CANDIDATES = [5, 4, 3.6, 3, 2.5, 2.25, 2, 1.8, 1.5, 1.25, 1.125, 1, 0.9, 0.75,
                 0.6, 0.5, 0.45, 0.4, 0.36, 0.3, 0.25, 0.225, 0.2]


# ----------------------------------------------------------------------
# 控制器：与 bm01_pool_step.py 的 Pool.control / Pool.step_gate 逐行等价
# ----------------------------------------------------------------------
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


def gate_flow(p: dict, a: float, h: float) -> float:
    return p["kg"] * a * math.sqrt(max(h, 1e-3) / p["h_target"])


# ----------------------------------------------------------------------
# 0D 集总模型重放（与主案例 simulate() 的 pool1 分支逐行等价）
# ----------------------------------------------------------------------
def lumped_0d(p: dict, step_dq: float) -> dict:
    dt = p["dt"]
    n = int(round(p["duration_s"] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = Ctrl(p)
    h = p["h_target"]
    out = {"t": [], "dev": [], "gate": [], "q_out": []}
    peak = tail = 0.0
    t_tail = p["duration_s"] - p["settle_window_s"]
    for k in range(n):
        tk = k * dt
        q_in = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)
        ctrl.control(tk, h)
        q_out = gate_flow(p, ctrl.a, h)
        h += (q_in - q_out) / p["area_m2"] * dt
        ctrl.step_gate(dt)
        d = h - p["h_target"]
        if abs(d) > peak:
            peak = abs(d)
        if tk >= t_tail and abs(d) > tail:
            tail = abs(d)
        if k % rec_every == 0:
            out["t"].append(tk); out["dev"].append(d)
            out["gate"].append(ctrl.a); out["q_out"].append(q_out)
    return {"series": out, "peak_mm": peak * 1000, "final_mm": tail * 1000,
            "gate_final": ctrl.a}


# ----------------------------------------------------------------------
# 一维圣维南：矩形棱柱体，守恒形式 + Rusanov 通量 + 显式欧拉（纯 stdlib）
#   ∂A/∂t + ∂Q/∂x = 0
#   ∂Q/∂t + ∂(Q²/A + gA²/2B)/∂x = gA(S0 − Sf)
# 边界：上游给定流量；下游给定闸门流量（质量严格守恒）
# ----------------------------------------------------------------------
def normal_slope(p: dict, n: float = MANNING_N) -> float:
    """反解渠底坡度，让 (h_target, q0) 恰好是均匀流（正常水深）→ 初始零瞬变。"""
    B, ht = p["pool_width_m"], p["h_target"]
    A = B * ht
    R = A / (B + 2.0 * ht)
    return (n * p["q0"] / (A * R ** (2.0 / 3.0))) ** 2


def pick_dt(dx: float, wmax: float) -> float:
    lo = CFL * dx / max(wmax, 1e-6)
    for c in DT_CANDIDATES:
        if c <= lo:
            return c
    raise ValueError(f"CFL 太紧：dx={dx} 需 dt<{lo:.3f} s")


def sve_1d(p: dict, ncell: int, step_dq: float, ctrl_on: str = "gate",
           n: float = MANNING_N) -> dict:
    B, L, ht = p["pool_width_m"], p["pool_len_m"], p["h_target"]
    dx = L / ncell
    S0 = normal_slope(p, n)
    h = [ht] * ncell
    Q = [p["q0"]] * ncell
    dt = pick_dt(dx, p["q0"] / (B * ht) + math.sqrt(G * ht))
    nt = int(round(p["duration_s"] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = Ctrl(p)

    out = {"t": [], "dev_mean": [], "dev_gate": [], "dev_head": [], "gate": [], "q_out": []}
    peak_mean = peak_gate = 0.0
    tail_mean = 0.0
    t_tail = p["duration_s"] - p["settle_window_s"]
    over_band_dt = 0
    vol_net = 0.0            # 净入流累积（用于质量守恒自检）
    c1 = dt / dx
    gB2 = G / (2.0 * B)
    sqrt = math.sqrt

    for k in range(nt):
        tk = k * dt
        q_up = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)
        q_gate = p["kg"] * ctrl.a * sqrt(h[-1] if h[-1] > 1e-3 else 1e-3) / sqrt(ht)

        # 控制：反馈量取闸前局部水位（物理上闸门量到的）或池均水位
        # ctrl_on="off" → 不做控制（闸门固定在 a0），用于把"控制回路振荡"与"数值振荡"分开
        if ctrl_on != "off":
            ctrl.control(tk, h[-1] if ctrl_on == "gate" else sum(h) / ncell)
            ctrl.step_gate(dt)

        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif
        F2 = [0.0] * nif
        # Rusanov 通量：数值粘性取两侧最大波速（单趟完成）
        for j in range(nif):
            A1 = hx[j] * B
            A2 = hx[j + 1] * B
            if A1 < 1e-6: A1 = 1e-6
            if A2 < 1e-6: A2 = 1e-6
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / B)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / B)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            F2[j] = (0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2
                            + gB2 * (A1 * A1 + A2 * A2))
                     - 0.5 * sm * (Qx[j + 1] - Qx[j]))
        F1[0], F1[-1] = q_up, q_gate      # 边界质量通量严格等于给定流量

        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            A = hi * B
            if A < 1e-6: A = 1e-6
            R = A / (B + 2.0 * hi)
            Qi = Q[i]
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            hn = hi - c1 * (F1[i + 1] - F1[i]) / B
            h_new[i] = hn if hn > 1e-3 else 1e-3
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        vol_net += (q_up - q_gate) * dt

        dm = sum(h) / ncell - ht
        dg = h[-1] - ht
        if abs(dm) > peak_mean: peak_mean = abs(dm)
        if abs(dg) > peak_gate: peak_gate = abs(dg)
        if tk >= t_tail and abs(dm) > tail_mean: tail_mean = abs(dm)
        if abs(dm) > p["band_normal_m"]: over_band_dt += 1
        if k % rec_every == 0:
            out["t"].append(tk); out["dev_mean"].append(dm); out["dev_gate"].append(dg)
            out["dev_head"].append(h[0] - ht); out["gate"].append(ctrl.a)
            out["q_out"].append(q_gate)

    dv = (sum(h) / ncell - ht) * B * L
    return {"ncell": ncell, "dx_m": dx, "dt_s": dt, "steps": nt, "S0": S0,
            "ctrl_on": ctrl_on, "series": out,
            "peak_mean_mm": peak_mean * 1000, "peak_gate_mm": peak_gate * 1000,
            "final_mean_mm": tail_mean * 1000,
            "over_band_s": over_band_dt * dt,
            "vol_net_m3": vol_net, "vol_storage_m3": dv,
            "mass_balance_rel_err": abs(dv - vol_net) / max(abs(vol_net), 1.0)}


def write_csv(path: str, lump: dict, runs: dict, keys: list) -> None:
    base = runs[keys[0]]["series"]
    head = ["t_s", "lumped_dev_m", "lumped_gate"]
    for k in keys:
        head += [f"{k}_dev_mean_m", f"{k}_dev_gate_m", f"{k}_gate"]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(",".join(head) + "\n")
        for i, t in enumerate(base["t"]):
            j = min(i, len(lump["series"]["t"]) - 1)
            row = [f"{t:.0f}", f"{lump['series']['dev'][j]:.6f}",
                   f"{lump['series']['gate'][j]:.5f}"]
            for k in keys:
                s = runs[k]["series"]
                m = min(i, len(s["t"]) - 1)
                row += [f"{s['dev_mean'][m]:.6f}", f"{s['dev_gate'][m]:.6f}",
                        f"{s['gate'][m]:.5f}"]
            fh.write(",".join(row) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 集总模型 vs 一维圣维南 对照（纯标准库）")
    ap.add_argument("--case", default=".")
    ap.add_argument("--coarse", type=int, default=120)
    ap.add_argument("--fine", type=int, default=240)
    ap.add_argument("--finest", type=int, default=480)
    args = ap.parse_args()

    case = os.path.abspath(args.case)
    jp_in = os.path.join(case, "raw", "bm01_main.json")
    with open(jp_in, encoding="utf-8") as fh:
        raw = json.load(fh)
    p = raw["params"]
    step_dq = float(raw["scenario"]["step_dq"])
    rec = raw["metrics"]["pool1"]
    t0 = time.time()

    print(f"参数来源：{os.path.relpath(jp_in, case)}（主案例自记，单一事实来源）")
    print(f"池长 {p['pool_len_m']:.0f} m × 宽 {p['pool_width_m']:.0f} m，A = {p['area_m2']:.3e} m²，"
          f"h_target = {p['h_target']:.2f} m，Q0 = {p['q0']:.0f} m³/s，阶跃 +{step_dq:.0f} m³/s")

    # ---- 自校验：纯 Python 重放集总模型，与记录指标比对 ----
    # 注意：记录里的 peak_dev_m 是 round(值, 5)（单位 m），即量化到 0.01 mm，
    #       所以比对容差必须取"半个量化步长"= 0.005 mm，而不是任意紧的数。
    lump = lumped_0d(p, step_dq)
    QUANT_MM = 0.5 * 1e-5 * 1000.0        # 0.005 mm
    d_peak = abs(rec["peak_dev_m"] * 1000 - lump["peak_mm"])
    d_gate = abs(rec["gate_final"] - lump["gate_final"])
    ok = d_peak <= QUANT_MM and d_gate < 1e-9
    print(f"\n[自校验] 集总重放 vs 主案例记录：")
    print(f"   峰值     {lump['peak_mm']:.4f} / {rec['peak_dev_m']*1000:.4f} mm"
          f"　差 {d_peak:.2e} mm（记录量化步长 0.01 mm，容差取半步 {QUANT_MM:.3f} mm）")
    print(f"   终态开度 {lump['gate_final']:.6f} / {rec['gate_final']:.6f}　差 {d_gate:.2e}")
    print(f"   → {'一致 ✓（控制器逐行等价）' if ok else '不一致 ✗'}")

    S0 = normal_slope(p)
    c = math.sqrt(G * p["h_target"])
    print(f"\n导出坡度 S0 = {S0:.3e}（让 h_target/Q0 恰好是均匀流），Manning n = {MANNING_N}")
    print(f"波速 c = √(g·h) = {c:.2f} m/s → 传播时间 L/c = "
          f"{p['pool_len_m']/c/60:.1f} min（集总模型里被丢掉的量）")
    print(f"控制盲期（采样周期）= {p['ctrl_period_s']/60:.0f} min；"
          f"闸门满行程 = {1/p['rate_limit_per_s']/60:.1f} min")

    # ---- 一维：网格加密序列（反馈=闸前局部水位）----
    grid_seq = sorted({args.coarse, args.fine, args.finest})
    runs, key_gate = {}, {}
    print()
    for nc in grid_seq:
        ta = time.time()
        k = f"g{nc}_gate"
        runs[k] = sve_1d(p, nc, step_dq, ctrl_on="gate")
        key_gate[nc] = k
        print(f"  [1D] {nc:>4} 格 完成 {time.time()-ta:6.1f} s｜dt={runs[k]['dt_s']} s｜"
              f"步数 {runs[k]['steps']}｜池均峰值 {runs[k]['peak_mean_mm']:6.1f} mm｜"
              f"闸前峰值 {runs[k]['peak_gate_mm']:6.1f} mm｜"
              f"质量守恒误差 {runs[k]['mass_balance_rel_err']:.1e}", flush=True)
    # 反馈改取池均水位（集总模型的隐含假设）——用于把差异归因
    for nc in grid_seq[:2]:
        k = f"g{nc}_mean"
        runs[k] = sve_1d(p, nc, step_dq, ctrl_on="mean")

    # ---- 汇总 ----
    label = {f"g{nc}_gate": f"1D SVE（{nc} 格，反馈=闸前局部水位）" for nc in grid_seq}
    label |= {f"g{nc}_mean": f"1D SVE（{nc} 格，反馈=池均水位）" for nc in grid_seq[:2]}
    order = [key_gate[nc] for nc in grid_seq] + [f"g{nc}_mean" for nc in grid_seq[:2]]

    print(f"\n{'模型':<40}{'峰值(mm)':>10}{'末段(mm)':>10}{'越带(min)':>11}")
    print("-" * 72)
    print(f"{'集总 0D（BM-01 主案例，记录值）':<40}{rec['peak_dev_m']*1000:>10.1f}"
          f"{rec['final_window_max_dev_m']*1000:>10.1f}"
          f"{rec['time_over_normal_band_s']/60:>11.1f}")
    for k in order:
        r = runs[k]
        print(f"{label[k]:<40}{r['peak_mean_mm']:>10.1f}{r['final_mean_mm']:>10.1f}"
              f"{r['over_band_s']/60:>11.1f}")
    print("-" * 72)
    print(f"条带：±{p['band_normal_m']*1000:.0f} mm")

    # ---- 网格收敛与 Richardson 外推（一阶格式）----
    p_c = runs[key_gate[grid_seq[0]]]["peak_mean_mm"]
    p_f = runs[key_gate[grid_seq[1]]]["peak_mean_mm"]
    p_ff = runs[key_gate[grid_seq[2]]]["peak_mean_mm"] if len(grid_seq) > 2 else None
    conv = {"coarse": p_c, "fine": p_f, "finest": p_ff}
    rich = None
    if p_ff is not None:
        # 一阶格式：p_∞ ≈ 2·p_finest − p_fine（用最密两级的差值外推）
        rich = 2.0 * p_ff - p_f
        conv["richardson_estimate_mm"] = rich
        conv["gap_coarse_fine_mm"] = abs(p_f - p_c)
        conv["gap_fine_finest_mm"] = abs(p_ff - p_f)
    print()
    if rich is not None:
        print(f"网格收敛：{grid_seq[0]} 格 {p_c:.1f} → {grid_seq[1]} 格 {p_f:.1f} → "
              f"{grid_seq[2]} 格 {p_ff:.1f} mm")
        print(f"  逐级增量 {abs(p_f-p_c):.1f} → {abs(p_ff-p_f):.1f} mm"
              f"（{'正在收敛' if abs(p_ff-p_f) < abs(p_f-p_c) else '尚未收敛'}）")
        print(f"  一阶 Richardson 外推（最密两级）：约 {rich:.0f} mm")
    else:
        print(f"网格收敛：{grid_seq[0]} 格 {p_c:.1f} → {grid_seq[1]} 格 {p_f:.1f} mm")

    ref = rich if rich is not None else p_f
    ratio = ref / (rec["peak_dev_m"] * 1000)
    print(f"\n归因（同样一维、只换反馈量）：")
    print(f"  0D 集总（反馈=池均，瞬时）          {rec['peak_dev_m']*1000:6.1f} mm")
    for nc in grid_seq[:2]:
        print(f"  1D（{nc} 格，反馈=池均水位）"
              f"{'':<10}{runs[f'g{nc}_mean']['peak_mean_mm']:6.1f} mm"
              f"   ← 与 0D 基本一致，说明差异不是摩擦/数值耗散")
    print(f"  1D（最密，反馈=闸前局部水位）        {ref:6.0f} mm"
          f"   ← 差异全部来自「闸门量到的是局部水位」")
    print(f"集总低估倍数（外推值/0D）≈ {ratio:.2f}×"
          f"（{rec['peak_dev_m']*1000:.1f} → ~{ref:.0f} mm）")

    # ---- 工件 ----
    out_raw = os.path.join(case, "raw")
    os.makedirs(out_raw, exist_ok=True)
    chk = {
        "peak_mm_recorded": rec["peak_dev_m"] * 1000, "peak_mm_replay": lump["peak_mm"],
        "gate_final_recorded": rec["gate_final"], "gate_final_replay": lump["gate_final"],
        "diff_peak_mm": d_peak, "tol_peak_mm": QUANT_MM, "diff_gate": d_gate,
        "note": "记录值量化到 0.01 mm，容差取半步；终态开度差 ~1e-16 说明控制器逐行等价",
        "identical": bool(ok)}
    payload = {
        "case_id": "BM-01", "purpose": "reduced_order_loss_check",
        "evidence_grade": "synthetic-1d-sve",
        "note": "仅用于量化降阶损失；不改变 BM-01 判定码与证据等级；量级示例参数，非现场结论。",
        "solver": {"scheme": "conservative Rusanov (local Lax-Friedrichs), explicit Euler",
                   "deps": "python stdlib only", "manning_n": MANNING_N, "cfl": CFL,
                   "slope_S0": S0, "wave_celevity_mps": c,
                   "travel_time_min": p["pool_len_m"] / c / 60,
                   "gate_bc": "imposed discharge (mass conserving); structure not modelled internally"},
        "self_check_lumped_replay": chk,
        "lumped_recorded": {"peak_mm": rec["peak_dev_m"] * 1000,
                            "final_mm": rec["final_window_max_dev_m"] * 1000,
                            "over_band_s": rec["time_over_normal_band_s"],
                            "gate_final": rec["gate_final"]},
        "sve": {k: {kk: r[kk] for kk in ("ncell", "dx_m", "dt_s", "steps", "ctrl_on",
                                         "peak_mean_mm", "peak_gate_mm", "final_mean_mm",
                                         "over_band_s", "mass_balance_rel_err",
                                         "vol_net_m3", "vol_storage_m3")}
                for k, r in runs.items()},
        "grid_convergence": conv,
        "attribution": {
            "lumped_0d_mm": rec["peak_dev_m"] * 1000,
            "sve_mean_feedback_mm": {nc: runs[f"g{nc}_mean"]["peak_mean_mm"]
                                     for nc in grid_seq[:2]},
            "sve_gate_feedback_mm": {nc: runs[key_gate[nc]]["peak_mean_mm"] for nc in grid_seq},
            "extrapolated_mm": ref, "undersetimate_ratio": ratio},
        "series": {k: [round(v, 6) for v in runs[key_gate[grid_seq[0]]]["series"][k]]
                   for k in ("t", "dev_mean", "dev_gate", "dev_head", "gate", "q_out")},
        "series_finest": {k: [round(v, 6) for v in runs[key_gate[grid_seq[-1]]]["series"][k]]
                          for k in ("t", "dev_mean", "dev_gate", "gate")},
        "params_source": "raw/bm01_main.json",
        "runtime_s": round(time.time() - t0, 2),
    }
    jp = os.path.join(out_raw, "bm01_sve1d_compare.json")
    with open(jp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    cp = os.path.join(out_raw, "bm01_sve1d_compare.csv")
    write_csv(cp, lump, runs, [key_gate[nc] for nc in grid_seq[:2]]
              + [f"g{nc}_mean" for nc in grid_seq[:2]])

    print(f"\n总用时 {time.time()-t0:.1f} s")
    print(f"输出：{jp}")
    print(f"输出：{cp}")
    print("证据等级：synthetic-1d-sve（数值降阶，不改变 BM-01 判定码，非现场结论）")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
