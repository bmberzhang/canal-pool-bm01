# -*- coding: utf-8 -*-
"""
P0-1' 缺口A 闭合实验：0D 蓄量模型 与「1D 解的全池瞬时空间平均」是否等价？
================================================================================
论文核心叙事是「0D 的全池均一水位假设 ↔ 对 1D 解做全池空间平均」。
此前只有「峰值相近（72.7 vs 74.8 mm）」这一旁证，从未做过
**同对象 / 同控制器 / 同时长 / 逐时刻** 的对照实验。

本脚本在同一物理对象（BM-01 单池，L=6000 m，B=50 m）、同一 PI 控制器
（kp=1.6, ki=0.0016，含死区 / 条件积分抗饱和 / 速率限幅）、同一 12 h 时长、
同一 +45 m³/s 上游阶跃下跑三条臂：

  arm-0D        : bm01_pool_step.py 的 0D 集总蓄量模型（反馈 = 池水位 h）
  arm-1D-mean   : 1D 圣维南，ctrl_on="mean"（反馈 = 全池瞬时空间平均 Σh/n）
  arm-1D-100m   : 1D 圣维南，ctrl_on="gate_avg" 闸前 100 m 平均（对照）
  arm-1D-point  : 1D 圣维南，ctrl_on="gate" 闸前单点（对照）

并做两件此前没做的事：
  (1) **闭合性检验**：对 1D-mean 解数值求 d(h̄)/dt，与 0D 的蓄量方程
      (Q_up − Q_gate)/(B·L) 逐时刻比对。矩形棱柱渠池下，0D 方程恰是
      1D 连续方程的全池空间积分，故该残差应 ≈ 机器精度 —— 这直接把
      「0D = 1D 的空间平均」从假设升级为**可验证的恒等式**。
  (2) **归因分解**：0D 与 1D 唯一的物理差别在于闸门出流率定用什么水位
      （0D 用 h̄，1D 的 Q_gate = kg·a·√h[-1] 用闸前单点 h[-1]）。于是
      记录 1D 解的「闸前单点 − 全池平均」偏差 (h[-1] − h̄)，量化非均匀度，
      即 0D 假设的**闭合误差来源**。

预登记判据（先写死，防事后拟合）：
  H0 等价：|peak_1D − peak_0D| / peak_0D ≤ 15%  且  判定码一致  且  反馈序列 RMSE ≤ 20 mm
  H1 不等价：以上任一不满足
输出：raw/bm01_p0p_0d_1d_equivalence.json
"""
import json, math, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2          # 1D 圣维南（纯 stdlib）
import bm01_pool_step as PS         # 0D 蓄量模型（需 numpy）

DURATION_S = 43200.0                # 12 h
DT_0D = 5.0                         # 0D 时间步（与既有 0D 主程序一致）
NCELL = 480                         # 1D 网格
STEP_DQ = 45.0
BAND_MM = 150.0                     # 判定阈值：peak > 150 mm → A3
SETTLE_WINDOW_S = 1800.0

P1D = {
    "pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0,
    "h_target": 7.0, "a0": 0.5, "kg": 360.0, "kp": 1.6, "ki": 0.0016,
    "ctrl_period_s": 300.0, "deadband_m": 0.01, "meas_resolution_m": 0.001,
    "rate_limit_per_s": 1.0 / 3000.0,
    "t_step_s": 900.0, "band_normal_m": 0.15, "band_alarm_m": 0.3,
    "settle_tol_m": 0.02, "settle_window_s": SETTLE_WINDOW_S,
    "duration_s": DURATION_S,
}


def verdict(peak_mm):
    return "A1(带内)" if peak_mm <= BAND_MM else "A3(越带)"


# ---------------- 0D 臂 ----------------
def run_0d():
    t0 = time.time()
    p = PS.Params(duration_s=DURATION_S, dt=DT_0D, t_step_s=900.0, step_dq=STEP_DQ)
    sc = PS.Scenario(key="main", label="0D 单池 + 本地水位 PI",
                     pools=1, step_dq=STEP_DQ, feedforward=False, note="")
    r = PS.simulate(p, sc)
    t = r["series"]["t"]            # numpy, 每 5 s
    dev = r["series"]["h1_dev"]     # m
    a = r["series"]["a1"]
    q_out = r["series"]["q1_out"]
    n = len(t)
    # 峰值 / 末段（口径与 1D 一致：全程峰值、末段窗口内峰值）
    peak = float(max(abs(x) for x in dev)) * 1000.0
    i_tail = int(round((DURATION_S - SETTLE_WINDOW_S) / DT_0D))
    tail = float(max(abs(x) for x in dev[i_tail:])) * 1000.0
    # 下采样到 30 s（与 1D 记录口径一致）
    rec = max(1, int(round(30.0 / DT_0D)))
    return {
        "t": [float(t[k]) for k in range(0, n, rec)],
        "dev_mm": [float(dev[k]) * 1000.0 for k in range(0, n, rec)],
        "gate": [float(a[k]) for k in range(0, n, rec)],
        "q_out": [float(q_out[k]) for k in range(0, n, rec)],
        "peak_mm": round(peak, 2), "final_mm": round(tail, 2),
        "judgement": verdict(peak), "wall_s": round(time.time() - t0, 1),
    }


# ---------------- 1D 臂 ----------------
def run_1d(mode, sensor_len=100.0):
    t0 = time.time()
    kw = {}
    if mode == "gate_avg":
        kw["sensor_len_m"] = sensor_len
    r = V2.sve_1d_v2(P1D, NCELL, STEP_DQ, ctrl_on=mode, m_slope=0.0, **kw)
    s = r["series"]
    return {
        "label": {"mean": "mean(全池)", "gate_avg": f"gate_avg({sensor_len:.0f}m)",
                  "gate": "gate(单点)"}[mode],
        "t": [float(x) for x in s["t"]],
        "dev_mean_mm": [float(x) * 1000.0 for x in s["dev_mean"]],
        "dev_gate_mm": [float(x) * 1000.0 for x in s["dev_gate"]],
        "dev_head_mm": [float(x) * 1000.0 for x in s["dev_head"]],
        "gate": [float(x) for x in s["gate"]],
        "q_out": [float(x) for x in s["q_out"]],
        "peak_mean_mm": round(r["peak_mean_mm"], 2),
        "final_mean_mm": round(r["final_mean_mm"], 2),
        "peak_gate_mm": round(r["peak_gate_mm"], 2),
        "over_band_s": round(r["over_band_s"], 1),
        "mass_balance_rel_err": r["mass_balance_rel_err"],
        "dt_s": r["dt_s"],
        "judgement": verdict(r["peak_mean_mm"]),
        "wall_s": round(time.time() - t0, 1),
    }


# ---------------- 工具 ----------------
def interp(xs, ys, x):
    """线性插值（xs 单调递增）。"""
    if x <= xs[0]: return ys[0]
    if x >= xs[-1]: return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x: lo = mid
        else: hi = mid
    f = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + f * (ys[hi] - ys[lo])


def compare(a, b, ta, tb):
    """把 a（时间轴 ta）插值到 b 的时间轴 tb 上，返回差异统计。"""
    ya = [interp(ta, a, x) for x in tb]
    n = len(ya)
    d = [ya[i] - b[i] for i in range(n)]
    rms = math.sqrt(sum(x * x for x in d) / n)
    mx = max(abs(x) for x in d)
    ma, mb = sum(ya) / n, sum(b) / n
    cov = sum((ya[i] - ma) * (b[i] - mb) for i in range(n)) / n
    va = math.sqrt(sum((x - ma) ** 2 for x in ya) / n)
    vb = math.sqrt(sum((x - mb) ** 2 for x in b) / n)
    corr = cov / (va * vb) if va * vb > 1e-12 else float("nan")
    return {"rms_mm": round(rms, 3), "max_abs_mm": round(mx, 3),
            "corr": round(corr, 6), "n": n}


def main():
    out = {"purpose": "P0-1' 0D vs 1D全池平均 等价性对照",
           "setup": {"duration_h": DURATION_S / 3600.0, "step_dq": STEP_DQ,
                     "band_mm": BAND_MM, "ncell": NCELL, "dt_0d": DT_0D,
                     "kp": P1D["kp"], "ki": P1D["ki"],
                     "prereg_criteria": "H0: |peak1D-peak0D|/peak0D<=15% AND 判定码一致 AND RMSE<=20mm"},
           "arms": {}}

    print("[P0-1'] 跑 0D 臂（12h, dt=5s）…", flush=True)
    r0 = run_0d()
    out["arms"]["0D"] = {k: v for k, v in r0.items() if k not in ("t", "dev_mm", "gate", "q_out")}
    print(f"    0D: peak={r0['peak_mm']} final={r0['final_mm']} -> {r0['judgement']} ({r0['wall_s']}s)", flush=True)
    out["arms"]["0D"]["series_n"] = len(r0["t"])

    for mode, sl in (("mean", 0.0), ("gate_avg", 100.0), ("gate", 0.0)):
        print(f"[P0-1'] 跑 1D 臂 {mode}…", flush=True)
        r1 = run_1d(mode, sl)
        key = r1["label"]
        out["arms"][key] = {k: v for k, v in r1.items()
                            if k not in ("t", "dev_mean_mm", "dev_gate_mm", "dev_head_mm", "gate", "q_out")}
        print(f"    {key}: peak={r1['peak_mean_mm']} final={r1['final_mean_mm']} -> {r1['judgement']} ({r1['wall_s']}s)", flush=True)

        if mode == "mean":
            # ---- (1) 轨迹对照：0D vs 1D-mean ----
            out["trajectory_0d_vs_1d_mean"] = {
                "level_dev": compare(r0["dev_mm"], r1["dev_mean_mm"], r0["t"], r1["t"]),
                "gate_cmd": compare(r0["gate"], r1["gate"], r0["t"], r1["t"]),
                "peak_0d_mm": r0["peak_mm"], "peak_1d_mm": r1["peak_mean_mm"],
                "peak_rel_diff_pct": round(abs(r1["peak_mean_mm"] - r0["peak_mm"]) / max(r0["peak_mm"], 1e-9) * 100.0, 2),
                "final_0d_mm": r0["final_mm"], "final_1d_mm": r1["final_mean_mm"],
                "verdict_0d": r0["judgement"], "verdict_1d": r1["judgement"],
                "verdict_match": r0["judgement"] == r1["judgement"],
            }

            # ---- (2) 闭合性检验：1D 的 d(h̄)/dt vs 0D 蓄量方程 ----
            B, L = P1D["pool_width_m"], P1D["pool_len_m"]
            area = B * L
            t = r1["t"]; hm = r1["dev_mean_mm"]; qo = r1["q_out"]
            resid_mm_per_min = []
            for i in range(1, len(t)):
                dt = t[i] - t[i - 1]
                if dt <= 0: continue
                dH_dt = (hm[i] - hm[i - 1]) / 1000.0 / dt            # m/s
                tmid = 0.5 * (t[i] + t[i - 1])
                q_up = P1D["q0"] + (STEP_DQ if tmid >= P1D["t_step_s"] else 0.0)
                q_g = 0.5 * (qo[i] + qo[i - 1])
                closure = (q_up - q_g) / area                        # 0D 预报的 dh̄/dt
                resid_mm_per_min.append((dH_dt - closure) * 1000.0 * 60.0)
            rms_res = math.sqrt(sum(x * x for x in resid_mm_per_min) / len(resid_mm_per_min))
            idx = sorted(range(len(resid_mm_per_min)), key=lambda i: -abs(resid_mm_per_min[i]))[:5]
            out["closure_check_1d_mean"] = {
                "note": "d(全池平均水位)/dt 与 0D 蓄量方程 (Q_up-Q_gate)/(B*L) 之差；矩形棱柱池应为机器精度",
                "rms_mm_per_min": round(rms_res, 6),
                "max_abs_mm_per_min": round(max(abs(x) for x in resid_mm_per_min), 6),
                "top5_residuals_mm_per_min": [round(resid_mm_per_min[i], 6) for i in idx],
                "mass_balance_rel_err": r1["mass_balance_rel_err"],
            }

            # ---- (3) 非均匀度：1D 的「闸前单点 − 全池平均」= 0D 假设的误差来源 ----
            diff = [r1["dev_gate_mm"][i] - r1["dev_mean_mm"][i] for i in range(len(t))]
            out["nonuniformity_1d"] = {
                "note": "Q_gate 用 h[-1]，0D 用 h̄；两者之差即 0D 闸门率定的闭合误差来源",
                "peak_single_minus_mean_mm": round(max(diff), 2),
                "min_single_minus_mean_mm": round(min(diff), 2),
                "peak_mean_mm": r1["peak_mean_mm"],
                "peak_gate_mm": r1["peak_gate_mm"],
                "gate_over_mean_ratio": round(r1["peak_gate_mm"] / max(r1["peak_mean_mm"], 1e-9), 4),
            }

    # ---- (4) 反事实合成：把 1D-mean 的水位序列喂给 0D 的闸门率定 ----
    # 若 0D 用 h̄ 而闸门率定也用 h̄，则判定码如何？（用 1D-mean 结果直接回答）
    out["counterfactual"] = {
        "question": "0D 与 1D 判定码差异是否完全来自「闸门率定取样位置」(h[-1] vs h̄)？",
        "peak_1D_feedback_mean_gate_local_mm": out["arms"]["mean(全池)"]["peak_mean_mm"],
        "peak_1D_feedback_mean_gate_local_verdict": out["arms"]["mean(全池)"]["judgement"],
        "peak_0D_mm": r0["peak_mm"],
        "peak_0D_verdict": r0["judgement"],
        "peak_1D_point_mm": out["arms"]["gate(单点)"]["peak_mean_mm"],
        "peak_1D_point_verdict": out["arms"]["gate(单点)"]["judgement"],
        "peak_1D_100m_mm": out["arms"]["gate_avg(100m)"]["peak_mean_mm"],
        "peak_1D_100m_verdict": out["arms"]["gate_avg(100m)"]["judgement"],
    }

    os.makedirs("raw", exist_ok=True)
    path = os.path.join("raw", "bm01_p0p_0d_1d_equivalence.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\n已写入 {path}")
    tc = out["trajectory_0d_vs_1d_mean"]
    print(f"结论：0D vs 1D-mean 峰值差 {tc['peak_rel_diff_pct']}% | RMSE {tc['level_dev']['rms_mm']} mm | "
          f"corr {tc['level_dev']['corr']} | 判定码一致={tc['verdict_match']}")


if __name__ == "__main__":
    main()
