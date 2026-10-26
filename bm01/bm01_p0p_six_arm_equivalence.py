# -*- coding: utf-8 -*-
"""
任务 六臂（双模型共识第 2 项 · 升级版）：0D ↔ 1D-全池平均 等价性 —— 矩形 vs 梯形双断面
================================================================================
把 P0-1'（缺口A）从"单断面"扩到"双断面"，检验等价性是否**断面稳健**。

六臂（同一 PI 控制器 kp=1.6/ki=0.0016，12 h，+45 m³/s，480 格）：
  A0-m0   0D 蓄量模型（矩形，A=L·b 常数）
  A1-m0   1D-mean（矩形）
  A0-m15  0D 蓄量模型（梯形 m=1.5，**变面积** dh/dt=(Qin−Qout)/(L·B(h))，B(h)=b+2mh）
  A1-m15  1D-mean（梯形 m=1.5）
  A1g-m15 1D 闸前单点（梯形）
  A1a-m15 1D 闸前 100m 平均（梯形）

闭合性检验两式：
  (i)  一阶：d(h̄)/dt  vs  (Qin−Qout)/(L·B(h̄))      —— 矩形下为恒等式；梯形下含二阶项
  (ii) 二阶项估计：mean(A(h)) − A(h̄) ≈ m·(h² 的空间方差)，用 h 的沿程方差近似（本文未存全场，
       故以 (i) 的残差作为二阶项的上界估计，并在结论中如实标注该近似）。

输出：raw/bm01_p0p_six_arm_equivalence.json
"""
import json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2
import bm01_pool_step as PS

DURATION_S = 43200.0
DT_0D = 5.0
NCELL = 480
STEP_DQ = 45.0
BAND_MM = 150.0
SETTLE_WINDOW_S = 1800.0
B0, L, HT = 50.0, 6000.0, 7.0

P1D = {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": HT, "a0": 0.5, "kg": 360.0,
       "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
       "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
       "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
       "settle_window_s": SETTLE_WINDOW_S, "duration_s": DURATION_S}


def verdict(peak_mm):
    return "A1(带内)" if peak_mm <= BAND_MM else "A3(越带)"


def run_0d_rect():
    """既有 0D 蓄量模型（矩形，A 常数）。"""
    p = PS.Params(duration_s=DURATION_S, dt=DT_0D, t_step_s=900.0, step_dq=STEP_DQ)
    sc = PS.Scenario(key="main", label="", pools=1, step_dq=STEP_DQ, feedforward=False, note="")
    r = PS.simulate(p, sc)
    t = r["series"]["t"]; dev = r["series"]["h1_dev"]; a = r["series"]["a1"]; qo = r["series"]["q1_out"]
    rec = max(1, int(round(30.0 / DT_0D)))
    n = len(t)
    peak = float(max(abs(x) for x in dev)) * 1000.0
    i_tail = int(round((DURATION_S - SETTLE_WINDOW_S) / DT_0D))
    tail = float(max(abs(x) for x in dev[i_tail:])) * 1000.0
    return {"t": [float(t[k]) for k in range(0, n, rec)],
            "dev_mm": [float(dev[k]) * 1000.0 for k in range(0, n, rec)],
            "gate": [float(a[k]) for k in range(0, n, rec)],
            "q_out": [float(qo[k]) for k in range(0, n, rec)],
            "peak_mm": round(peak, 2), "final_mm": round(tail, 2), "judgement": verdict(peak)}


def run_0d_trap(m):
    """梯形变面积 0D：dV/dt = Qin − Qout，V = L·A(h)，A(h)=h(b+mh)。
       → dh/dt = (Qin − Qout) / (L·B(h))，B(h)=b+2mh。闸门率定用同一 h（均一水位假设）。"""
    sec = V2.Section(B0, m)
    S0 = V2.section_normal_slope(sec, HT, P1D["q0"], V2.MANNING_N)
    h = HT; a = 0.5; a_cmd = 0.5; integral = 0.0; last_t = -1e18
    dt = DT_0D
    ts, dev, gate, qout = [], [], [], []
    k = 0
    while True:
        tk = k * dt
        q_up = P1D["q0"] + (STEP_DQ if tk >= 900.0 else 0.0)
        q_g = 360.0 * a * math.sqrt(max(h, 1e-3)) / math.sqrt(HT)
        # 控制器（与 1D Ctrl 完全一致）
        if tk - last_t >= 300.0:
            last_t = tk
            meas = round(h / 0.001) * 0.001
            e = meas - HT
            at_hi = a >= 1.0 - 1e-12 and e > 0
            at_lo = a <= 1e-12 and e < 0
            if abs(e) >= 0.01 and not (at_hi or at_lo):
                integral += e * 300.0
            a_cmd = min(max(0.5 + 1.6 * e + 0.0016 * integral, 0.0), 1.0)
        # 速率限幅
        step = (1.0 / 3000.0) * dt
        d = a_cmd - a
        a = a + d if abs(d) < step else a + (step if d > 0 else -step)
        a = min(max(a, 0.0), 1.0)
        # 蓄量推进（变面积）
        area_surf = L * (B0 + 2.0 * m * h)
        h += (q_up - q_g) / area_surf * dt
        if k % max(1, int(round(30.0 / dt))) == 0:
            ts.append(tk); dev.append(h - HT); gate.append(a); qout.append(q_g)
        k += 1
        if tk >= DURATION_S: break
    peak = max(abs(x) for x in dev) * 1000.0
    i_tail = int(round((DURATION_S - SETTLE_WINDOW_S) / 30.0))
    tail = max(abs(x) for x in dev[i_tail:]) * 1000.0
    return {"t": ts, "dev_mm": [x * 1000.0 for x in dev], "gate": gate, "q_out": qout,
            "peak_mm": round(peak, 2), "final_mm": round(tail, 2), "judgement": verdict(peak),
            "S0": S0, "area_surf_m2_initial": L * B0}


def run_1d(mode, m, sensor_len=100.0):
    kw = {}
    if mode == "gate_avg": kw["sensor_len_m"] = sensor_len
    s0 = None if m == 0.0 else V2.section_normal_slope(V2.Section(B0, m), HT, P1D["q0"], V2.MANNING_N)
    r = V2.sve_1d_v2(P1D, NCELL, STEP_DQ, ctrl_on=mode, m_slope=m, s0=s0, **kw)
    s = r["series"]
    return {"t": [float(x) for x in s["t"]],
            "dev_mean_mm": [float(x) * 1000.0 for x in s["dev_mean"]],
            "dev_gate_mm": [float(x) * 1000.0 for x in s["dev_gate"]],
            "gate": [float(x) for x in s["gate"]],
            "q_out": [float(x) for x in s["q_out"]],
            "peak_mean_mm": round(r["peak_mean_mm"], 2),
            "final_mean_mm": round(r["final_mean_mm"], 2),
            "peak_gate_mm": round(r["peak_gate_mm"], 2),
            "mass_balance_rel_err": r["mass_balance_rel_err"],
            "judgement": verdict(r["peak_mean_mm"])}


def interp(xs, ys, x):
    if x <= xs[0]: return ys[0]
    if x >= xs[-1]: return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        m_ = (lo + hi) // 2
        if xs[m_] <= x: lo = m_
        else: hi = m_
    f = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + f * (ys[hi] - ys[lo])


def compare(a, b, ta, tb):
    ya = [interp(ta, a, x) for x in tb]
    n = len(ya); d = [ya[i] - b[i] for i in range(n)]
    rms = math.sqrt(sum(x * x for x in d) / n); mx = max(abs(x) for x in d)
    ma, mb = sum(ya) / n, sum(b) / n
    cov = sum((ya[i] - ma) * (b[i] - mb) for i in range(n)) / n
    va = math.sqrt(sum((x - ma) ** 2 for x in ya) / n); vb = math.sqrt(sum((x - mb) ** 2 for x in b) / n)
    return {"rms_mm": round(rms, 3), "max_abs_mm": round(mx, 3),
            "corr": round(cov / (va * vb), 6) if va * vb > 1e-12 else None, "n": n}


def closure(t, hdev_mm, qout, m):
    """d(h̄)/dt vs (Qin−Qout)/(L·B(h̄))；返回 mm/min 残差统计。"""
    res = []
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        if dt <= 0: continue
        hbar = HT + 0.5 * (hdev_mm[i] + hdev_mm[i - 1]) / 1000.0
        area = L * (B0 + 2.0 * m * hbar)
        tmid = 0.5 * (t[i] + t[i - 1])
        q_up = P1D["q0"] + (STEP_DQ if tmid >= 900.0 else 0.0)
        q_g = 0.5 * (qout[i] + qout[i - 1])
        closure_rate = (q_up - q_g) / area
        dHdt = (hdev_mm[i] - hdev_mm[i - 1]) / 1000.0 / dt
        res.append((dHdt - closure_rate) * 1000.0 * 60.0)
    return {"rms_mm_per_min": round(math.sqrt(sum(x * x for x in res) / len(res)), 6),
            "max_abs_mm_per_min": round(max(abs(x) for x in res), 6)}


def main():
    out = {"purpose": "六臂：0D↔1D全池平均等价性（矩形 vs 梯形）",
           "setup": {"duration_h": 12.0, "step_dq": STEP_DQ, "ncell": NCELL,
                     "band_mm": BAND_MM, "kp": 1.6, "ki": 0.0016,
                     "prereg": "H0: |Δpeak|/peak<=15% AND 判定码一致 AND RMSE<=20mm"},
           "arms": {}, "comparisons": {}}
    print("[六臂] 0D 矩形 / 0D 梯形…", flush=True)
    a0m0 = run_0d_rect(); a0m15 = run_0d_trap(1.5)
    out["arms"]["A0-m0"] = {k: v for k, v in a0m0.items() if k not in ("t", "dev_mm", "gate", "q_out")}
    out["arms"]["A0-m15"] = {k: v for k, v in a0m15.items() if k not in ("t", "dev_mm", "gate", "q_out")}
    print(f"    A0-m0 peak={a0m0['peak_mm']} {a0m0['judgement']} | A0-m15 peak={a0m15['peak_mm']} {a0m15['judgement']}", flush=True)

    for key, mode, m, sl in (("A1-m0", "mean", 0.0, 0.0), ("A1-m15", "mean", 1.5, 0.0),
                             ("A1g-m15", "gate", 1.5, 0.0), ("A1a-m15", "gate_avg", 1.5, 100.0)):
        print(f"[六臂] 跑 1D {key}…", flush=True)
        r = run_1d(mode, m, sl)
        out["arms"][key] = {k: v for k, v in r.items()
                            if k not in ("t", "dev_mean_mm", "dev_gate_mm", "gate", "q_out")}
        print(f"    {key}: peak={r['peak_mean_mm']} final={r['final_mean_mm']} → {r['judgement']}", flush=True)
        if key == "A1-m0":
            out["comparisons"]["m0_0D_vs_1Dmean"] = {
                "traj": compare(a0m0["dev_mm"], r["dev_mean_mm"], a0m0["t"], r["t"]),
                "peak_0d": a0m0["peak_mm"], "peak_1d": r["peak_mean_mm"],
                "peak_rel_diff_pct": round(abs(r["peak_mean_mm"] - a0m0["peak_mm"]) / a0m0["peak_mm"] * 100, 2),
                "verdict_0d": a0m0["judgement"], "verdict_1d": r["judgement"],
                "verdict_match": a0m0["judgement"] == r["judgement"],
                "closure": closure(r["t"], r["dev_mean_mm"], r["q_out"], 0.0)}
        if key == "A1-m15":
            out["comparisons"]["m15_0D_vs_1Dmean"] = {
                "traj": compare(a0m15["dev_mm"], r["dev_mean_mm"], a0m15["t"], r["t"]),
                "peak_0d": a0m15["peak_mm"], "peak_1d": r["peak_mean_mm"],
                "peak_rel_diff_pct": round(abs(r["peak_mean_mm"] - a0m15["peak_mm"]) / max(a0m15["peak_mm"], 1e-9) * 100, 2),
                "verdict_0d": a0m15["judgement"], "verdict_1d": r["judgement"],
                "verdict_match": a0m15["judgement"] == r["judgement"],
                "closure": closure(r["t"], r["dev_mean_mm"], r["q_out"], 1.5)}

    out["cross_section_summary"] = {
        "m0_peak_rel_diff_pct": out["comparisons"]["m0_0D_vs_1Dmean"]["peak_rel_diff_pct"],
        "m0_rmse_mm": out["comparisons"]["m0_0D_vs_1Dmean"]["traj"]["rms_mm"],
        "m15_peak_rel_diff_pct": out["comparisons"]["m15_0D_vs_1Dmean"]["peak_rel_diff_pct"],
        "m15_rmse_mm": out["comparisons"]["m15_0D_vs_1Dmean"]["traj"]["rms_mm"],
        "verdicts": {k: out["arms"][k]["judgement"] for k in
                     ("A0-m0", "A1-m0", "A0-m15", "A1-m15", "A1g-m15", "A1a-m15")},
    }
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_p0p_six_arm_equivalence.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n已写入 {jp}")
    s = out["cross_section_summary"]
    print(f"矩形：Δpeak {s['m0_peak_rel_diff_pct']}% RMSE {s['m0_rmse_mm']}mm ｜ 梯形：Δpeak {s['m15_peak_rel_diff_pct']}% RMSE {s['m15_rmse_mm']}mm")
    print(f"判定码：{s['verdicts']}")


if __name__ == "__main__":
    main()
