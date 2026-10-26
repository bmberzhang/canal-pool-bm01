# -*- coding: utf-8 -*-
"""
(W0) C3 余量谱：把 C3 从"一个正类"做成"有分离度的判据"
================================================================================
三模型第三轮一致指出：C3 目前只有 **1 个正类**（L=2000/ΔQ=45，ΔQ/ΔQ* = 1.062，余量仅 6.2%），
"单点分离不构成判据"。要求先出**余量谱**：
  - 8 格 T_band/T_slew（以及等价量 ΔQ/ΔQ*）全列；
  - **密扫边界**：在 L=2000 上把 ΔQ 从 20 扫到 60（步长 2），找出"三档全 A3"的真实转变点；
  - 比较真实转变点与解析阈值 ΔQ* = √(k_g·r·A_surf·band)；
  - 若正类与最近负类的间距 < 10%，正文只能写"诊断不等式"，不得写"判定阈值"。

设置：固定物理网格 dx = 25 m（与 E-MC 一致），12 h，矩形，三档取样。
输出：raw/bm01_pred_spectrum.json
"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

DX_M = 25.0
L = 2000.0
NCELL = int(round(L / DX_M))
BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
DQS = list(range(8, 42, 2))           # 8,10,...,40（覆盖真实转变区 (10,20]，见首轮扫描）
ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def base_p():
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def main():
    t0 = time.time()
    A_surf = B0 * L
    dq_star = math.sqrt(KG * RATE * A_surf * 0.15)
    print("解析阈值 ΔQ* = %.3f m³/s（A_surf=%.0f m², 条带 150 mm）" % (dq_star, A_surf), flush=True)

    rows = []
    for dq in DQS:
        arms = {}
        for name, kw in ARMS:
            r = V2.sve_1d_v2(base_p(), NCELL, float(dq), m_slope=0.0, **kw)
            arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                          "final_mean_mm": round(r["final_mean_mm"], 2)}
        n_a3 = sum(1 for a in ARMS if arms[a[0]]["peak_mean_mm"] > BAND_MM)
        T_band = A_surf * 0.15 / dq
        T_slew = dq / (KG * RATE)
        rows.append({"dq": dq, "peaks_mm": {k: v["peak_mean_mm"] for k, v in arms.items()},
                     "n_A3_arms": n_a3, "ALL_A3": n_a3 == 3,
                     "T_band_s": round(T_band, 1), "T_slew_s": round(T_slew, 1),
                     "T_band_over_T_slew": round(T_band / T_slew, 4),
                     "dq_over_dqstar": round(dq / dq_star, 4)})
        print("[SPEC] L=%.0f ΔQ=%4.1f | point=%8.1f avg100=%8.1f mean=%8.1f | A3档=%d %s | "
              "T_band/T_slew=%.3f  ΔQ/ΔQ*=%.3f"
              % (L, dq, arms["point"]["peak_mean_mm"], arms["avg100"]["peak_mean_mm"],
                 arms["mean"]["peak_mean_mm"], n_a3, "ALL_A3" if n_a3 == 3 else "     ",
                 T_band / T_slew, dq / dq_star), flush=True)

    all_a3 = [r["dq"] for r in rows if r["ALL_A3"]]
    not_a3 = [r["dq"] for r in rows if not r["ALL_A3"]]
    onset = min(all_a3) if all_a3 else None
    last_ok = max(not_a3) if not_a3 else None

    # 分离度：以 ΔQ/ΔQ* 度量
    pos = [r["dq_over_dqstar"] for r in rows if r["ALL_A3"]]
    neg = [r["dq_over_dqstar"] for r in rows if not r["ALL_A3"]]
    sep = None
    if pos and neg:
        sep = {"pos_min": round(min(pos), 4), "neg_max": round(max(neg), 4),
               "gap_abs": round(min(pos) - max(neg), 4),
               "gap_rel_pct": round((min(pos) / max(neg) - 1.0) * 100.0, 2)}

    out = {"purpose": "W0 C3 余量谱：密扫边界，检验 C3 是判据还是仅是诊断不等式",
           "config": {"L_m": L, "ncell": NCELL, "dx_m": DX_M, "duration_h": 12.0, "band_mm": BAND_MM,
                      "dqs": DQS, "A_surf_m2": A_surf, "dq_star_m3s": round(dq_star, 3)},
           "rows": rows, "all_A3_dqs": all_a3, "not_all_A3_dqs": not_a3,
           "onset_dq_all_A3": onset, "last_dq_not_all_A3": last_ok,
           "separation_in_dq_over_dqstar": sep,
           "wall_total_s": round(time.time() - t0, 1)}

    if onset is not None and last_ok is not None:
        pct_err = abs(onset - dq_star) / dq_star * 100.0
        wide = sep and sep["gap_rel_pct"] >= 10.0
        out["verdict"] = (
            "实测转变点落在 ΔQ ∈ (%d, %d]，即 %.1f m³/s（解析 ΔQ*=%.2f，"
            "偏差 %.1f%%）。正类 ΔQ/ΔQ* 最小值 %.3f、最近负类 %.3f，相对间距 %.1f%% → %s。"
            % (last_ok, onset, (last_ok + onset) / 2.0, dq_star, pct_err,
               sep["pos_min"], sep["neg_max"], sep["gap_rel_pct"],
               "间距 ≥10%，C3 可作为**判定阈值**写入正文" if wide
               else "间距 <10%，C3 只能写为**诊断不等式**（与观测相容），不得写为硬阈值"))
    else:
        out["verdict"] = "扫描区间内未出现（或全程出现）全 A3，无法定位转变点，需扩大扫描区间。"

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_pred_spectrum.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n全 A3 的 ΔQ:", all_a3)
    print("非全 A3 的 ΔQ:", not_a3)
    print("分离度:", sep)
    print("verdict:", out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
