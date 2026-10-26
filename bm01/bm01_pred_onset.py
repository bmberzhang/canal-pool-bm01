# -*- coding: utf-8 -*-
"""
(W0-b) 因一的经验起爆点：在多个池长上找"三档全 A3"的真实边界
================================================================================
W0 首轮扫描发现：L=2000 的真实起爆点在 ΔQ ∈ (14, 16] m³/s，
而解析式 C3 预测 ΔQ* = 42.4 m³/s —— **高估 2.8 倍**。
⇒ C3 不能作为判定阈值（三模型已预警"单点正类不构成判据"，此处被证实），
  而且它连"必要条件"都不是（起爆格不满足 C3 不等式）。

本脚本取第二、第三个池长的起爆点，检验有无经验规律：
  - 若 ΔQ_onset ∝ √(A_surf)，则与 C3 同标度、只是系数不同（可改写成带系数的诊断式）；
  - 若标度不符，则认定"因一目前只能操作性识别（三档全失败），无可靠解析判据"，如实写进局限。

设置：固定 dx = 25 m，12 h，矩形，三档取样（point / avg100 / mean）。
输出：raw/bm01_pred_onset.json
"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

DX_M = 25.0
BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
CASES = ((2000.0, [10, 12, 13, 14, 15, 16]),          # 收窄复跑，确认首轮起爆点
         (6000.0, [14, 16, 20, 24, 28, 32, 36, 40, 45, 50]),
         (12000.0, [24, 32, 40, 45, 50, 60, 70]))
ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def base_p(L):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def main():
    t0 = time.time()
    out = {"purpose": "W0-b 多池长的因一起爆点：检验 C3 的标度是否正确",
           "config": {"dx_m": DX_M, "duration_h": 12.0, "band_mm": BAND_MM,
                      "kg": KG, "rate_limit_per_s": RATE},
           "by_length": {}, "C3_reference": {}, "wall_total_s": None}

    for L, dqs in CASES:
        ncell = int(round(L / DX_M))
        A_surf = B0 * L
        dq_star = math.sqrt(KG * RATE * A_surf * 0.15)
        rows = []
        for dq in dqs:
            arms = {}
            for name, kw in ARMS:
                r = V2.sve_1d_v2(base_p(L), ncell, float(dq), m_slope=0.0, **kw)
                arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                              "final_mean_mm": round(r["final_mean_mm"], 2),
                              "early_peak_mm": round(r["early_peak_mm"], 2)}
            n_a3 = sum(1 for a in ARMS if arms[a[0]]["peak_mean_mm"] > BAND_MM)
            rows.append({"dq": dq, "peaks_mm": {k: v["peak_mean_mm"] for k, v in arms.items()},
                         "mean_arm_peak_mm": arms["mean"]["peak_mean_mm"],
                         "mean_arm_final_mm": arms["mean"]["final_mean_mm"],
                         "n_A3_arms": n_a3, "ALL_A3": n_a3 == 3})
            print("[ONSET] L=%6.0f ΔQ=%4.1f | point=%8.1f avg100=%8.1f mean=%8.1f | A3=%d %s"
                  % (L, dq, arms["point"]["peak_mean_mm"], arms["avg100"]["peak_mean_mm"],
                     arms["mean"]["peak_mean_mm"], n_a3,
                     "ALL_A3" if n_a3 == 3 else "     "), flush=True)
        a3 = [r["dq"] for r in rows if r["ALL_A3"]]
        ok = [r["dq"] for r in rows if not r["ALL_A3"]]
        onset_bound = [max(ok), min(a3)] if (ok and a3) else None
        out["by_length"]["L%d" % int(L)] = {
            "pool_len_m": L, "ncell": ncell, "A_surf_m2": A_surf,
            "C3_dq_star_m3s": round(dq_star, 2),
            "onset_bracket_m3s": onset_bound,
            "onset_mid_m3s": round(sum(onset_bound) / 2.0, 2) if onset_bound else None,
            "ratio_C3_over_actual": round(dq_star / (sum(onset_bound) / 2.0), 2) if onset_bound else None,
            "rows": rows}
        print("\n>>> L=%.0f：全 A3 起爆区间 = %s；C3 预测 ΔQ*=%.1f → 比值 %.2f\n"
              % (L, onset_bound, dq_star,
                 (dq_star / (sum(onset_bound) / 2.0)) if onset_bound else float('nan')), flush=True)
        if not (ok and a3):
            print("    （该池长在扫描区间内未完整跨越起爆点，需扩区间）\n", flush=True)

    # 标度检验：ΔQ_onset vs √A_surf
    pts = [(v["A_surf_m2"], v["onset_mid_m3s"]) for v in out["by_length"].values()
           if v["onset_mid_m3s"]]
    if len(pts) >= 2:
        pts.sort()
        ratios = [(math.sqrt(a), q) for a, q in pts]
        k = sum(q / s for s, q in ratios) / len(ratios)
        resid = [abs(q / s - k) / k for s, q in ratios]
        out["scaling_check"] = {
            "model": "ΔQ_onset = k·√(A_surf·band·k_g·r)  (即 C3 形式)",
            "fitted_k": round(k, 4), "rel_residuals": [round(x, 4) for x in resid],
            "sqrt_A_points": [(round(s, 1), q) for s, q in ratios],
            "consistent": bool(max(resid) < 0.15)}
    out["wall_total_s"] = round(time.time() - t0, 1)

    onset_txt = "; ".join("L=%.0f: ΔQ∈%s" % (v["pool_len_m"], v["onset_bracket_m3s"])
                          for v in out["by_length"].values() if v["onset_bracket_m3s"])
    if len(pts) >= 2 and out["scaling_check"]["consistent"]:
        out["verdict"] = ("因一起爆点存在且可定位（%s）。但与 √A_surf 的标度检验**通过**"
                          "（拟合系数 k=%.3f，最大相对残差 %.1f%%）→ C3 的**形式**正确、"
                          "**系数**错（原式隐含 k=1，实测 k≈%.3f）→ 应改写为带系数的诊断式，"
                          "并注明系数由本案例族标定。" % (onset_txt, out["scaling_check"]["fitted_k"],
                                                         max(resid) * 100, k))
    else:
        out["verdict"] = ("因一起爆点存在且可定位（%s），但 C3 的闭式（隐含系数 k=1）"
                          "预测值高估数倍，且 √A_surf 标度检验%s → "
                          "**不得把 C3 写成判定阈值**；因一目前只能操作性识别"
                          "（三档全失败），如实写入局限。" % (onset_txt,
                          "未通过" if len(pts) >= 2 else "样本不足无法检验"))

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_pred_onset.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
