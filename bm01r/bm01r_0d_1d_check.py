# -*- coding: utf-8 -*-
"""
(W2) BM-01R 的 0D/1D 对照：排除"0D 与 1D 的差异来自水面线形状"
================================================================================
三模型第三轮指出：E-PRED-R 用"0D 判 A1、1D 单点判 A3"来论证 BM-01R 是取样翻转，
但**0D 与 1D 全池平均的等价性（F1）此前只在矩形/梯形小案例上验证过**；
BM-01R 是 25.6 km 真实梯形渠段，"0D 漏掉水面线形状"这一刀必须先堵上。

做法：在 BM-01R 的**真实几何与真实参数**上，跑三个反馈通道：
  point   = 闸前单点（B0–B4 当年用的就是它）
  avg100  = 闸前 100 m 平均
  mean    = 全池空间平均
并与已有的 0D 集总结果（peak 45.31 mm → A1）逐项对照。

参数（严格取自 bm01r_lit.py，未改动）：b=20.75 m, m=2.25, h_target=4.50 m, L=25567 m,
q0=135 m³/s, ΔQ=+45, n=0.0150, kg=270 (=135/0.5), kp=1.60, ki=1.60e-3,
rate_limit=0.02/60 per s, ctrl_period=300 s, 24 h, ncell=256（与 B 组一致）。

输出：raw/bm01r_0d_1d_check.json
"""
import json, os, sys, math, time
sys.path.insert(0, r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收")
import bm01_sve1d_v2 as SVE2

HERE = os.path.dirname(os.path.abspath(__file__))
B_M, M_SLOPE, H_T, L_M = 20.75, 2.25, 4.50, 25567.0
Q0, DQ, NMAN, KG = 135.0, 45.0, 0.0150, 270.0
NCELL = 256
BAND_MM = 150.0
ZERO_D_PEAK_MM = 45.31        # raw/bm01r_literature.json → runs.A1.peak_mm
ZERO_D_CODE = "A1"

ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def P():
    return {"pool_width_m": B_M, "pool_len_m": L_M, "h_target": H_T, "q0": Q0, "kg": KG,
            "duration_s": 86400.0, "t_step_s": 900.0, "settle_window_s": 1800.0,
            "band_normal_m": 0.15, "a0": 0.50, "kp": 1.60, "ki": 1.60e-3,
            "ctrl_period_s": 300.0, "deadband_m": 0.010, "meas_resolution_m": 0.001,
            "rate_limit_per_s": 0.02 / 60.0, "band_alarm_m": 0.30, "settle_tol_m": 0.020}


def main():
    t0 = time.time()
    sec = SVE2.Section(B_M, M_SLOPE)
    s0 = SVE2.section_normal_slope(sec, H_T, Q0, NMAN)

    # C3 复核（用 BM-01R 真实 kg=270，而不是 BM-01 的 360）
    A_surf = (B_M + 2.0 * M_SLOPE * H_T) * L_M
    r_rate = 0.02 / 60.0
    dq_star = math.sqrt(KG * r_rate * A_surf * 0.15)
    T_band = A_surf * 0.15 / DQ
    T_slew = DQ / (KG * r_rate)

    out = {"purpose": "W2：BM-01R 真实梯形几何上的 0D/1D 对照，排除水面线形状差异",
           "config": {"b_m": B_M, "m_slope": M_SLOPE, "h_target_m": H_T, "L_m": L_M,
                      "q0": Q0, "dq": DQ, "n": NMAN, "kg": KG,
                      "rate_limit_per_s": r_rate, "ncell": NCELL, "dx_m": L_M / NCELL,
                      "duration_h": 24.0, "S0": s0, "band_mm": BAND_MM},
           "zero_d_reference": {"peak_mm": ZERO_D_PEAK_MM, "code": ZERO_D_CODE,
                                "source": "raw/bm01r_literature.json runs.A1"},
           "C3_with_real_kg": {"A_surf_m2": round(A_surf, 0), "T_band_s": round(T_band, 1),
                               "T_slew_s": round(T_slew, 1), "dq_star_m3s": round(dq_star, 1),
                               "dq_over_dqstar": round(DQ / dq_star, 4),
                               "predict_capacity_fail": bool(T_band < T_slew)},
           "arms": {}}

    print("BM-01R 真实几何 0D/1D 对照（ncell=%d, dx=%.1f m, 24 h, +%.0f m³/s）"
          % (NCELL, L_M / NCELL, DQ), flush=True)
    print("C3 复核（真实 kg=%.0f）：ΔQ*=%.1f m³/s，ΔQ=45 → ΔQ/ΔQ*=%.3f → %s"
          % (KG, dq_star, DQ / dq_star,
             "预测容量失败" if T_band < T_slew else "预测非容量失败"), flush=True)

    for name, kw in ARMS:
        r = SVE2.sve_1d_v2(P(), NCELL, DQ, m_slope=M_SLOPE, n=NMAN, s0=s0, **kw)
        pk = r["peak_mean_mm"]
        out["arms"][name] = {"peak_mean_mm": round(pk, 2),
                            "final_mean_mm": round(r["final_mean_mm"], 2),
                            "peak_gate_mm": round(r["peak_gate_mm"], 2),
                            "early_peak_mm": round(r["early_peak_mm"], 2),
                            "judgement": "A1" if pk <= BAND_MM else "A3",
                            "mass_err": float("%.2e" % r["mass_balance_rel_err"]),
                            "wall_s": round(time.time() - t0, 1)}
        print("[0D/1D] %-7s peak_mean=%8.2f mm  peak_gate=%8.2f mm  %s  (%.0fs)"
              % (name, pk, r["peak_gate_mm"], out["arms"][name]["judgement"],
                 out["arms"][name]["wall_s"]), flush=True)

    mean_pk = out["arms"]["mean"]["peak_mean_mm"]
    diff = mean_pk - ZERO_D_PEAK_MM
    rel = diff / ZERO_D_PEAK_MM * 100.0
    out["zero_d_vs_1d_mean"] = {"zero_d_mm": ZERO_D_PEAK_MM, "one_d_mean_mm": mean_pk,
                                "abs_diff_mm": round(diff, 2), "rel_diff_pct": round(rel, 2),
                                "codes_agree": out["arms"]["mean"]["judgement"] == ZERO_D_CODE}
    out["wall_total_s"] = round(time.time() - t0, 1)

    cond = abs(rel) <= 10.0 and out["arms"]["mean"]["judgement"] == ZERO_D_CODE
    if cond and out["arms"]["point"]["judgement"] == "A3":
        out["verdict"] = (
            "0D ≡ 1D 全池平均在 BM-01R 的**真实梯形几何**上成立：0D 集总 %s（%.2f mm）"
            "vs 1D 全池平均 %s（%.2f mm），相差 %.2f%%（%.2f mm）→ "
            "'0D 漏掉水面线形状'这一刀被排除。"
            "同时 1D 闸前单点判 %s（%.2f mm）而全池平均判 %s → "
            "**BM-01R 确为取样翻转**，且是 F2 受控反例在真实渠段上的第二次独立命中。"
            "C3（用真实 kg=%.0f）给出 ΔQ/ΔQ* = %.3f → 非容量失败，与之一致。"
            % (ZERO_D_CODE, ZERO_D_PEAK_MM, out["arms"]["mean"]["judgement"], mean_pk,
               abs(rel), abs(diff), out["arms"]["point"]["judgement"],
               out["arms"]["point"]["peak_mean_mm"], out["arms"]["mean"]["judgement"],
               KG, DQ / dq_star))
    else:
        out["verdict"] = ("0D 与 1D 全池平均在 BM-01R 上未达预期一致（相差 %.2f%%），"
                          "须逐项排查（水面线形状 / 时窗 / 网格），不得急于下结论。" % rel)

    os.makedirs(os.path.join(HERE, "raw"), exist_ok=True)
    jp = os.path.join(HERE, "raw", "bm01r_0d_1d_check.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
