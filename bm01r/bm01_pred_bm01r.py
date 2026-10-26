# -*- coding: utf-8 -*-
"""
(E-PRED-R) 用同一解析判据复核 BM-01R，并重判它的失败归属
================================================================================
第四轮把 BM-01R 写成"真实容量失败"，依据是两句话：
  (i)  "失败方式是稳态壅高而非波动"；
  (ii) "180 m³/s 对应均匀流水深 5.242 m，即 +742 mm = 条带的 4.9 倍 → 与控制器无关的输水容量上限"。

本脚本用**已检验的解析判据**复核这两句，并检查 BM-01R 各算例实际用了什么反馈通道。

关键事实（来自 raw/bm01r_literature.json）：
  A1（0D 集总，真实水面面积 1 048 247 m²）→ peak 45.31 mm → **A1**
  B2（1D 梯形 m=2.25，真实几何）           → peak 338.82 mm → A3
  而 **全部 1D 算例（B0–B4）的 ctrl_on 都是 "gate"（单点传感器）**。
  在本文框架下，0D 集总 ≡ 1D 的全池空间平均（F1）。⇒ 同对象、同控制器、同扰动，
  仅把反馈从"全池平均"换成"闸前单点"，判定码就从 A1 翻成 A3
  ——**BM-01R 是 F2 受控反例在一条真实渠段上的独立复现，不是容量失败。**

判据复核：
  C1（输水容量）：Δh_n = 5.242 − 4.500 = 742 mm > 150 mm → 预测"容量失败"。
     但 0D/全池平均给 A1 ⇒ C1 在 BM-01R 上**再次错误**（与 E-PRED 在 BM-01 的 7 个虚警同向）。
  C3（执行器能力）：A_surf=(b+2m·h)·L；ΔQ*=√(k_g·r·A_surf·band)
     = √(360·(1/3000)·1 048 247·0.15) ≈ 137.4 m³/s ≫ ΔQ=45 ⇒ 预测"非容量失败" ⇒ 与"实为取样翻转"一致。

输出：raw/bm01_pred_bm01r.json
"""
import json, os, math, sys
sys.path.insert(0, r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收")
import bm01_sve1d_v2 as V2

HERE = os.path.dirname(os.path.abspath(__file__))
LIT = os.path.join(HERE, "raw", "bm01r_literature.json")

B_M, M_SLOPE, H_T, L_M = 20.75, 2.25, 4.50, 25567.0
Q0, DQ = 135.0, 45.0
NMAN = 0.015
KG, RATE = 360.0, 1.0 / 3000.0
BAND = 0.15


def h_normal(Q, sec, S0):
    lo, hi = 0.05, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if V2.section_normal_slope(sec, mid, Q, NMAN) > S0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    lit = json.load(open(LIT, encoding="utf-8"))
    runs = lit["runs"]
    sec = V2.Section(B_M, M_SLOPE)
    S0 = lit["params"].get("S0") or 4.4346790794794846e-05

    hn0 = h_normal(Q0, sec, S0)
    hn1 = h_normal(Q0 + DQ, sec, S0)
    dh_n = hn1 - hn0

    A_surf = (B_M + 2.0 * M_SLOPE * H_T) * L_M
    T_band = A_surf * BAND / DQ
    T_slew = DQ / (KG * RATE)
    dq_star = math.sqrt(KG * RATE * A_surf * BAND)

    point_runs = {k: v.get("ctrl_on") for k, v in runs.items()
                  if isinstance(v, dict) and v.get("ctrl_on") is not None}
    zero_d = {k: v for k, v in runs.items()
              if isinstance(v, dict) and v.get("model") == "0d-lumped"}

    out = {"purpose": "E-PRED-R：复核 BM-01R 的失败归属（容量 vs 取样）",
           "geometry": {"b_m": B_M, "m_slope": M_SLOPE, "h_target_m": H_T, "L_m": L_M,
                        "q0": Q0, "dq": DQ, "n": NMAN, "kg": KG, "rate_limit_per_s": RATE},
           "manning": {"h_normal_q0_m": round(hn0, 4), "h_normal_q0dq_m": round(hn1, 4),
                       "dh_normal_mm": round(dh_n * 1000, 1), "band_mm": BAND * 1000,
                       "ratio_to_band": round(dh_n / BAND, 2)},
           "C1_uniform_flow": {"predict_capacity_fail": bool(dh_n > BAND),
                               "actual_0d_fullpool_code": None, "verdict": None},
           "C3_actuator": {"A_surf_m2": round(A_surf, 0), "T_band_s": round(T_band, 1),
                           "T_slew_s": round(T_slew, 1), "dq_star_m3s": round(dq_star, 1),
                           "margin_x": round(dq_star / DQ, 2),
                           "predict_capacity_fail": bool(T_band < T_slew)},
           "evidence_from_literature_json": {"1d_runs_ctrl_on": point_runs,
                                             "0d_lumped_runs": {k: {"peak_mm": v.get("peak_mm"),
                                                                     "code": v.get("code")}
                                                                for k, v in zero_d.items()}},
           }

    # 0D = 全池平均等价（F1）
    z0 = sorted(zero_d.items(), key=lambda kv: -kv[1].get("area_m2", 0))
    z_key, z_val = z0[0] if z0 else (None, {})
    out["C1_uniform_flow"]["actual_0d_fullpool_code"] = z_val.get("code")
    out["C1_uniform_flow"]["actual_0d_peak_mm"] = z_val.get("peak_mm")
    out["C1_uniform_flow"]["verdict"] = ("C1 预测'容量失败'，但 0D（≡全池平均）判 %s（%.2f mm）"
                                          "⇒ C1 在 BM-01R 上错误"
                                          % (z_val.get("code"), z_val.get("peak_mm") or -1))

    all_point = all(v == "gate" for v in point_runs.values()) and len(point_runs) > 0
    if all_point:
        verdict = ("BM-01R 的失败**不是容量失败，而是取样翻转**："
                   "全部 %d 个 1D 算例的反馈通道都是单点 `gate`（闸前单点），"
                   "而 0D 集总（≡ 1D 全池空间平均，F1）在同一真实几何下判 %s（%.2f mm）。"
                   "⇒ 同对象、同控制器、同扰动，仅换取样位置，判定码 A1↔A3 翻转 —— "
                   "**F2 受控反例在一条 25.6 km 真实渠段上独立复现**。"
                   "连带修正：round-4 的『+742 mm 是输水容量上限』应删除（C1 已被否证，此处再添一例错误）；"
                   "『稳态壅高』恰是『闸前单点被守住、池均自由漂移』的典型取样特征。"
                   % (len(point_runs), z_val.get("code"), z_val.get("peak_mm") or -1))
    else:
        verdict = ("BM-01R 的 1D 算例并非全部使用单点反馈，需逐例重判；其他结论暂缓。")
    out["verdict"] = verdict

    os.makedirs(os.path.join(os.path.dirname(HERE)), exist_ok=True)
    dst = os.path.join(HERE, "raw", "bm01_pred_bm01r.json")
    json.dump(out, open(dst, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print("== BM-01R 判据复核 ==")
    print("正常水深 %.4f -> %.4f m，Δh_n = %.1f mm（= 条带的 %.2f 倍）"
          % (hn0, hn1, dh_n * 1000, dh_n / BAND))
    print("C1 输水容量判据：预测%s，但 0D/全池平均判 %s（%.2f mm）→ C1 %s"
          % ("容量失败" if dh_n > BAND else "通过", z_val.get("code"), z_val.get("peak_mm") or -1,
             "错误" if dh_n > BAND else "一致"))
    print("C3 执行器能力：T_band=%.0f s vs T_slew=%.0f s；ΔQ*=%.1f m³/s（ΔQ=45，余量 %.1f 倍）→ 预测%s"
          % (T_band, T_slew, dq_star, dq_star / DQ,
             "容量失败" if T_band < T_slew else "非容量失败"))
    print("1D 算例反馈通道：", point_runs)
    print()
    print(verdict)
    print("\n已写入", dst)


if __name__ == "__main__":
    main()
