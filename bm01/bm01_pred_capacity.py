# -*- coding: utf-8 -*-
"""
(E-PRED) 两因框架里"因一"到底是什么：纯解析判据的预测效度检验
================================================================================
第四轮把"因一"写成"输水容量（阶跃后所需均匀流水深是否越带）"。
但这**从未被检验过**，而且有一个明显矛盾：
  均匀流水深变化 Δh_n 只取决于 ΔQ，与池长 L 无关；
  而 E-MC 实测里 L=2000/ΔQ=45 全部 A3，L≥6000/ΔQ=45 却翻转或通过。
  ⇒ 如果 Δh_n 真是因一，就无法解释"同 ΔQ、不同 L"的差异。

本实验用**纯算术**（不跑求解器）比较两个候选解析判据，与 E-MC 实测逐格对照：

  C1（输水容量/均匀流可达性）：Δh_n(ΔQ) = h_n(Q0+ΔQ) − h_n(Q0) > 0.15 m ?
     h_n 由曼宁公式反解（实测断面）。预测"越带者必全 A3"。

  C3（执行器能力/蓄量可吸收性）：
     阶跃填满可用条带蓄量的时间  T_band = B·L·band / ΔQ
     执行器追上阶跃所需的最短时间 T_slew = ΔQ / (k_g · r)，  r = 速率限制（满行程 1/3000 每秒）
     若 T_band < T_slew ⇒ 执行器物理上追不上 ⇒ 容量型失败
     等价阈值：  ΔQ² > k_g · r · A_surf · band       （A_surf = B·L）
     **与控制器增益无关，与取样位置无关** —— 正是"因一"要求的性质。

判定"预测效度"：
  - 若 C1 预测错误被数据否证 → "因一 = 输水容量"的说法必须撤回；
  - 若 C3 8 格全中 → 因一应改述为"执行器能力（可控性）限制"，且给出可算的阈值。

输出：raw/bm01_pred_capacity.json
"""
import json, os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BAND = 0.15
B0 = 50.0
KG = 360.0
RATE_LIMIT = 1.0 / 3000.0     # 满行程 3000 s
NMAN = V2.MANNING_N


def h_normal(Q, sec, S0):
    """给定流量与纵坡，二分反解正常水深。
    注意：S( h ) = (nQ/(A·R^{2/3}))² 随 h **单调递减**（A 随 h 增大）。
    """
    lo, hi = 0.05, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if V2.section_normal_slope(sec, mid, Q, NMAN) > S0:
            lo = mid          # S 太大 ⇒ h 还不够
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    sec = V2.Section(B0, 0.0)
    S0 = V2.section_normal_slope(sec, 7.0, 180.0, NMAN)
    h0 = h_normal(180.0, sec, S0)
    print("正常水深 h_n(180) = %.6f m（应=7.000000，验证反解）" % h0, flush=True)

    emc = json.load(open(os.path.join("raw", "bm01_mc_multicase.json"), encoding="utf-8"))
    RES = emc["results"]

    out = {"purpose": "E-PRED：检验两因框架中『因一』的解析判据是否预测 E-MC 实测",
           "h_n_180_m": round(h0, 6), "S0": S0, "band_m": BAND,
           "kg": KG, "rate_limit_per_s": RATE_LIMIT,
           "cells": [], "criteria": {}}

    cells = []
    for lk in sorted(RES.keys(), key=lambda s: int(s[1:])):
        L = RES[lk]["pool_len_m"]
        ncell = RES[lk]["ncell"]
        for dk in sorted(RES[lk]["dq_cases"].keys(), key=lambda s: float(s.replace("dq", ""))):
            arms = RES[lk]["dq_cases"][dk]
            dq = float(dk.replace("dq", ""))
            def pk(a):
                return arms[a]["peak_mean_mm"] if a in arms else None
            # 三档峰值（E-MC 的键名 point/avg100/mean）
            p_point = pk("point"); p_avg = pk("avg100"); p_mean = pk("mean")
            n_a3 = sum(1 for x in (p_point, p_avg, p_mean) if x is not None and x > BAND * 1000)
            actual = "ALL_A3" if n_a3 == 3 else ("FLIP" if p_mean is not None and p_mean <= BAND * 1000
                                                 and n_a3 >= 1 else "ALL_A1")
            # C1
            dh_n = h_normal(180.0 + dq, sec, S0) - h0
            c1_fail = dh_n > BAND
            # C3
            A_surf = B0 * L
            T_band = A_surf * BAND / dq
            T_slew = dq / (KG * RATE_LIMIT)
            c3_fail = T_band < T_slew
            dq_star = math.sqrt(KG * RATE_LIMIT * A_surf * BAND)
            cells.append({"L_m": L, "dq": dq, "ncell": ncell,
                          "peaks_mm": {"point": p_point, "avg100": p_avg, "mean": p_mean},
                          "actual": actual, "n_A3_arms": n_a3,
                          "C1_dh_normal_mm": round(dh_n * 1000, 2), "C1_predict_fail": c1_fail,
                          "C3_T_band_s": round(T_band, 1), "C3_T_slew_s": round(T_slew, 1),
                          "C3_predict_fail": c3_fail, "C3_dq_star_m3s": round(dq_star, 2),
                          "C3_margin_pct": round((dq / dq_star - 1) * 100, 1)})
            print("[PRED] L=%6.0f dq=%4.0f | 实测 %-7s (A3档数=%d) | C1 Δh_n=%7.1f mm %s | "
                  "C3 T_band=%7.1fs T_slew=%6.1fs %s (ΔQ*=%.1f)"
                  % (L, dq, actual, n_a3, dh_n * 1000, "FAIL" if c1_fail else "ok",
                     T_band, T_slew, "FAIL" if c3_fail else "ok", dq_star), flush=True)

    n_all_a3 = sum(1 for c in cells if c["actual"] == "ALL_A3")
    def score(key):
        tp = sum(1 for c in cells if c[key] and c["actual"] == "ALL_A3")
        fp = sum(1 for c in cells if c[key] and c["actual"] != "ALL_A3")
        fn = sum(1 for c in cells if (not c[key]) and c["actual"] == "ALL_A3")
        return {"hit": tp, "false_alarm": fp, "miss": fn, "exact": (fp == 0 and fn == 0)}

    s1, s3 = score("C1_predict_fail"), score("C3_predict_fail")
    out["criteria"] = {"C1_uniform_flow": s1, "C3_actuator": s3}
    out["n_cells"] = len(cells)
    out["n_all_A3"] = n_all_a3

    if s3["exact"] and not s1["exact"]:
        verdict = ("E-PRED 结论：**因一不是『输水容量』，而是『执行器能力（可控性）限制』**。"
                   "C1（均匀流水深越带）把全部 4 个 ΔQ=45 格都判『容量失败』，"
                   "但实测只有 L=2000 那一格全 A3 → C1 被数据否证（虚警 %d）。"
                   "C3（T_band < T_slew，即 ΔQ² > k_g·r·A_surf·band）在 %d 格上**全中**"
                   "（命中 %d、虚警 %d、漏判 %d），且阈值 ΔQ* = √(k_g·r·A_surf·band) 是**闭式解析**、"
                   "与控制器增益和取样位置都无关 —— 正符合『因一』要求的性质。"
                   "⇒ 两因框架须改述：因一 = 执行器能力限制（池子相对扰动太小 / 执行器太慢），"
                   "因二 = 反馈取样尺度。"
                   % (s1["false_alarm"], len(cells), s3["hit"], s3["false_alarm"], s3["miss"]))
    elif s3["exact"] and s1["exact"]:
        verdict = "C1 与 C3 在本数据上均完全命中，需更大扫描范围区分二者。"
    else:
        verdict = ("两个解析判据都未能完全解释实测（C1 虚警 %d 漏判 %d；C3 虚警 %d 漏判 %d）"
                   "→ 因一暂无可靠解析判据，须如实声明为局限。"
                   % (s1["false_alarm"], s1["miss"], s3["false_alarm"], s3["miss"]))
    out["cells"] = cells
    out["verdict"] = verdict

    os.makedirs("raw", exist_ok=True)
    json.dump(out, open(os.path.join("raw", "bm01_pred_capacity.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n" + verdict)


if __name__ == "__main__":
    main()
