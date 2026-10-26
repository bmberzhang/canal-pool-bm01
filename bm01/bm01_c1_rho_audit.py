# -*- coding: utf-8 -*-
"""
(C1) rho 单参数标度律的「执行器速率稳健性审计」
================================================================================
背景：三模型第四轮把此列为**最高优先动作**，理由——
  台账 §1.7 证据：极限环周期 = 池体共振周期的 3.7–6.1 倍，**由执行器时间尺度设定**
  ⇒ 系统至少有两个独立时间尺度（控制采样 T_ctrl、执行器速率）。
  rho = T_ctrl·c/L 只含"控制采样/波行时间"，若它把系统压成单参数，
  则"L=2000 与 L=6000 曲线在 rho 上重合"可能只是巧合。
  **审稿人可一句"你的 rho 里藏着第三个时间尺度"打穿。**

【预注册判据】（先于计算写死，写在脚本开头）
  对每个**固定 rho 档**，在 4 个执行器速率倍数 {0.5, 1, 2, 4} 上重算起爆与跃升倍数 R：
    - 若 4 个速率下该档的"形态标签"（边际型 R<3 / 失稳型 R≥3）**全部相同**，且起爆点存在
      ⇒ 该档 **稳健**；
    - 若任一档在速率变化下**标签翻转**，或起爆点在某速率下消失/退化为第三区制
      ⇒ **审计失败**（rho 单参数标度作废，只能报告"两池长经验重合"）。
  审计通过条件：**全部档位稳健**。

设置：与 W4''-fix 同构（240/80 格、12 h、kp=1.6 ki=0.0016、矩形、全池平均反馈）。
  默认速率 = 1/3000 m/s（既有 RATE）；倍数 0.5/1/2/4 即 1/6000 … 4/3000。
档位：L=2000（rho 0.44 / 0.88 / 1.32）+ L=6000（rho 0.44 / 0.88）——覆盖两池长的关键区制。
输出：raw/bm01_c1_rho_audit.json
"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BAND_MM = 150.0
B0, HT, Q0, KG = 50.0, 7.0, 180.0, 360.0
RATE0 = 1.0 / 3000.0
C_EFF = math.sqrt(9.81 * HT) + Q0 / (B0 * HT)
DELTA = 2.0
SEARCH_LO, SEARCH_HI = 0.5, 120.0

LEVELS = [(2000.0, 100.0), (2000.0, 200.0), (2000.0, 300.0), (6000.0, 300.0), (6000.0, 600.0)]
RATE_MULT = [0.5, 1.0, 2.0, 4.0]


def base_p(L, period, rate):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": Q0, "h_target": HT, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": period, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": rate, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def mean_peak(L, ncell, period, rate, dq):
    return V2.sve_1d_v2(base_p(L, period, rate), ncell, dq, m_slope=0.0, ctrl_on="mean")["peak_mean_mm"]


def find_onset(L, ncell, period, rate):
    if mean_peak(L, ncell, period, rate, SEARCH_HI) <= BAND_MM:
        return None
    if mean_peak(L, ncell, period, rate, SEARCH_LO) > BAND_MM:
        return SEARCH_LO  # 第三区制：零扰动即失稳
    lo, hi = SEARCH_LO, SEARCH_HI
    while hi - lo > 1.0:
        mid = 0.5 * (lo + hi)
        if mean_peak(L, ncell, period, rate, mid) > BAND_MM:
            hi = mid
        else:
            lo = mid
    return hi


def main():
    t0 = time.time()
    rows = []
    for (L, period) in LEVELS:
        ncell = int(round(L / 25.0))
        rho = period / (L / C_EFF)
        for rm in RATE_MULT:
            rate = RATE0 * rm
            t1 = time.time()
            onset = find_onset(L, ncell, period, rate)
            base_mp = mean_peak(L, ncell, period, rate, SEARCH_LO)
            level_valid = base_mp <= BAND_MM
            if onset is None:
                R = None
                mp_on = mp_prev = None
                onset_txt = "无起爆(ΔQ≤120 内池均恒 <150)"
            elif abs(onset - SEARCH_LO) < 1e-9:
                R = None
                mp_on = mp_prev = None
                onset_txt = "第三区制(零扰动即失稳)"
            else:
                pv = max(SEARCH_LO, onset - DELTA)
                mp_on = mean_peak(L, ncell, period, rate, onset)
                mp_prev = mean_peak(L, ncell, period, rate, pv)
                R = mp_on / mp_prev if mp_prev > 0 else float('inf')
                onset_txt = "%.1f" % onset
            label = None if R is None else ("失稳型" if R >= 3.0 else "边际型")
            rows.append({"L_m": L, "ctrl_period_s": period, "rho": round(rho, 3),
                         "rate_mult": rm, "rate_limit_per_s": rate,
                         "onset_dq": None if onset is None else round(onset, 2),
                         "baseline_peak_mm": round(base_mp, 2), "level_valid": bool(level_valid),
                         "mean_peak_onset_mm": None if mp_on is None else round(mp_on, 2),
                         "mean_peak_prev_mm": None if mp_prev is None else round(mp_prev, 2),
                         "jump_ratio_R": None if R is None else round(R, 3),
                         "form_label": label, "onset_txt": onset_txt})
            print("[C1] L=%5.0f T=%5.0fs rho=%.2f rate×%.1f | 起爆 %s | 基线 %7.1f mm | R=%s | %s"
                  % (L, period, rho, rm, onset_txt, base_mp,
                     "None" if R is None else "%.2f" % R, label or ("无效档" if not level_valid else "?")),
                  flush=True)

    # 按 rho 分组，检查标签在速率维度是否恒定
    by_rho = {}
    for r in rows:
        by_rho.setdefault(r["rho"], []).append(r)
    audits = []
    for rho, rs in sorted(by_rho.items()):
        labels = set(x["form_label"] for x in rs)
        robust = (len(labels) == 1) and (None not in labels)
        audits.append({"rho": rho, "labels": sorted([str(x) for x in labels]),
                       "onset_dq_by_rate": {str(x["rate_mult"]): x["onset_dq"] for x in rs},
                       "R_by_rate": {str(x["rate_mult"]): x["jump_ratio_R"] for x in rs},
                       "robust": bool(robust)})
        print("[C1-AUDIT] rho=%.2f → 标签集合 %s | %s" % (rho, sorted([str(l) for l in labels]),
              "稳健" if robust else "**翻转/退化**"), flush=True)

    audit_pass = all(a["robust"] for a in audits)
    out = {"purpose": "C1：rho 单参数标度律的执行器速率稳健性审计",
           "preregistered_criterion": ("固定 rho 下，4 个执行器速率倍数 {0.5,1,2,4} 上，"
                                       "形态标签（边际型 R<3 / 失稳型 R≥3）须全部相同且起爆点存在；"
                                       "任一档翻转或退化 ⇒ 审计失败，rho 单参数标度作废。"),
           "settings": {"cells": "L/25 m", "duration_s": 43200, "kp": 1.6, "ki": 0.0016,
                        "feedback": "mean（全池平均）", "rate0_per_s": RATE0,
                        "rates_mult": RATE_MULT},
           "rows": rows, "by_rho_audit": audits,
           "audit_pass": bool(audit_pass),
           "verdict": ("**C1 审计" + ("通过**：rho 分区对执行器速率稳健（4 档速率下同 rho 的形态标签不变），"
                        "rho 作为单参数标度律的可信度得到支持。"
                        if audit_pass else
                        "**失败**：存在固定 rho 档在改变执行器速率时形态标签翻转或起爆点退化 ⇒ "
                        "rho=T_ctrl·c/L 不是充分单参数，第三时间尺度（执行器速率）未被包含，"
                        "只能报告『L=2000 与 L=6000 在给定执行器参数下的经验重合』。")),
           "wall_total_s": round(time.time() - t0, 1)}

    json.dump(out, open("raw/bm01_c1_rho_audit.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 raw/bm01_c1_rho_audit.json（%.0fs）" % out["wall_total_s"])


if __name__ == "__main__":
    main()
