# -*- coding: utf-8 -*-
"""
(E-TUNE-D) 谷底加密扫描：堵住 R1 的最后 4.4%
================================================================================
上一轮 E-TUNE 的结论是"局部档全平面最小峰值 156.6 mm > 条带 150 mm"，
但**只差 4.4%**，且网格是 2 倍步长。审稿人下一刀必然是：
  "你把 kp 从 0.20 到 6.40 只取了 6 个点，谷底很可能在 0.8~1.6 之间被你跨过去了；
   加密之后局部档是不是就掉进 150 mm 条带、翻转也就消失了？"

本实验就是这个疑问的**直接、可证伪回答**：
  - 在 E-TUNE 已定位的谷底邻域（kp 0.10~3.20，1.2~1.4 倍步长；ki 0.0001~0.0008）做细扫描；
  - 组数 12 kp × 6 ki = 72 组 × 2 档 = 144 次仿真；
  - 判据同样是"峰值全池平均 ≤ 150 mm 记 A1"。
输出判定：
  - 若 gate_avg 在**所有 72 组**仍 >150 mm，且谷底最小值与 E-TUNE 一致（≈156 mm）
    → 谷底是**平台**不是尖点，R1 彻底死亡，正文可写"任意整定下局部档都判 A3"；
  - 若出现 <150 mm 的格 → 必须诚实记录并下调结论（翻转存在整定通道）。

设置：BM-01，240 格（dx=25 m），12 h，矩形，+45 m³/s（与 E-TUNE 完全同构，可比）。
输出：raw/bm01_tune_valley.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

KPS = [0.10, 0.30, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 2.00, 2.40, 2.80, 3.20]
KIS = [0.0001, 0.0002, 0.0003, 0.0004, 0.0006, 0.0008]
NCELL = 240
DQ = 45.0
BAND_MM = 150.0

ARMS = (("gate_avg", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def P_of(kp, ki):
    return {"pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": kp, "ki": ki, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def main():
    t0 = time.time()
    maps = {}
    for name, kw in ARMS:
        grid, best, best_loc = {}, 1e18, None
        n_a1 = 0
        for kp in KPS:
            row = {}
            for ki in KIS:
                t1 = time.time()
                r = V2.sve_1d_v2(P_of(kp, ki), NCELL, DQ, m_slope=0.0, **kw)
                pk = r["peak_mean_mm"]
                jd = "A1" if pk <= BAND_MM else "A3"
                if jd == "A1":
                    n_a1 += 1
                if pk < best:
                    best, best_loc = pk, (kp, ki)
                row["ki%.5f" % ki] = {"peak_mean_mm": round(pk, 2), "judgement": jd,
                                      "final_mean_mm": round(r["final_mean_mm"], 2)}
                print("[VALLEY] %-8s kp=%.2f ki=%.5f -> %8.1f mm %s (%.0fs)"
                      % (name, kp, ki, pk, jd, time.time() - t1), flush=True)
            grid["kp%.2f" % kp] = row
        maps[name] = {"grid": grid, "a1_cells": n_a1, "total_cells": len(KPS) * len(KIS),
                      "min_peak_all_tunings_mm": round(best, 2),
                      "min_peak_location": [best_loc[0], best_loc[1]]}
        print("\n>>> %s：A1 格 %d/%d；全平面最小峰值 %.2f mm（kp=%.2f, ki=%.5f）\n"
              % (name, n_a1, len(KPS) * len(KIS), best, best_loc[0], best_loc[1]), flush=True)

    ga, mn = maps["gate_avg"], maps["mean"]
    overlap = sum(1 for kp in KPS for ki in KIS
                  if ga["grid"]["kp%.2f" % kp]["ki%.5f" % ki]["judgement"] == "A1"
                  and mn["grid"]["kp%.2f" % kp]["ki%.5f" % ki]["judgement"] == "A1")
    # 谷底平台性：统计 gate_avg 在最小峰值 1.10 倍以内的格数，看是尖点还是平台
    floor = ga["min_peak_all_tunings_mm"]
    plateau = sum(1 for kp in KPS for ki in KIS
                  if ga["grid"]["kp%.2f" % kp]["ki%.5f" % ki]["peak_mean_mm"] <= floor * 1.10)
    margin = (floor - BAND_MM) / BAND_MM * 100.0

    out = {"purpose": "E-TUNE-D 谷底加密扫描：检验 R1 是否只是'增益网格太粗'造成的假象",
           "config": {"ncell": NCELL, "dx_m": 6000.0 / NCELL, "duration_h": 12.0,
                      "step_dq": DQ, "band_mm": BAND_MM, "kps": KPS, "kis": KIS},
           "maps": maps,
           "overlap_A1_cells_both_arms": overlap,
           "gate_avg_floor_mm": floor,
           "gate_avg_floor_margin_pct": round(margin, 2),
           "gate_avg_plateau_cells_within_10pct_floor": plateau,
           "total_cells": len(KPS) * len(KIS),
           "wall_total_s": round(time.time() - t0, 1)}

    if ga["a1_cells"] == 0:
        verdict = ("R1 彻底证伪（加密后依然）：局部档在 %d 组增益里 A1 格 = 0，"
                   "全平面最小峰值 %.1f mm，仍比 %.0f mm 条带高 %.1f%%；"
                   "最小峰值邻域 %.1f 倍以内有 %d/%d 格 → 谷底是**平台**而非尖点，"
                   "'加密网格就能掉进条带'的通道不存在。"
                   % (len(KPS) * len(KIS), floor, BAND_MM, margin, 1.10, plateau,
                      len(KPS) * len(KIS)))
    else:
        verdict = ("警告：加密后局部档出现 %d/%d 个 A1 格（最小 %.1f mm）→ "
                   "翻转确实存在整定通道，正文必须下调为'默认整定下翻转、但可通过重整定部分消除'。"
                   % (ga["a1_cells"], len(KPS) * len(KIS), floor))
    out["verdict"] = verdict

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_tune_valley.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))
    print("gate_avg A1=%d/%d  最小 %.2f mm（余量 %.2f%%）  重叠 A1=%d"
          % (ga["a1_cells"], len(KPS) * len(KIS), floor, margin, overlap))
    print("verdict:", verdict)


if __name__ == "__main__":
    main()
