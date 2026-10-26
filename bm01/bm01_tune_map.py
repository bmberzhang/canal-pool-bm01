# -*- coding: utf-8 -*-
"""
(E-TUNE) 增益平面扫描：闭合"翻转只是整定未覆盖工作域"这一刀
================================================================================
审稿人/仲裁方提的最硬质疑（R1）：
  "同一控制器对某些池稳定、对另一些池振荡，是**整定覆盖域**问题，
   不自动等于取样伪影 —— 也许只要重新整定，翻转就消失了。"

可证伪的实验设计：
  对两组取样（闸前 100 m 平均 gate_avg ／ 全池平均 mean）**在同一张增益平面 (kp, ki) 上**分别扫描，
  画出各自的"判定码地图"。判据：
   - 若存在**一片可观区域**使两组都判 A1 → R1 成立，翻转确可由重新整定消除；
   - 若 mean 的 A1 区域很大、而 gate_avg 的 A1 区域**为空或几乎为空**（即使把 kp 压到 0.125×、
     ki 压到 0.125× 都过不去）→ 翻转**不是整定不足**造成的，R1 被证伪。

设置：BM-01，240 格（dx=25 m），12 h，矩形，+45 m³/s。
输出：raw/bm01_tune_map.json（含两张 6×5 判定码地图 + 各自 A1 格数与峰值最小值）
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

KPS = [0.20, 0.40, 0.80, 1.60, 3.20, 6.40]
KIS = [0.0002, 0.0004, 0.0008, 0.0016, 0.0032]
NCELL = 240
DQ = 45.0
BAND_MM = 150.0


def P_of(kp, ki):
    return {"pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": kp, "ki": ki, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


ARMS = (("gate_avg", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def main():
    t0 = time.time()
    maps = {}
    for name, kw in ARMS:
        grid = {}
        best = 1e18
        for kp in KPS:
            row = {}
            for ki in KIS:
                t1 = time.time()
                r = V2.sve_1d_v2(P_of(kp, ki), NCELL, DQ, m_slope=0.0, **kw)
                pk = r["peak_mean_mm"]
                jd = "A1" if pk <= BAND_MM else "A3"
                best = min(best, pk)
                row["ki%.4f" % ki] = {"peak_mean_mm": round(pk, 2), "judgement": jd,
                                      "final_mean_mm": round(r["final_mean_mm"], 2)}
                print("[TUNE] %-8s kp=%.2f ki=%.4f -> %8.1f mm %s (%.0fs)"
                      % (name, kp, ki, pk, jd, time.time() - t1), flush=True)
            grid["kp%.2f" % kp] = row
        n_a1 = sum(1 for kp in KPS for ki in KIS
                   if grid["kp%.2f" % kp]["ki%.4f" % ki]["judgement"] == "A1")
        maps[name] = {"grid": grid, "a1_cells": n_a1, "total_cells": len(KPS) * len(KIS),
                      "min_peak_all_tunings_mm": round(best, 2),
                      "min_peak_location": min(
                          ((kp, ki) for kp in KPS for ki in KIS),
                          key=lambda t: grid["kp%.2f" % t[0]]["ki%.4f" % t[1]]["peak_mean_mm"])}
        mn = maps[name]["min_peak_location"]
        print("\n>>> %s：A1 格 %d/%d；全平面最小峰值 %.1f mm（kp=%.2f, ki=%.4f）\n"
              % (name, n_a1, len(KPS) * len(KIS), best, mn[0], mn[1]), flush=True)

    ga, mn = maps["gate_avg"], maps["mean"]
    overlap = sum(1 for kp in KPS for ki in KIS
                  if ga["grid"]["kp%.2f" % kp]["ki%.4f" % ki]["judgement"] == "A1"
                  and mn["grid"]["kp%.2f" % kp]["ki%.4f" % ki]["judgement"] == "A1")
    out = {"purpose": "E-TUNE 增益平面扫描：检验'翻转是否只是整定未覆盖工作域'",
           "config": {"ncell": NCELL, "dx_m": 6000.0 / NCELL, "duration_h": 12.0,
                      "step_dq": DQ, "band_mm": BAND_MM, "kps": KPS, "kis": KIS},
           "maps": maps,
           "overlap_A1_cells_both_arms": overlap,
           "total_cells": len(KPS) * len(KIS),
           "wall_total_s": round(time.time() - t0, 1)}

    # 判据
    if overlap == 0:
        verdict = ("R1 被证伪：在 30 个增益组合里【不存在】任何一个组合使闸前局部档也判 A1；"
                   "全池平均档 A1 格 %d/30。即使把 kp 压到 0.125×、ki 压到 0.125×，"
                   "局部档最小峰值仍为 %.1f mm（条带 %.0f mm）→ 翻转不是整定不足造成的。"
                   % (mn["a1_cells"], ga["min_peak_all_tunings_mm"], BAND_MM))
    elif overlap <= 3:
        verdict = ("R1 大致维持但边界极窄：仅 %d/30 个增益组合使两组都判 A1 → "
                   "重新整定虽存在理论通道，但稳定域窄到工程上不具可操作性。" % overlap)
    else:
        verdict = ("R1 成立：%d/30 个增益组合使两组都判 A1 → 翻转可由重新整定消除，"
                   "正文不得主张'取样伪影'。" % overlap)
    out["verdict"] = verdict

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_tune_map.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))
    print("重叠 A1 格 =", overlap, "/", len(KPS) * len(KIS))
    print("verdict:", verdict)


if __name__ == "__main__":
    main()
