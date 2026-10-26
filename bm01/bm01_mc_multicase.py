# -*- coding: utf-8 -*-
"""
(E-MC) 多案例泛化：判定码翻转是否**系统复现**，还是只是 BM-01 的个案轶事？
================================================================================
审稿人最可能的一句话："That's a single synthetic case study."

本实验用同一 1D 求解器扫 L × ΔQ 矩阵，每格报三个取样档的判定码：
  取样档 = {闸前单点 gate / 闸前 100 m 平均 gate_avg / 全池平均 mean}
  判定码 = peak_mean_mm <= 150 mm ? A1 : A3

网格策略：**固定物理网格 dx = 25 m**（ncell = L/25），
使不同池长的空间分辨率一致、闸前 100 m 窗恒为 4 格 —— 避免"网格变了解释不清楚"。

输出：raw/bm01_mc_multicase.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

DX_M = 25.0
LENGTHS = [2000.0, 6000.0, 12000.0, 18000.0]
DQS = [10.0, 45.0]
BAND_MM = 150.0
DURATION_S = 43200.0   # 12 h（论文已确立：4h 会假收敛）


def base_p(L):
    return {"pool_len_m": L, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": DURATION_S}


def code(peak_mm):
    return "A1" if peak_mm <= BAND_MM else "A3"


ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def main():
    t0 = time.time()
    res = {}
    for L in LENGTHS:
        ncell = int(round(L / DX_M))
        key = "L%d" % int(L)
        res[key] = {"pool_len_m": L, "ncell": ncell, "dx_m": round(L / ncell, 3),
                    "trapezoidal_note": "矩形 m=0", "dq_cases": {}}
        for dq in DQS:
            arms = {}
            for name, kw in ARMS:
                t1 = time.time()
                r = V2.sve_1d_v2(base_p(L), ncell, dq, m_slope=0.0, **kw)
                arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                              "final_mean_mm": round(r["final_mean_mm"], 2),
                              "early_peak_mm": round(r["early_peak_mm"], 2),
                              "judgement": code(r["peak_mean_mm"]),
                              "mass_err": float("%.2e" % r["mass_balance_rel_err"]),
                              "wall_s": round(time.time() - t1, 1)}
                print("[MC] L=%6.0f ncell=%4d dq=%5.1f %-7s peak=%8.1f mm %s (%.0fs)"
                      % (L, ncell, dq, name, arms[name]["peak_mean_mm"],
                         arms[name]["judgement"], arms[name]["wall_s"]), flush=True)
            arms["flip_point_vs_mean"] = (arms["mean"]["judgement"] == "A1"
                                          and arms["point"]["judgement"] == "A3")
            arms["flip_avg100_vs_mean"] = (arms["mean"]["judgement"] == "A1"
                                           and arms["avg100"]["judgement"] == "A3")
            arms["ratio_point_over_mean"] = round(
                arms["point"]["peak_mean_mm"] / max(arms["mean"]["peak_mean_mm"], 1e-9), 2)
            res[key]["dq_cases"]["dq%.0f" % dq] = arms

    # ---- 汇总：翻转是否系统复现 ----
    flips_p, flips_a, cells = 0, 0, 0
    for L in LENGTHS:
        for dq in DQS:
            a = res["L%d" % int(L)]["dq_cases"]["dq%.0f" % dq]
            cells += 1
            flips_p += 1 if a["flip_point_vs_mean"] else 0
            flips_a += 1 if a["flip_avg100_vs_mean"] else 0
    out = {
        "purpose": "E-MC 多案例泛化：判定码翻转是否系统复现",
        "method": "同一 1D 求解器；固定 dx=25 m；12 h；+ΔQ 阶跃；矩形断面",
        "lengths_m": LENGTHS, "dqs": DQS, "dx_m": DX_M, "duration_h": DURATION_S / 3600.0,
        "band_mm": BAND_MM, "results": res,
        "summary": {
            "cells": cells,
            "point_A3_cells": sum(1 for L in LENGTHS for dq in DQS
                                  if res["L%d" % int(L)]["dq_cases"]["dq%.0f" % dq]["point"]["judgement"] == "A3"),
            "avg100_A3_cells": sum(1 for L in LENGTHS for dq in DQS
                                   if res["L%d" % int(L)]["dq_cases"]["dq%.0f" % dq]["avg100"]["judgement"] == "A3"),
            "mean_A1_cells": sum(1 for L in LENGTHS for dq in DQS
                                 if res["L%d" % int(L)]["dq_cases"]["dq%.0f" % dq]["mean"]["judgement"] == "A1"),
            "flip_point_vs_mean_cells": flips_p,
            "flip_avg100_vs_mean_cells": flips_a,
        },
        "wall_total_s": round(time.time() - t0, 1),
    }
    s = out["summary"]
    out["verdict"] = ("翻转系统性复现：%d/%d 格 闸前单点 A3 而全池平均 A1；%d/%d 格 闸前100m平均 A3 而全池平均 A1"
                      % (s["flip_point_vs_mean_cells"], cells, s["flip_avg100_vs_mean_cells"], cells))

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_mc_multicase.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))
    print("summary:", json.dumps(s, ensure_ascii=False))
    print("verdict:", out["verdict"])


if __name__ == "__main__":
    main()
