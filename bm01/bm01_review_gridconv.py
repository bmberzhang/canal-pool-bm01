# -*- coding: utf-8 -*-
"""回应审稿 3.2：valley-floor cell 的网格收敛 + 时间步收敛。

对 gate_avg 档谷底 (kp=3.20, ki=0.0003) 与 mean 档最佳 (kp=2.00, ki=0.0008)，
跑 ncell = 120/240/480/960（dx = 50/25/12.5/6.25 m）网格收敛，
以及 CFL 减半的时间步收敛，回答：
  "154.23 mm 只比 150 mm 高 2.8%，是否被网格/时间步误差翻转？"

输出：raw/bm01_review_gridconv.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

def P_of(kp, ki):
    return {"pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": kp, "ki": ki, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}

CELLS = [("gate_avg_floor", "gate_avg", 3.20, 0.0003),
         ("mean_best", "mean", 2.00, 0.0008)]

def run():
    out = {"purpose": "valley-floor grid & timestep convergence (review 3.2)",
           "cells": []}
    t0 = time.time()
    for name, arm, kp, ki in CELLS:
        rec = {"name": name, "arm": arm, "kp": kp, "ki": ki,
               "grid": [], "cfl_half": None}
        for nc in (120, 240, 480, 960):
            r = V2.sve_1d_v2(P_of(kp, ki), nc, 45.0, ctrl_on=arm, sensor_len_m=100.0)
            rec["grid"].append({
                "ncell": nc, "dx_m": round(r["dx_m"], 3),
                "dt_s": round(r["dt_s"], 4),
                "peak_mean_mm": round(r["peak_mean_mm"], 3),
                "peak_gate_mm": round(r["peak_gate_mm"], 3),
                "mass_balance_rel_err": r["mass_balance_rel_err"],
            })
        # CFL 减半（时间步收敛）：取最细网格 960 格
        base = V2.sve_1d_v2(P_of(kp, ki), 960, 45.0, ctrl_on=arm, sensor_len_m=100.0)
        half = V2.sve_1d_v2(P_of(kp, ki), 960, 45.0, ctrl_on=arm, sensor_len_m=100.0, cfl=V2.CFL / 2.0)
        rec["cfl_half"] = {
            "base_cfl": V2.CFL, "base_peak_mm": round(base["peak_mean_mm"], 3),
            "half_cfl": round(V2.CFL / 2.0, 3), "half_peak_mm": round(half["peak_mean_mm"], 3),
            "delta_mm": round(half["peak_mean_mm"] - base["peak_mean_mm"], 3),
        }
        out["cells"].append(rec)
        print(f"[{name}] done, grid peaks:", [g['peak_mean_mm'] for g in rec['grid']],
              "cfl_half delta:", rec['cfl_half']['delta_mm'], 'mm')
    out["elapsed_s"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "raw", "bm01_review_gridconv.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("SAVED raw/bm01_review_gridconv.json  (%.1fs)" % out["elapsed_s"])

if __name__ == "__main__":
    run()
