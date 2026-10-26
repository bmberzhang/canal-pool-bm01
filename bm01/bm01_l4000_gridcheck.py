# -*- coding: utf-8 -*-
"""补 L=4000 m（0.66）的三网格核对：dx = 50/25/12.5 m（ncell = 80/160/320）。
目的：新增的 0.66 触发点与论文既有三个关键算例一样，判定在网格上不变。"""
import sys, os, json, time
sys.path.insert(0, r"D:/自主水网测试案例/BME-01_中线渠池_上游阶跃可吸收".replace("BME", "BM"))
import bm01_sve1d_v2 as V2

L, DQ, DUR = 4000.0, 45.0, 43200.0
BAND = 150.0


def P():
    return {"pool_len_m": L, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": DUR}


ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))

out = {}
for dx in (50.0, 25.0, 12.5):
    ncell = int(round(L / dx))
    rec = {"dx_m": dx, "ncell": ncell}
    for name, kw in ARMS:
        r = V2.sve_1d_v2(P(), ncell, DQ, m_slope=0.0, **kw)
        rec[name] = {"peak_mm": round(r["peak_mean_mm"], 2),
                     "verdict": "A1" if r["peak_mean_mm"] <= BAND else "A3"}
        print("[L4000] dx=%5.1f n=%3d %-7s %8.2f mm %s" % (dx, ncell, name, rec[name]["peak_mm"], rec[name]["verdict"]), flush=True)
    out["dx%.1f" % dx] = rec

ok = all(out[k]["mean"]["verdict"] == "A1" and out[k]["point"]["verdict"] == "A3" for k in out)
print("\n三网格判定一致（mean A1 / point A3）:", ok)
json.dump(out, open("raw/bm01_l4000_gridcheck.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("[WRITE] raw/bm01_l4000_gridcheck.json")
