# -*- coding: utf-8 -*-
"""
补 6.2 触发点表的中间档：L = 4000 m（T_c/(L/c) = 300/454.5 ≈ 0.66）
=====================================================================
审稿风险："0.44 翻转、1.32 崩坏，那 0.6 / 0.8 呢？" —— 表里只有 0.22/0.44/1.32，
读者无法判断 0.66 落在哪一侧。本脚本按与 bm01_mc_multicase.py 完全相同的设置
（固定 dx = 25 m → ncell = 160，12 h，ΔQ = 45 m³/s，三臂 gate/avg100/mean）
补跑 L = 4000 m，并同时重跑 L = 6000 m 主算例做回归自检
（应复现 pool-mean 75.05 / gate 470.78 mm，否则本次结果不可用）。
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

DX_M = 25.0
BAND_MM = 150.0
DURATION_S = 43200.0
DQ = 45.0
C_CEL = 8.8           # 论文 6.2 用的对流波速 c


def base_p(L):
    return {"pool_len_m": L, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": DURATION_S}


def code(pk):
    return "A1" if pk <= BAND_MM else "A3"


ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))

out = {}
for L in (6000.0, 4000.0):
    ncell = int(round(L / DX_M))
    rec = {"pool_len_m": L, "ncell": ncell, "dx_m": L / ncell,
           "L_over_c_s": round(L / C_CEL, 1),
           "ratio": round(300.0 / (L / C_CEL), 2), "dq": DQ}
    for name, kw in ARMS:
        t1 = time.time()
        r = V2.sve_1d_v2(base_p(L), ncell, DQ, m_slope=0.0, **kw)
        rec[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                     "verdict": code(r["peak_mean_mm"]),
                     "mass_err": float("%.2e" % r["mass_balance_rel_err"]),
                     "wall_s": round(time.time() - t1, 1)}
        print("[TRIG] L=%5.0f ncell=%3d ratio=%.2f %-7s peak=%8.2f mm %s (%.0fs)"
              % (L, ncell, rec["ratio"], name, rec[name]["peak_mean_mm"],
                 rec[name]["verdict"], rec[name]["wall_s"]), flush=True)
    rec["flip_point_vs_mean"] = (rec["mean"]["verdict"] == "A1" and rec["point"]["verdict"] == "A3")
    rec["flip_avg100_vs_mean"] = (rec["mean"]["verdict"] == "A1" and rec["avg100"]["verdict"] == "A3")
    out["L%d" % int(L)] = rec
    print("[TRIG] flip(point vs mean) =", rec["flip_point_vs_mean"], flush=True)

# 回归自检：主算例必须复现
chk = out["L6000"]
ok = abs(chk["mean"]["peak_mean_mm"] - 75.05) < 0.5 and abs(chk["point"]["peak_mean_mm"] - 470.78) < 3.0
print("\n[REGRESSION] L=6000 reproduce: mean=%.2f (expect 75.05), gate=%.2f (expect 470.78) -> %s"
      % (chk["mean"]["peak_mean_mm"], chk["point"]["peak_mean_mm"], "PASS" if ok else "FAIL"), flush=True)

os.makedirs("raw", exist_ok=True)
with open("raw/bm01_l4000_trigger.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("[WRITE] raw/bm01_l4000_trigger.json", flush=True)
