# -*- coding: utf-8 -*-
"""核查 BM-01R 池均臂 mass_balance_rel_err = 1.56% 是"分母太小"还是真守恒缺陷。
输出：绝对误差、相对净入流、相对总过流（吞吐量）。"""
import sys, os, time
sys.path.insert(0, r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收")
import bm01_sve1d_v2 as SVE2

B_M, M_SLOPE, H_T, L_M = 20.75, 2.25, 4.50, 25567.0
Q0, DQ, NMAN, KG, NCELL = 135.0, 45.0, 0.0150, 270.0, 256


def P():
    return {"pool_width_m": B_M, "pool_len_m": L_M, "h_target": H_T, "q0": Q0, "kg": KG,
            "duration_s": 86400.0, "t_step_s": 900.0, "settle_window_s": 1800.0,
            "band_normal_m": 0.15, "a0": 0.50, "kp": 1.60, "ki": 1.60e-3,
            "ctrl_period_s": 300.0, "deadband_m": 0.010, "meas_resolution_m": 0.001,
            "rate_limit_per_s": 0.02 / 60.0, "band_alarm_m": 0.30, "settle_tol_m": 0.020}


sec = SVE2.Section(B_M, M_SLOPE)
s0 = SVE2.section_normal_slope(sec, H_T, Q0, NMAN)
for name, kw in (("mean", dict(ctrl_on="mean")), ("point", dict(ctrl_on="gate"))):
    t1 = time.time()
    r = SVE2.sve_1d_v2(P(), NCELL, DQ, m_slope=M_SLOPE, n=NMAN, s0=s0, **kw)
    dv, vn = r["vol_storage_m3"], r["vol_net_m3"]
    T = 86400.0
    gross = (Q0 + DQ) * T          # 总入流体积
    print("%-6s 池内蓄量增量 %12.1f m³  净入流 %12.1f m³  绝对差 %10.1f m³" % (name, dv, vn, dv - vn))
    print("       相对净入流 %8.4f%%   相对总过流(%.2e m³) %.2e   峰值 %.2f mm  (%.0fs)"
          % (abs(dv - vn) / max(abs(vn), 1) * 100, gross, abs(dv - vn) / gross,
             r["peak_mean_mm"], time.time() - t1))
