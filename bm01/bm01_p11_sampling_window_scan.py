# -*- coding: utf-8 -*-
"""
P1-1' 采样窗长连续扫描：区分"尺度效应(单调衰减)" vs "模态节点命中(阶跃式突降)"
==============================================================================
这是冲一区的机制鉴别关键实验。用 sve_1d_v2 的 gate_avg 档 + 不同 sensor_len_m，
对 m=0 / m=1.5 两个断面，480格，12h，扫描六档采样窗长：
  0(单点 gate) / 100 / 500 / 1500 / 3000 / 6000(全池 mean)
观察峰值 peak_mean_mm 随窗长的变化：单调下降→尺度效应；阶跃突降→模态节点命中。
"""
import json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

P = {
    "pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0,
    "h_target": 7.0, "a0": 0.5, "kg": 360.0, "kp": 1.6, "ki": 0.0016,
    "ctrl_period_s": 300.0, "deadband_m": 0.01, "meas_resolution_m": 0.001,
    "rate_limit_per_s": 1.0 / 3000.0,
    "t_step_s": 900.0, "band_normal_m": 0.15, "band_alarm_m": 0.3,
    "settle_tol_m": 0.02, "settle_window_s": 1800.0,
    "duration_s": 43200.0,   # 12h
}
STEP_DQ = 45.0
BAND_MM = 150.0

# 采样窗长档（m）：0=单点gate，6000=全池mean
WINDOWS = [0.0, 100.0, 500.0, 1500.0, 3000.0, 6000.0]


def sec_s0(m):
    sec = V2.Section(P["pool_width_m"], m)
    return V2.section_normal_slope(sec, P["h_target"], P["q0"], V2.MANNING_N)


def run(m, sensor_len):
    p = dict(P)
    s0 = None if m == 0.0 else sec_s0(m)
    t0 = time.time()
    if sensor_len <= 0.0:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="gate", m_slope=m, s0=s0)
        label = "gate(单点)"
    elif sensor_len >= 6000.0:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="mean", m_slope=m, s0=s0)
        label = "mean(全池)"
    else:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="gate_avg", m_slope=m, s0=s0,
                         sensor_len_m=sensor_len)
        label = f"gate_avg({sensor_len:.0f}m)"
    return label, r, time.time() - t0


def main():
    out = {"purpose": "P1-1_sampling_window_scan",
           "step_dq": STEP_DQ, "band_mm": BAND_MM, "duration": "12h",
           "grid": 480, "windows": WINDOWS, "results": {}}
    for m in (0.0, 1.5):
        key = f"m={m}"
        out["results"][key] = {}
        print(f"[P1-1'] m={m} 采样窗长连续扫描（480格 12h）…", flush=True)
        for w in WINDOWS:
            label, r, wall = run(m, w)
            judge = "A1(带内)" if r["peak_mean_mm"] <= BAND_MM else "A3(越带)"
            out["results"][key][f"{w:.0f}m"] = {
                "label": label,
                "early_peak_mm": round(r["early_peak_mm"], 2),
                "peak_mean_mm": round(r["peak_mean_mm"], 2),
                "final_mean_mm": round(r["final_mean_mm"], 2),
                "over_band_s": round(r["over_band_s"], 1),
                "judgement": judge, "wall_s": round(wall, 1),
            }
            print(f"    {w:>5.0f}m [{label}]: peak={r['peak_mean_mm']:>6.1f} "
                  f"final={r['final_mean_mm']:>6.1f} -> {judge}", flush=True)

    os.makedirs("raw", exist_ok=True)
    outpath = os.path.join("raw", "bm01_p11_sampling_window_scan.json")
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\n已写入 {outpath}")


if __name__ == "__main__":
    main()
