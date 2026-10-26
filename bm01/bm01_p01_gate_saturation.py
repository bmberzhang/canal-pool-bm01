# -*- coding: utf-8 -*-
"""
P0-1：闸门开度时程检查 —— H6 机制的生死判
==========================================
判断在"采样窗长扫描"中，闸门开度 a(t) 是否发生饱和（a→0 或 a→1）。
- 若饱和 → H6 机制改为"饱和极限环 + 窗口依赖"，机制变弱但仍有工程价值；
- 若不饱和 → H6 "相位延迟调制稳定性" 机制存活。

对每个采样窗长档，记录：
  a_min, a_max, 饱和时长占比（a<1e-9 或 a>1-1e-9 的时间占比）
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
WINDOWS = [0.0, 100.0, 500.0, 1500.0, 3000.0, 6000.0]


def sec_s0(m):
    sec = V2.Section(P["pool_width_m"], m)
    return V2.section_normal_slope(sec, P["h_target"], P["q0"], V2.MANNING_N)


def gate_saturation_analysis(m, sensor_len):
    """跑一个场景，分析闸门开度时程的饱和情况。"""
    p = dict(P)
    s0 = None if m == 0.0 else sec_s0(m)
    if sensor_len <= 0.0:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="gate", m_slope=m, s0=s0)
        label = "gate(0m)"
    elif sensor_len >= 6000.0:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="mean", m_slope=m, s0=s0)
        label = "mean(6000m)"
    else:
        r = V2.sve_1d_v2(p, 480, STEP_DQ, ctrl_on="gate_avg", m_slope=m, s0=s0,
                         sensor_len_m=sensor_len)
        label = f"gate_avg({sensor_len:.0f}m)"

    gate = r["series"]["gate"]
    n = len(gate)
    a_min = min(gate)
    a_max = max(gate)
    # 饱和时长占比（a 接近 0 或 1）
    sat_lo = sum(1 for a in gate if a <= 1e-9)
    sat_hi = sum(1 for a in gate if a >= 1.0 - 1e-9)
    sat_frac = (sat_lo + sat_hi) / n
    return label, r, a_min, a_max, sat_lo / n, sat_hi / n, sat_frac


def main():
    out = {"purpose": "P0-1_gate_saturation_check", "duration": "12h",
           "grid": 480, "results": {}}
    for m in (0.0, 1.5):
        key = f"m={m}"
        out["results"][key] = {}
        print(f"[P0-1] m={m} 闸门开度饱和检查…", flush=True)
        for w in WINDOWS:
            t0 = time.time()
            label, r, a_min, a_max, f_lo, f_hi, sat = gate_saturation_analysis(m, w)
            out["results"][key][f"{w:.0f}m"] = {
                "label": label, "a_min": round(a_min, 4), "a_max": round(a_max, 4),
                "sat_lo_frac": round(f_lo, 4), "sat_hi_frac": round(f_hi, 4),
                "sat_total_frac": round(sat, 4),
                "peak_mean_mm": round(r["peak_mean_mm"], 2),
                "wall_s": round(time.time() - t0, 1),
            }
            verdict = "饱和" if sat > 0.01 else "未饱和"
            print(f"    {w:>5.0f}m [{label}]: a∈[{a_min:.3f},{a_max:.3f}] "
                  f"饱和占比{sat*100:.1f}% -> {verdict} | peak={r['peak_mean_mm']:.1f}mm",
                  flush=True)

    os.makedirs("raw", exist_ok=True)
    outpath = os.path.join("raw", "bm01_p01_gate_saturation.json")
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\n已写入 {outpath}")


if __name__ == "__main__":
    main()
