# -*- coding: utf-8 -*-
"""
R-03 / R-04 补实验：延长仿真核查 + 采样尺度物理反事实（优化版）
================================================================
对应审阅清单两条 P0 硬伤。本版用 480 格做延长仿真（480 是论文"最细网格收敛值"口径，
且单 case 12h 仅 ~36s，可快速出全量数据），960 格仅作"数值病态放大"的 4h 补充观察。

R-04 表8延长仿真核查：
  对表8"收敛=否"的3行 + 1对照行，480格延长 4h→8h→12h，区分 ①收敛 / ②单调漂移 / ③振荡。
R-03 采样尺度物理反事实：
  固定480格，扫 ctrl_on∈{gate, gate_avg, mean}，看"反馈测点采样尺度"如何翻转判定码。

只用 Python 标准库，import bm01_sve1d_v2。
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
    "duration_s": 14400.0,
}
STEP_DQ = 45.0
SETTLE_TOL_M = P["settle_tol_m"]


def sec_s0(m):
    sec = V2.Section(P["pool_width_m"], m)
    return V2.section_normal_slope(sec, P["h_target"], P["q0"], V2.MANNING_N)


def run_one(nc, step_dq, ctrl_on, m, dur):
    p = dict(P); p["duration_s"] = dur
    s0 = None if m == 0.0 else sec_s0(m)
    t0 = time.time()
    r = V2.sve_1d_v2(p, nc, step_dq, ctrl_on=ctrl_on, m_slope=m, s0=s0)
    r["wall_s"] = round(time.time() - t0, 1)
    return r


def classify(r):
    final = r["final_mean_mm"] / 1000.0
    peak = r["peak_mean_mm"] / 1000.0
    if final <= SETTLE_TOL_M:
        return "①收敛"
    if peak > final * 1.2:
        return "③振荡"
    return "②单调漂移"


def main():
    out = {
        "purpose": "R03_sampling_scale_counterfactual + R04_table8_extended_replay",
        "step_dq": STEP_DQ, "band_normal_mm": P["band_normal_m"] * 1000,
        "settle_tol_m": SETTLE_TOL_M, "settle_window_s": P["settle_window_s"],
        "R04_extended": {}, "R03_sampling_scale": {}, "R04_960_4h_pathology": {},
    }

    # ---- R-04：480 格延长仿真（论文最细网格口径）----
    r04 = [
        ("m=0|gate@480",     480, "gate",     0.0),   # 表8判"否"
        ("m=1.5|gate@480",   480, "gate",     1.5),   # 表8判"否"
        ("m=1.5|gate_avg@480", 480, "gate_avg", 1.5), # 表8判"否"（疑点行）
        ("m=0|gate_avg@480", 480, "gate_avg", 0.0),   # 对照：表8判"是"
    ]
    for label, nc, chan, m in r04:
        print(f"[R-04] {label} 480格延长仿真…", flush=True)
        row = {"grid": nc, "ctrl_on": chan, "m_slope": m, "durations": {}}
        for dur in (14400.0, 28800.0, 43200.0):
            r = run_one(nc, STEP_DQ, chan, m, dur)
            row["durations"][f"{int(dur/3600)}h"] = {
                "early_peak_mm": round(r["early_peak_mm"], 2),
                "peak_mean_mm": round(r["peak_mean_mm"], 2),
                "final_mean_mm": round(r["final_mean_mm"], 2),
                "over_band_s": round(r["over_band_s"], 1),
                "class": classify(r), "wall_s": r["wall_s"],
            }
            print(f"    {int(dur/3600)}h: peak={r['peak_mean_mm']:.1f} "
                  f"final={r['final_mean_mm']:.1f} -> {classify(r)}", flush=True)
        out["R04_extended"][label] = row

    # ---- R-04 补充：960 格 4h 数值病态观察（只跑 4h，快）----
    for label, nc, chan, m in [
        ("m=0|gate@960", 960, "gate", 0.0),
        ("m=1.5|gate@960", 960, "gate", 1.5),
        ("m=1.5|gate_avg@960", 960, "gate_avg", 1.5),
    ]:
        print(f"[R-04b] {label} 960格 4h 病态观察…", flush=True)
        r = run_one(nc, STEP_DQ, chan, m, 14400.0)
        out["R04_960_4h_pathology"][label] = {
            "peak_mean_mm": round(r["peak_mean_mm"], 2),
            "final_mean_mm": round(r["final_mean_mm"], 2),
            "peak_gate_mm": round(r["peak_gate_mm"], 2),
            "over_band_s": round(r["over_band_s"], 1),
            "wall_s": r["wall_s"],
        }
        print(f"    peak={r['peak_mean_mm']:.1f} gate_peak={r['peak_gate_mm']:.1f} "
              f"final={r['final_mean_mm']:.1f}", flush=True)

    # ---- R-03：采样尺度物理反事实（480格，三档测点）----
    print("[R-03] 采样尺度反事实（480格）…", flush=True)
    for m in (0.0, 1.5):
        key = f"m={m}"
        out["R03_sampling_scale"][key] = {}
        for chan in ("gate", "gate_avg", "mean"):
            r = run_one(480, STEP_DQ, chan, m, 14400.0)
            band_mm = P["band_normal_m"] * 1000
            judge = "A1(带内)" if r["peak_mean_mm"] <= band_mm else "A3(越带)"
            out["R03_sampling_scale"][key][chan] = {
                "early_peak_mm": round(r["early_peak_mm"], 2),
                "peak_mean_mm": round(r["peak_mean_mm"], 2),
                "final_mean_mm": round(r["final_mean_mm"], 2),
                "over_band_s": round(r["over_band_s"], 1),
                "judgement": judge, "wall_s": r["wall_s"],
            }
            print(f"    {key} {chan}: peak={r['peak_mean_mm']:.1f} final={r['final_mean_mm']:.1f} -> {judge}",
                  flush=True)

    os.makedirs("raw", exist_ok=True)
    outpath = os.path.join("raw", "bm01_r03_r04_extended.json")
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\n已写入 {outpath}")


if __name__ == "__main__":
    main()
