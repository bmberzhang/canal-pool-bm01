# -*- coding: utf-8 -*-
"""A4 机理时程图：BM-01 主算例 gate 反馈 vs mean 反馈的完整时序对比。
提取：闸门开度 gate、下泄流量 q_out、池均 dev_mean、闸前 dev_gate、上游 dev_head、
容积平衡残差（pool mean 由体积平衡给出，此处以 dev_mean 为池均，残差用 dev_gate 与 dev_mean 差示意）。
输出 SVG 供 Chrome 渲染 PNG。"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "raw", "bm01_review_mech_ts.json")


def P():
    return {'pool_len_m': 6000.0, 'pool_width_m': 50.0, 'q0': 180.0, 'h_target': 7.0,
            'a0': 0.5, 'kg': 360.0, 'kp': 1.6, 'ki': 0.0016, 'ctrl_period_s': 300.0,
            'deadband_m': 0.01, 'meas_resolution_m': 0.001, 'rate_limit_per_s': 1.0/3000.0,
            'duration_s': 43200.0, 't_step_s': 900.0, 'band_normal_m': 0.15,
            'band_alarm_m': 0.30, 'settle_tol_m': 0.02, 'settle_window_s': 1800.0}


def run(ctrl_on, ncell=240):
    r = V2.sve_1d_v2(P(), ncell, 45.0, ctrl_on=ctrl_on)
    return r


def main():
    t0 = time.time()
    data = {}
    for arm in ("gate", "mean"):
        r = run(arm)
        s = r["series"]
        # 下采样到 15 min 分辨率（series 可能很密）
        data[arm] = {
            "t_h": [x/3600.0 for x in s["t"]],
            "dev_mean_mm": [x*1000.0 for x in s["dev_mean"]],
            "dev_gate_mm": [x*1000.0 for x in s["dev_gate"]],
            "dev_head_mm": [x*1000.0 for x in s["dev_head"]],
            "gate_a": s["gate"],
            "q_out": s["q_out"],
            "peak_mean_mm": round(r["peak_mean_mm"], 2),
            "peak_gate_mm": round(r["peak_gate_mm"], 2),
        }
        print("%s: peak_mean=%.1f peak_gate=%.1f  n=%d  (%.0fs)"
              % (arm, r["peak_mean_mm"], r["peak_gate_mm"], len(s["t"]), time.time()-t0), flush=True)
    json.dump(data, open(OUT, "w", encoding="utf-8"), ensure_ascii=False)
    print("saved", OUT)


if __name__ == "__main__":
    main()
