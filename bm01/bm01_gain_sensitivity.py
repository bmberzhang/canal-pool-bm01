# -*- coding: utf-8 -*-
"""
BM-01 增益敏感性扫描（补审阅"漏点10"：kp=1.6/ki=0.0016 无整定依据）

用 bm01_sve1d_core.sve_1d 核心（纯 Python，快）在 kp×ki 平面上扫描，
输出每个 (kp,ki) 组合下的峰值偏差(mm)、越带时长(s)、末段收敛(mm)、判定码 A1/A3，
从而量化"判定码对控制器增益的稳健性"。

运行：
  python bm01_gain_sensitivity.py
输出：
  raw/bm01_gain_sensitivity.json  +  stdout 摘要

作者注：本脚本是论文投稿前的补充实验，非原始 BM-01 主结果。
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_sve1d_core import sve_1d

# 基准参数（与 bm01_sve1d_core.py __main__ 一致）
P = {
    "pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0,
    "h_target": 7.0, "a0": 0.5, "kg": 360.0,
    "kp": 1.6, "ki": 0.0016,
    "ctrl_period_s": 300.0, "deadband_m": 0.01, "meas_resolution_m": 0.001,
    "rate_limit_per_s": 1.0 / 3000.0, "duration_s": 14400.0,
    "t_step_s": 900.0, "band_normal_m": 0.15, "band_alarm_m": 0.3,
    "settle_tol_m": 0.02, "settle_window_s": 1800.0,
}

DQ = 45.0          # 主案例阶跃 +45 m³/s
NC = 480           # 最细网格
BAND_MM = 150.0    # ±150mm 条带
SETTLE_MM = 20.0   # R3 收敛阈值

# 增益网格（对数间隔，覆盖基准点 1.6 / 0.0016 附近）
KP_LIST = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
KI_LIST = [0.0002, 0.0004, 0.0008, 0.0016, 0.0032, 0.0064, 0.0128]

def judge(peak_mm, over_s, final_mm):
    r1 = peak_mm <= BAND_MM
    r2 = over_s <= 60.0          # R2 饱和（此处用越带时长近似饱和，保守）
    r3 = final_mm <= SETTLE_MM
    return "A1" if (r1 and r2 and r3) else "A3", r1, r2, r3

def run(kp, ki):
    p = dict(P); p["kp"] = kp; p["ki"] = ki
    r = sve_1d(p, NC, DQ, ctrl_on="gate_avg")
    code, r1, r2, r3 = judge(r["peak_mean_mm"], r["over_band_s"], r["final_mean_mm"])
    return {
        "kp": kp, "ki": ki,
        "peak_mean_mm": round(r["peak_mean_mm"], 2),
        "peak_gate_mm": round(r["peak_gate_mm"], 2),
        "over_band_s": round(r["over_band_s"], 1),
        "final_mean_mm": round(r["final_mean_mm"], 2),
        "mass_balance_rel_err": r["mass_balance_rel_err"],
        "code": code, "R1": r1, "R2": r2, "R3": r3,
    }

def main():
    rows = []
    for kp in KP_LIST:
        for ki in KI_LIST:
            row = run(kp, ki)
            rows.append(row)
            print(f"kp={kp:5.2f} ki={ki:8.5f} -> {row['code']}  "
                  f"峰{row['peak_mean_mm']:7.1f}mm 越带{row['over_band_s']:6.0f}s "
                  f"末段{row['final_mean_mm']:5.1f}mm")

    n_a1 = sum(1 for r in rows if r["code"] == "A1")
    print(f"\n共 {len(rows)} 组合：A1={n_a1} ({n_a1/len(rows)*100:.0f}%)  "
          f"A3={len(rows)-n_a1} ({ (len(rows)-n_a1)/len(rows)*100:.0f}%)")

    # 基准点位置
    base = next(r for r in rows if r["kp"] == 1.6 and r["ki"] == 0.0016)
    print(f"基准点 kp=1.6 ki=0.0016 -> {base['code']}（峰 {base['peak_mean_mm']}mm）")

    os.makedirs("raw", exist_ok=True)
    out = {
        "grid": {"kp": KP_LIST, "ki": KI_LIST, "dq": DQ, "ncell": NC,
                 "band_mm": BAND_MM, "settle_mm": SETTLE_MM},
        "results": rows,
        "summary": {"n_total": len(rows), "n_a1": n_a1, "n_a3": len(rows)-n_a1,
                    "base_code": base["code"], "base_peak_mm": base["peak_mean_mm"]},
    }
    with open("raw/bm01_gain_sensitivity.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n已写入 raw/bm01_gain_sensitivity.json")

if __name__ == "__main__":
    main()
