# -*- coding: utf-8 -*-
"""
BM-01 增益敏感性扫描 —— 0D 主模型口径（修正版）

上一版 bm01_gain_sensitivity.py 误用了 1D 圣维南求解器（sve_1d + gate_avg），
得到的是"一维矩形闸前局部越带"场景，与论文主案例（0D 蓄量模型单池+45m³/s 判 A1）口径不符。

本版直接复用 bm01_pool_step 的 0D 主模型（Params/simulate/judge），
在 kp×ki 平面上扫描主案例（+45 m³/s），量化"判定码对控制器增益的稳健性"。

运行：python bm01_gain_sensitivity_0d.py
输出：raw/bm01_gain_sensitivity_0d.json
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, Scenario, simulate, judge, build_scenarios

KP_LIST = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
KI_LIST = [0.0002, 0.0004, 0.0008, 0.0016, 0.0032, 0.0064, 0.0128]

def run(kp, ki):
    p = Params()
    p.kp = kp
    p.ki = ki
    # 主案例：单池 + 本地 PID，+45 m³/s
    sc = build_scenarios(p)[0]   # scenarios[0] = main
    r = simulate(p, sc)
    j = judge(p, sc, r)
    m = r["metrics"] if "metrics" in r else {}
    # 提取峰值偏差与末段偏差
    peak = max(abs(x) for x in r["series"]["h1_dev"]) if "h1_dev" in r["series"] else None
    # judge 结果
    return {
        "kp": kp, "ki": ki,
        "code": j["code"],
        "peak_m": round(peak, 4) if peak is not None else None,
    }

def main():
    # 先取主案例 scenario 基准（kp/ki 默认）
    base_sc = build_scenarios(Params())[0]
    rows = []
    for kp in KP_LIST:
        for ki in KI_LIST:
            row = run(kp, ki)
            rows.append(row)
            print(f"kp={kp:5.2f} ki={ki:8.5f} -> {row['code']}  峰 {row['peak_m']*1000 if row['peak_m'] else 0:7.1f}mm")

    n_a1 = sum(1 for r in rows if r["code"] == "A1")
    print(f"\n共 {len(rows)} 组合：A1={n_a1} ({n_a1/len(rows)*100:.0f}%)  A3={len(rows)-n_a1}")

    base = next(r for r in rows if r["kp"] == 1.6 and r["ki"] == 0.0016)
    print(f"基准点 kp=1.6 ki=0.0016 -> {base['code']}（峰 {base['peak_m']*1000 if base['peak_m'] else 0:.1f}mm）")

    os.makedirs("raw", exist_ok=True)
    out = {"grid": {"kp": KP_LIST, "ki": KI_LIST, "step_dq": base_sc.step_dq},
           "results": rows,
           "summary": {"n_total": len(rows), "n_a1": n_a1, "n_a3": len(rows)-n_a1,
                       "base_code": base["code"], "base_peak_m": base["peak_m"]}}
    with open("raw/bm01_gain_sensitivity_0d.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n已写入 raw/bm01_gain_sensitivity_0d.json")

if __name__ == "__main__":
    main()
