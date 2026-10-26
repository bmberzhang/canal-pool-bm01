# -*- coding: utf-8 -*-
"""
任务 G（Q&A 中被审稿方判定为"全文命门"）：非单调性的网格稳健性
================================================================================
质询点：若"中间窗长最不稳（argmax 在 1500m）"只是 480 格下的数值伪影，
则本文唯一的"偏离教科书描述函数"的现象就不成立 → 全文降级。

三网格 × 六窗长（m=0 矩形 / 12h / +45 m³/s）：
  grids   = 240, 480, 960
  windows = 0, 100, 500, 1500, 3000, 6000 (m)

判据（采纳 Q&A 第 5 轮 B 修订后的版本，替换其最初临时定的 50% 阈值）：
  G-b-1【形状判据·主】argmax 窗长所属档位在三个网格下必须同档：
        档位划分 = {0m} / {100,500,1500m} / {3000,6000m}
  G-b-2【幅度判据·辅】以 480 格为基准，960 与 240 格峰值相对差 ≤ 30%
        （Richardson 外推的"趋势可信"带）；超 30% 时正文禁止报绝对峰值，
        只报"同档 argmax 存在"。

输出：raw/bm01_g_grid_robustness.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

P = {"pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
     "kg": 360.0, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
     "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
     "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
     "settle_window_s": 1800.0, "duration_s": 43200.0}
GRIDS = [240, 480, 960]
WINDOWS = [0.0, 100.0, 500.0, 1500.0, 3000.0, 6000.0]
STEP_DQ = 45.0
BAND_MM = 150.0


def band_of(win):
    if win <= 0.5: return "A:{0m}"
    if win <= 1500.5: return "B:{100,500,1500m}"
    return "C:{3000,6000m}"


def run(grid, win):
    if win <= 0.0:
        r = V2.sve_1d_v2(P, grid, STEP_DQ, ctrl_on="gate", m_slope=0.0)
    elif win >= 6000.0:
        r = V2.sve_1d_v2(P, grid, STEP_DQ, ctrl_on="mean", m_slope=0.0)
    else:
        r = V2.sve_1d_v2(P, grid, STEP_DQ, ctrl_on="gate_avg", m_slope=0.0, sensor_len_m=win)
    return r


def main():
    t0 = time.time()
    res = {}
    for g in GRIDS:
        res[str(g)] = {}
        for w in WINDOWS:
            t1 = time.time()
            r = run(g, w)
            res[str(g)][f"{w:.0f}m"] = {
                "peak_mean_mm": round(r["peak_mean_mm"], 2),
                "final_mean_mm": round(r["final_mean_mm"], 2),
                "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3",
                "wall_s": round(time.time() - t1, 1),
            }
            print(f"[G] grid={g} win={w:>5.0f}m: peak={r['peak_mean_mm']:>7.1f} "
                  f"final={r['final_mean_mm']:>7.1f} ({time.time()-t1:.0f}s)", flush=True)

    out = {"purpose": "G: 非单调性的网格稳健性", "grids": GRIDS, "windows": WINDOWS,
           "m_slope": 0.0, "duration_h": 12.0, "step_dq": STEP_DQ, "results": res,
           "band_partition": {"A": "{0m}", "B": "{100,500,1500m}", "C": "{3000,6000m}"}}

    # 判据
    argmax = {}
    for g in GRIDS:
        wbest = max(WINDOWS, key=lambda w: res[str(g)][f"{w:.0f}m"]["peak_mean_mm"])
        argmax[str(g)] = {"argmax_win_m": wbest, "argmax_band": band_of(wbest),
                          "argmax_peak_mm": res[str(g)][f"{wbest:.0f}m"]["peak_mean_mm"]}
    out["argmax"] = argmax
    bands = {argmax[str(g)]["argmax_band"] for g in GRIDS}
    out["G_b1_shape_pass"] = len(bands) == 1
    base = {w: res["480"][f"{w:.0f}m"]["peak_mean_mm"] for w in WINDOWS}
    rel = {}
    for g in ("240", "960"):
        rel[g] = {f"{w:.0f}m": round(abs(res[g][f"{w:.0f}m"]["peak_mean_mm"] - base[w]) / max(base[w], 1e-9) * 100, 2)
                  for w in WINDOWS}
    out["G_b2_rel_diff_pct_vs_480"] = rel
    out["G_b2_amplitude_pass"] = all(v <= 30.0 for g in rel for v in rel[g].values())
    out["verdict"] = ("G 通过（形状同档 + 幅度≤30%）" if (out["G_b1_shape_pass"] and out["G_b2_amplitude_pass"])
                      else "G 部分通过（形状同档，幅度超带）" if out["G_b1_shape_pass"]
                      else "G 失败（argmax 跨档 → 非单调性疑为网格伪影）")
    out["wall_total_s"] = round(time.time() - t0, 1)

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_g_grid_robustness.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n已写入 {jp}（{out['wall_total_s']}s）")
    print(f"argmax: {argmax}")
    print(f"G-b1 形状通过={out['G_b1_shape_pass']} | G-b2 幅度通过={out['G_b2_amplitude_pass']} | {out['verdict']}")


if __name__ == "__main__":
    main()
