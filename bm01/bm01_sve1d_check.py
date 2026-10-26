# -*- coding: utf-8 -*-
"""
BM-01 一维对照的稳健性检查 —— 闸前的持续摆动是物理的还是数值的？
=================================================================
主对照跑出来：一维模型的闸前局部水位在 4 h 末段仍在 ±200 mm 上下持续摆动，
而池均水位只摆 ±60 mm。在把它写成结论之前，必须排除三种"假象"：

  检查 1（开环）：控制器关掉（闸门固定 a0=0.5）。此时只有"来水大、水位单向上抬"，
                  不该有任何摆动。→ 用来证明格式本身不产生振荡。
                  注意：开环水位是**单调上升**的，所以必须**先去趋势**再看残差，
                  否则量到的是"漂移"而不是"振荡"。
  检查 2（CFL） ：CFL 0.40→0.10（时间步更小）。若是格式引起，结果会明显变。
  检查 3（网格） ：120 / 480 两套网格的振荡周期是否一致。

单位提醒：本模块内部 dev_* 都是**米**，输出一律换算成 mm。

用法：python bm01_sve1d_check.py
"""

from __future__ import annotations

import json
import math
import os

import bm01_sve1d_compare as M


def tail_stats(t, y_m, lo, hi):
    """末段统计（输入米，输出毫米）。先去线性趋势，再量残差，避免把漂移当振荡。"""
    idx = [i for i in range(len(y_m)) if lo <= t[i] <= hi]
    xs = [t[i] for i in idx]
    ys = [y_m[i] * 1000.0 for i in idx]          # → mm
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    k = sxy / sxx if sxx else 0.0
    b = my - k * mx
    res = [y - (k * x + b) for x, y in zip(xs, ys)]
    # 残差上穿零的周期
    ups = [xs[i] for i in range(1, n) if res[i - 1] < 0 <= res[i]]
    per = (sum(ups[i + 1] - ups[i] for i in range(len(ups) - 1)) / (len(ups) - 1)
           if len(ups) >= 3 else None)
    return {"resid_half_range_mm": round((max(res) - min(res)) / 2.0, 1),
            "resid_rms_mm": round((sum(r * r for r in res) / n) ** 0.5, 1),
            "drift_mm_per_h": round(k * 3600.0, 1),
            "abs_half_range_mm": round((max(ys) - min(ys)) / 2.0, 1),
            "period_s": round(per, 1) if per else None}


def main() -> int:
    case = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(case, "raw", "bm01_main.json"), encoding="utf-8") as fh:
        raw = json.load(fh)
    p, dq = raw["params"], float(raw["scenario"]["step_dq"])
    lo, hi = p["duration_s"] - 3600.0, p["duration_s"]      # 末段 60 min

    c = math.sqrt(M.G * p["h_target"])
    print(f"波速 c = {c:.2f} m/s；传播时间 L/c = {p['pool_len_m']/c/60:.1f} min；"
          f"控制周期 = {p['ctrl_period_s']:.0f} s；末段窗口 = 最后 60 min")
    print("（残差 = 去掉线性趋势后的摆动；漂移 = 趋势本身）\n")

    res = {}

    def show(label, r):
        per = f"{r['period_s']/60:.1f}" if r["period_s"] else "—"
        print(f"  {label:<34} 残差半幅 {r['resid_half_range_mm']:>7.1f} mm｜"
              f"残差RMS {r['resid_rms_mm']:>6.1f}｜漂移 {r['drift_mm_per_h']:>7.1f} mm/h｜"
              f"周期 {per:>5} min")

    print("检查 1：开环（控制器关闭，闸门固定 0.50）——格式本身会不会振？")
    for nc in (120, 240):
        r = M.sve_1d(p, nc, dq, ctrl_on="off")
        st_g = tail_stats(r["series"]["t"], r["series"]["dev_gate"], lo, hi)
        st_m = tail_stats(r["series"]["t"], r["series"]["dev_mean"], lo, hi)
        res[f"openloop_{nc}_gate"] = st_g
        res[f"openloop_{nc}_mean"] = st_m
        show(f"开环 {nc} 格 · 闸前", st_g)
        show(f"开环 {nc} 格 · 池均", st_m)
    print()

    print("检查 2：CFL 敏感性（闭环，反馈=闸前）")
    base = None
    for cfl in (0.40, 0.10):
        old = M.CFL
        M.CFL = cfl
        try:
            r = M.sve_1d(p, 120, dq, ctrl_on="gate")
        finally:
            M.CFL = old
        st = tail_stats(r["series"]["t"], r["series"]["dev_gate"], lo, hi)
        res[f"cfl{cfl}"] = st
        show(f"闭环 120 格 · 闸前（CFL={cfl}）", st)
        if cfl == 0.40:
            base = {"peak_mean_mm": r["peak_mean_mm"], "st": st}
        else:
            res["cfl_peak_shift_mm"] = abs(r["peak_mean_mm"] - base["peak_mean_mm"])
    print()

    print("检查 3：网格加密（闭环，反馈=闸前）——周期随网格变不变？")
    for nc in (120, 240, 480):
        r = M.sve_1d(p, nc, dq, ctrl_on="gate")
        s = r["series"]
        st_g = tail_stats(s["t"], s["dev_gate"], lo, hi)
        st_m = tail_stats(s["t"], s["dev_mean"], lo, hi)
        st_h = tail_stats(s["t"], s["dev_head"], lo, hi)
        res[f"closed_{nc}_gate"] = st_g
        res[f"closed_{nc}_mean"] = st_m
        res[f"closed_{nc}_head"] = st_h
        res[f"closed_{nc}_peak_mean_mm"] = r["peak_mean_mm"]
        show(f"闭环 {nc} 格 · 闸前", st_g)
        show(f"闭环 {nc} 格 · 池均", st_m)
        show(f"闭环 {nc} 格 · 池首", st_h)
    print()

    # ---- 判读 ----
    ol = max(res["openloop_120_gate"]["resid_rms_mm"], res["openloop_240_gate"]["resid_rms_mm"])
    cl = res["closed_480_gate"]["resid_rms_mm"]
    per = [res[f"closed_{nc}_gate"]["period_s"] for nc in (240, 480)]
    cfl_shift = res["cfl_peak_shift_mm"]
    peak120 = res["closed_120_peak_mean_mm"]

    v1 = ol < 0.15 * cl
    v2 = cfl_shift < 0.02 * peak120
    v3 = all(per) and abs(per[0] - per[1]) / max(per) < 0.10
    print("--- 判读 ---")
    print(f"检查 1 开环闸前残差 RMS {ol:.1f} mm vs 闭环 {cl:.1f} mm → "
          f"{'开环基本不动，摆动来自控制回路（不是格式）' if v1 else '开环也在摆，需怀疑格式'}")
    print(f"检查 2 CFL 0.40→0.10，池均峰值变化 {cfl_shift:.2f} mm（{cfl_shift/peak120*100:.2f}%）→ "
          f"{'对时间步不敏感，不是格式引起' if v2 else '对时间步敏感'}")
    if all(per):
        print(f"检查 3 闸前摆动周期：240 格 {per[0]/60:.1f} min vs 480 格 {per[1]/60:.1f} min → "
              f"{'两网格一致，属物理（控制回路）周期' if v3 else '两网格不一致'}")
    else:
        print("检查 3 周期未能稳定识别，需人工看图")
    print(f"空间分布：闸前残差 {res['closed_480_gate']['resid_rms_mm']:.1f} mm、"
          f"池均 {res['closed_480_mean']['resid_rms_mm']:.1f} mm、"
          f"池首 {res['closed_480_head']['resid_rms_mm']:.1f} mm → "
          f"{'全池都在摆（不是闸门边界单元的局部伪影）' if res['closed_480_head']['resid_rms_mm'] > 0.3*res['closed_480_gate']['resid_rms_mm'] else '摆动集中在闸前，需怀疑边界处理'}")

    out = {"purpose": "robustness: is the gate oscillation physical (control loop) or numerical?",
           "units_note": "内部 dev_* 为米，本文件统计输出为毫米；残差已去线性趋势",
           "wave_celevity_mps": c, "travel_time_min": p["pool_len_m"] / c / 60,
           "ctrl_period_s": p["ctrl_period_s"], "tail_window_s": hi - lo,
           "results": res,
           "verdict": {"openloop_resid_rms_mm": ol, "closedloop_resid_rms_mm": cl,
                       "cfl_peak_shift_mm": cfl_shift,
                       "period_240_s": per[0], "period_480_s": per[1],
                       "checks_pass": {"numerics_clean": bool(v1),
                                       "dt_insensitive": bool(v2),
                                       "grid_independent_period": bool(v3)}},
           "note": "证据等级 synthetic-1d-sve；不改变 BM-01 判定码与证据等级。"}
    fp = os.path.join(case, "raw", "bm01_sve1d_robustness.json")
    with open(fp, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    print(f"\n输出：{fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
