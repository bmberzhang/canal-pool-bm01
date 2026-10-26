# -*- coding: utf-8 -*-
"""
BM-01R 附录：长时仿真，验证"渠道输水能力"这条硬约束
================================================================================
短时（24 h）仿真显示 B1/B2 的池均水位偏差仍在单调上涨、远未收敛。
本脚本把真实渠段工况跑到 10 天，检查水位是否收敛到 Manning 均匀流水深。

判据（这是一条**独立于控制器**的物理约束）：
  渠底坡度固定时，通过的流量越大，所需水深越大。
  真实渠段 S0 = 4.435e-5 下：

      Q = 135 m³/s → h_n = 4.500 m （设计水深，S0 就是这么反解的）
      Q = 180 m³/s → h_n = 5.242 m （+742 mm）

  闸门能改变**过闸流量**，但改变不了**渠道本身的输水能力**。
  0D 集总模型把过流能力全部放在闸门上（Q = kg·a·√(h/ht)），
  所以它看不到这条约束——这正是 0D 与 1D 分歧的物理根源。

用法
----
  python bm01r_converge.py
  python bm01r_converge.py --days 10 --ncell 128
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bm01r_lit as R

G = 9.81


def manning_q(h: float, b: float, m: float, n: float, S0: float) -> float:
    A = h * (b + m * h)
    P = b + 2.0 * h * math.sqrt(1.0 + m * m)
    return (1.0 / n) * A * (A / P) ** (2.0 / 3.0) * math.sqrt(S0)


def normal_depth(Q: float, b: float, m: float, n: float, S0: float) -> float:
    lo, hi = 0.05, 40.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if manning_q(mid, b, m, n, S0) < Q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01R 长时收敛验证")
    ap.add_argument("--days", type=float, default=10.0)
    ap.add_argument("--ncell", type=int, default=128)
    ap.add_argument("--out", default=HERE)
    args = ap.parse_args()

    c = R.make_lit()
    b, m, ht, n, q0 = c["bottom_m"], c["m_slope"], c["h_target"], c["n"], c["q0"]
    sec = R.SVE2.Section(b, m)
    S0 = R.SVE2.section_normal_slope(sec, ht, q0, n)
    dur = args.days * 86400.0

    print("=" * 78)
    print(f"BM-01R 附录：长时收敛验证（{args.days:.0f} 天，{args.ncell} 格）")
    print("=" * 78)

    h_design = normal_depth(q0, b, m, n, S0)
    h_new = normal_depth(q0 + c["step_dq"], b, m, n, S0)
    print(f"\n【硬约束】渠底坡度 S0 = {S0:.4e} 固定，Manning 均匀流水深：")
    print(f"  Q = {q0:.0f} m³/s → h_n = {h_design:.3f} m")
    print(f"  Q = {q0 + c['step_dq']:.0f} m³/s → h_n = {h_new:.3f} m"
          f"  （设计水深之上 +{(h_new - ht) * 1000:.0f} mm）")
    print(f"  水位条带 ±{c['band_normal_m'] * 1000:.0f} mm ⇒ 缺口是条带的 "
          f"{(h_new - ht) / c['band_normal_m']:.1f} 倍")

    print(f"\n【仿真】真实渠段 · 梯形 m={m} · +{c['step_dq']:.0f} m³/s · {args.days:.0f} 天")
    t0 = time.time()
    r = R.run1d(c, args.ncell, m_slope=m, duration_s=dur)
    print(f"  完成：{r['steps']} 步 × {r['ncell']} 格，dt={r['dt_s']} s，"
          f"耗时 {time.time() - t0:.0f} s")
    print(f"  质量守恒 {r['mass_balance_rel_err']:.1e}")

    s = r["series"]
    t, dm, dg = s["t"], s["dev_mean"], s["dev_gate"]
    print(f"\n【时序】每 12 h 采样（池均 / 闸前，单位 mm）")
    marks = []
    for step in range(int(args.days * 2) + 1):
        tq = 43200.0 * step
        i = min(range(len(t)), key=lambda k: abs(t[k] - tq)) if t else 0
        if tq > t[-1] + 1:
            break
        print(f"  t = {t[i] / 3600:7.1f} h　池均 {dm[i] * 1000:9.1f}　闸前 {dg[i] * 1000:9.1f}")
        marks.append({"t_h": round(t[i] / 3600, 1),
                      "dev_mean_mm": round(dm[i] * 1000, 1),
                      "dev_gate_mm": round(dg[i] * 1000, 1)})

    print(f"\n【收敛判据】")
    print(f"  池均偏差峰值 {r['peak_mean_mm']:.1f} mm｜闸前偏差峰值 {r['peak_gate_mm']:.1f} mm")
    print(f"  末段（最后 30 min）池均最大偏差 {r['final_mean_mm']:.1f} mm　"
          f"※ final_mean_mm 是末段窗口内的最大值，不是末段均值（字段名沿用求解器约定）")
    print(f"  越正常条带累计 {r['over_band_s']:.0f} s（全程 {dur:.0f} s）")
    trend = dm[-1] * 1000 - dm[len(dm) // 2] * 1000
    print(f"  后半程池均变化 {trend:+.1f} mm ⇒ "
          f"{'已趋稳' if abs(trend) < 20 else '仍在漂移，说明平衡尚未到达'}")

    out = {
        "appendix": "long-run convergence check",
        "days": args.days, "ncell": args.ncell, "dt_s": r["dt_s"],
        "S0": S0, "manning": {
            "q_design_m3s": q0, "h_normal_design_m": round(h_design, 4),
            "q_stepped_m3s": q0 + c["step_dq"], "h_normal_stepped_m": round(h_new, 4),
            "rise_required_mm": round((h_new - ht) * 1000, 1),
            "band_mm": c["band_normal_m"] * 1000,
            "ratio_to_band": round((h_new - ht) / c["band_normal_m"], 2),
        },
        "metrics": {
            "peak_mean_mm": r["peak_mean_mm"], "peak_gate_mm": r["peak_gate_mm"],
            "final_mean_mm": r["final_mean_mm"], "over_band_s": r["over_band_s"],
            "mass_balance_rel_err": r["mass_balance_rel_err"],
        },
        "samples": marks,
        "note": ("本脚本只做收敛性检验，不产生判定码。判定见 bm01r_lit.py。"
                 "参数为公开文献值（literature-parameterized），sealed: false。"),
    }
    p = os.path.join(os.path.abspath(args.out), "raw", "bm01r_convergence.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
