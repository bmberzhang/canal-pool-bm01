# -*- coding: utf-8 -*-
"""
BM-01 附加诊断：前馈到底该「放多少」——增益 × 时机 二维扫描
================================================================================
问题
----
`chain_feedforward` 的前馈是「流量等价」的：来水多 72，就把下游闸门按能多放
72 的量开出去（`ff_gain = 1.0`）。直觉上应该少放一点（比如只放 36），因为池子
本身还有点调蓄能力。

本脚本把 `ff_gain` 和 `ff_delay_s` 一起扫，看 |峰值| 和振荡量最小的点在哪。

口径
----
复用 `bm01_ff_delay_sweep.run` / `.stats`（与 `bm01_pool_step.simulate` 逐行等价），
参数、控制器、判定码全不动，只改前馈通道的增益与参考时刻。
不构成工程结论，只把「少放一点行不行」变成一个数字。

自检
----
`ff_gain = 1.0, delay = 0, source = gate1` 必须复现主案例已知的 66.7 mm；
否则说明 run/stats 被改坏了，结果不可信。

用法
----
  python bm01_ff_gain_delay_sweep.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, build_scenarios
from bm01_ff_delay_sweep import run, stats

KNOWN_PEAK_MM = 66.7          # bm01_chain_feedforward.json 的已知值，用作自检锚点


def wobble(p: Params, h2):
    """振荡指标：阶跃后的过零次数 + 总变差(mm)。

    总变差 = 相邻采样点水位变化绝对值之和，曲线越平滑越小；
    比「越带时长」更敏感，能抓住「没超限但来回晃」的工况。
    """
    k0 = int(p.t_step_s / p.dt)
    y = h2[k0:]
    tv = float(np.abs(np.diff(y)).sum()) * 1000.0
    s = np.sign(y)
    s = s[s != 0]
    cross = int(np.sum(s[1:] * s[:-1] < 0)) if len(s) > 1 else 0
    return cross, tv


def sweep(p: Params, sc, src: str, gains, delays):
    """返回 {(gain, delay): (mx, cross, tv, over, sw)}。"""
    out = {}
    for g in gains:
        p.ff_gain = g
        for d in delays:
            h2, a2 = run(p, sc, float(d), src)
            _, _, mx, over, _, sw = stats(p, h2, a2)
            cross, tv = wobble(p, h2)
            out[(g, d)] = (mx, cross, tv, over, sw)
    return out


def main() -> int:
    p = Params()
    sc = next(s for s in build_scenarios(p) if s.key == "chain_feedforward")

    # ---- 自检：锚点必须复现 ----
    base_gain, base_p, base_tau = p.ff_gain, p.tau_ff_s, p.tau_chain_s
    h2, a2 = run(p, sc, 0.0, "gate1")
    _, _, mx0, _, _, _ = stats(p, h2, a2)
    assert abs(mx0 - KNOWN_PEAK_MM) < 0.15, \
        "自检失败：gain=1.0/delay=0/gate1 应为 %.1f mm，实得 %.1f mm" % (KNOWN_PEAK_MM, mx0)
    print("自检 OK：gain=1.0 delay=0 源=Q_G1 → %.1f mm（已知 %.1f）" % (mx0, KNOWN_PEAK_MM))

    # ---- 无前馈基准 ----
    sc_off = next(s for s in build_scenarios(p) if s.key == "chain_local_only")
    h2o, a2o = run(p, sc_off, 0.0)
    _, _, mx_off, over_off, _, sw_off = stats(p, h2o, a2o)
    cr_off, tv_off = wobble(p, h2o)
    print("无前馈基准：|峰| %.1f mm  过零 %d  总变差 %.0f mm  越带 %.0f s"
          % (mx_off, cr_off, tv_off, over_off))

    gains = [round(0.2 + 0.1 * i, 1) for i in range(13)]      # 0.2 ~ 1.4
    delays = [0, 120, 240, 300, 360, 480, 600, 900]
    print("阶跃 +%.0f m3/s　tau_chain_s=%.0f s　tau_ff_s=%.0f s　"
          "正常条带 ±%.2f m" % (sc.step_dq, p.tau_chain_s, p.tau_ff_s, p.band_normal_m))

    for src, title in [("inflow", "源 = 上游来水（传感器放在源头）"),
                       ("gate1", "源 = 上游闸门出流 Q_G1（现状接线）")]:
        print()
        print("=" * 100)
        print("[%s]　格子 = |峰值| mm（阶跃 +%.0f m3/s）" % (title, sc.step_dq))
        print("=" * 100)
        res = sweep(p, sc, src, gains, delays)
        print("%8s" % "gain\\delay" + "".join("%9s" % d for d in delays))
        print("-" * 100)
        best = None
        for g in gains:
            row = ""
            for d in delays:
                mx = res[(g, d)][0]
                row += "%9.1f" % mx
                if best is None or mx < best[0]:
                    best = (mx, g, d)
            print("%8.1f" % g + row)
        print("-" * 100)
        print("最小 |峰值| = %.1f mm　在 gain=%.1f, delay=%.0f s" % best)

        # 最优点附近的振荡量
        mx, g, d = best
        cross, tv, over, sw = res[(g, d)][1], res[(g, d)][2], res[(g, d)][3], res[(g, d)][4]
        print("该点：过零 %d 次　总变差 %.0f mm　越带 %.0f s　闸门行程 %.4f"
              % (cross, tv, over, sw))

        # 一维细扫：最优时机下，峰值/振荡随 gain 的变化
        print()
        print("细扫（delay=%d s 固定）：gain 从 0.2 到 1.4" % d)
        print("%8s %10s %10s %10s %10s %10s" %
              ("gain", "正向峰", "负向峰", "|峰|", "过零", "总变差"))
        print("-" * 100)
        p.ff_gain = 1.0
        for g in [round(0.2 + 0.05 * i, 2) for i in range(25)]:
            p.ff_gain = g
            h2, a2 = run(p, sc, float(d), src)
            pos, neg, mx, over, _, _ = stats(p, h2, a2)
            cross, tv = wobble(p, h2)
            print("%8.2f %10.1f %10.1f %10.1f %10d %10.0f" % (g, pos, neg, mx, cross, tv))

    # ---- 粗表步长 0.1 会跳过最优点，细扫定位（现状接线：Q_G1, delay=0）----
    print()
    print("=" * 100)
    print("[细扫定位] 源 = Q_G1, delay = 0：gain 0.70 ~ 1.00 步长 0.02")
    print("=" * 100)
    print("%8s %10s %10s %10s %8s %10s %9s %8s" %
          ("gain", "正向峰", "负向峰", "|峰|", "过零", "总变差", "越带(s)", "J"))
    print("-" * 100)
    fine = {}
    for i in range(16):
        g = round(0.70 + 0.02 * i, 2)
        p.ff_gain = g
        h2, a2 = run(p, sc, 0.0, "gate1")
        pos, neg, mx, over, _, _ = stats(p, h2, a2)
        cross, tv = wobble(p, h2)
        fine[g] = (pos, neg, mx, cross, tv, over)
    # 合并指标 J = 0.5·(|峰|/基准) + 0.5·(总变差/基准)，基准取现状 gain=1.0
    r_pk = fine[1.00][2]
    r_tv = fine[1.00][4]
    for g in fine:
        pos, neg, mx, cross, tv, over = fine[g]
        J = 0.5 * mx / r_pk + 0.5 * tv / r_tv
        print("%8.2f %10.1f %10.1f %10.1f %8d %10.0f %9.0f %8.3f" %
              (g, pos, neg, mx, cross, tv, over, J))
    print("-" * 100)
    g_pk = min(fine, key=lambda k: fine[k][2])
    g_tv = min(fine, key=lambda k: fine[k][4])
    g_bal = min(fine, key=lambda k: abs(fine[k][0] + fine[k][1]))   # 正负峰平衡
    g_J = min(fine, key=lambda k: 0.5 * fine[k][2] / r_pk + 0.5 * fine[k][4] / r_tv)
    print("|峰值| 最小：gain=%.2f → %.1f mm（现状 gain=1.0 为 %.1f mm，降幅 %.0f%%）"
          % (g_pk, fine[g_pk][2], fine[1.00][2],
             100 * (1 - fine[g_pk][2] / fine[1.00][2])))
    print("总变差最小：gain=%.2f → %.0f mm（现状 %.0f mm，降幅 %.0f%%）"
          % (g_tv, fine[g_tv][4], fine[1.00][4],
             100 * (1 - fine[g_tv][4] / fine[1.00][4])))
    print("正负峰最平衡：gain=%.2f（正 %.1f / 负 %.1f）"
          % (g_bal, fine[g_bal][0], fine[g_bal][1]))
    print("合并指标最小：gain=%.2f（|峰| %.1f mm，总变差 %.0f mm）"
          % (g_J, fine[g_J][2], fine[g_J][4]))
    print("推荐工作点 gain=%.2f：|峰| %.1f mm，过零 %d 次，总变差 %.0f mm，越带 %.0f s"
          % (g_J, fine[g_J][2], fine[g_J][3], fine[g_J][4], fine[g_J][5]))
    for g_half in [0.50, 1.00]:                      # 用户提的「多放 36」与现状
        p.ff_gain = g_half
        h2h, a2h = run(p, sc, 0.0, "gate1")
        _, _, mxh, _, _, _ = stats(p, h2h, a2h)
        crh, tvh = wobble(p, h2h)
        print("对照 gain=%.2f（多放 %.0f m3/s）：|峰| %.1f mm，过零 %d 次，总变差 %.0f mm"
              % (g_half, g_half * sc.step_dq, mxh, crh, tvh))

    # ---- 最优点是否随阶跃幅度漂移（决定能不能直接抄这个数）----
    from dataclasses import replace
    print()
    print("=" * 100)
    print("[敏感性] 阶跃幅度变化时，最优 gain 会不会漂")
    print("=" * 100)
    print("%10s %14s %10s %14s %10s" %
          ("阶跃(m3/s)", "最优gain", "|峰|最优", "gain=1.0时|峰|", "降幅"))
    print("-" * 100)
    for dq in [36.0, 45.0, 72.0, 120.0, 200.0]:
        sc_i = replace(sc, step_dq=dq)
        p.ff_gain = 1.0
        h2b, a2b = run(p, sc_i, 0.0, "gate1")
        _, _, mx_1, _, _, _ = stats(p, h2b, a2b)
        sc_b = None
        for i in range(21):
            g = round(0.60 + 0.02 * i, 2)
            p.ff_gain = g
            h2, a2 = run(p, sc_i, 0.0, "gate1")
            _, _, mx, _, _, _ = stats(p, h2, a2)
            if sc_b is None or mx < sc_b[1]:
                sc_b = (g, mx)
        print("%10.0f %14.2f %10.1f %14.1f %9.0f%%" %
              (dq, sc_b[0], sc_b[1], mx_1, 100 * (1 - sc_b[1] / mx_1)))

    p.ff_gain, p.tau_ff_s, p.tau_chain_s = base_gain, base_p, base_tau

    # ---- 最优增益对传播滞后 τ 的依赖：前馈在「水还没到」的窗口里抽空池子，
    #      抽空量 ∝ ff_gain·ΔQ·τ/A。若成立，τ 越大最优增益应越小。----
    print()
    print("=" * 100)
    print("[敏感性] 渠段传播滞后 tau_chain_s 变化时，最优 gain 会不会漂")
    print("=" * 100)
    print("%12s %12s %12s %12s %12s" %
          ("tau_chain(s)", "最优gain", "|峰|最优", "gain=1.0时|峰|", "降幅"))
    print("-" * 100)
    for tau in [150.0, 300.0, 450.0, 600.0, 900.0]:
        p.tau_chain_s = tau
        p.ff_gain = 1.0
        h2b, a2b = run(p, sc, 0.0, "gate1")
        _, _, mx_1, _, _, _ = stats(p, h2b, a2b)
        best_i = None
        for i in range(21):
            g = round(0.50 + 0.05 * i, 2)
            p.ff_gain = g
            h2, a2 = run(p, sc, 0.0, "gate1")
            _, _, mx, _, _, _ = stats(p, h2, a2)
            if best_i is None or mx < best_i[1]:
                best_i = (g, mx)
        print("%12.0f %12.2f %12.1f %12.1f %11.0f%%" %
              (tau, best_i[0], best_i[1], mx_1, 100 * (1 - best_i[1] / mx_1)))
    p.tau_chain_s = base_tau

    print()
    print("注：等效流量 = ff_gain × 阶跃幅度。gain=0.5 → 按「多放 36」设计；"
          "gain=%.2f → 按「多放 %.0f」设计（阶跃 %.0f 时）。"
          % (g_pk, g_pk * sc.step_dq, sc.step_dq))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
