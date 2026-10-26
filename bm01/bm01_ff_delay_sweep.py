# -*- coding: utf-8 -*-
"""
BM-01 附加诊断：前馈到底该「多早」动作才最优？
================================================================================
问题
----
`chain_feedforward` 里的前馈是**实时**的：

    下游闸门在 t 时刻，按**当前**上游闸门流量变化算前馈量。

但下游池在 t 时刻真正感受到的入流，是 **τ 之前**的上游出流
（串联用 `tau_chain_s = 300 s` 的一阶滞后描述）。

也就是说：前馈**早动了 τ**。水还没到，闸门先开大 → 池子被抽空（负向偏差）。
再等到反馈环发现水位掉了、把闸门关回来，振荡就起来了。

本脚本扫「前馈纯延迟」ff_delay，看 |峰值| 最小的点在哪。
预期：最优 ff_delay ≈ tau_chain_s，而不是现在的 0。

口径
----
与主案例同一套参数、同一套控制器（kp/ki/速率限/控制周期全不动），
**只改前馈信号的时序**。不改 BM-01 的判定码与证据等级，不构成工程结论。
本文件只为把「池首早知道入流」这件事到底值多少变成一个数字。

用法
----
  python bm01_ff_delay_sweep.py
"""

import math
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, Pool, build_scenarios


def run(p: Params, sc, ff_delay_s: float, ff_source: str = "gate1",
        ctrl_period_p2=None, trace=None):
    """与 bm01_pool_step.simulate 逐行等价，只把前馈参考量整体后移 ff_delay_s。

    ff_source = "gate1"  前馈源取上游闸门出流 Q_G1（现状；它本身滞后于来水阶跃）
              = "inflow" 前馈源取上游来水 q_up（相当于把传感器一直放到源头）

    ctrl_period_p2  只改下游池(pool2)的控制周期，上游池保持 p.ctrl_period_s。
                    `Params` 是全局的，直接改 p.ctrl_period_s 会同时加快两个池，
                    于是传到下游的扰动本身也变了 —— 那就不是单变量对照了。
                    默认 None = 逐行等价原行为。
    trace           传入 dict 时回填上游池时序 {"h1":..., "a1":...}，用于验证
                    「上游池确实没被改动」。
    """
    n = int(round(p.duration_s / p.dt)) + 1
    t = np.arange(n) * p.dt

    pools = [Pool(p) for _ in range(sc.pools)]
    if ctrl_period_p2 is not None and sc.pools == 2:
        # Pool 持有的是 Params 引用；给 pool2 单换一份改了控制周期的副本，其余字段全同
        pools[1].p = replace(p, ctrl_period_s=ctrl_period_p2)
    need_up = trace is not None
    h1 = np.zeros(n) if need_up else None
    a1 = np.zeros(n) if need_up else None
    lag_state = p.q0
    ff_state = 0.0
    h2 = np.zeros(n)
    a2 = np.zeros(n)
    q1o = np.zeros(n)
    nd = int(round(ff_delay_s / p.dt))

    for k in range(n):
        tk = float(t[k])
        q_up = p.q0 + (sc.step_dq if tk >= p.t_step_s else 0.0)
        pools[0].q_in = q_up

        if sc.pools == 2:
            if ff_source == "gate1":
                # 基线必须与 bm01_pool_step.simulate 一致：它读的是上一时刻的 q1_out。
                # q1o[k] 在本轮循环末尾才赋值，直接读 k 会拿到 0（初始化值），
                # 于是 ff_target=(0-180)/360=-0.5，闸门被关死，得到一个假峰。
                kd = k - 1 - nd
                q1_ref = q1o[kd] if kd >= 0 else p.q0
            else:
                # 源头信号是解析可算的，直接回看 nd 步
                q1_ref = p.q0 + (sc.step_dq if (tk - ff_delay_s) >= p.t_step_s else 0.0)
            ff_target = (p.ff_gain * (q1_ref - p.q0) / p.kg) if sc.feedforward else 0.0
            alpha = p.dt / max(p.tau_ff_s, p.dt)
            ff_state += alpha * (ff_target - ff_state)

        if sc.pools == 1:
            pools[0].control(tk)
        else:
            pools[0].control(tk)
            pools[1].control(tk, feedforward_open=ff_state)

        if sc.pools == 2:
            beta = p.dt / max(p.tau_chain_s, p.dt)
            lag_state += beta * (q1o[k - 1] - lag_state) if k > 0 else 0.0
            pools[1].q_in = lag_state

        q1 = pools[0].dynamics(p.dt)
        q1o[k] = q1
        for pl in pools:
            pl.step_gate(p.dt)
            pl.account_saturation(p.dt)
        if sc.pools == 2:
            pools[1].dynamics(p.dt)
            h2[k] = pools[1].h - p.h_target
            a2[k] = pools[1].a
        if need_up:
            h1[k] = pools[0].h - p.h_target
            a1[k] = pools[0].a

    if need_up:
        trace["h1"], trace["a1"] = h1, a1
    return h2, a2


def stats(p: Params, h2, a2):
    pos = float(h2.max()) * 1000
    neg = float(h2.min()) * 1000
    mx = float(np.abs(h2).max()) * 1000
    over = float(np.sum(np.abs(h2) > p.band_normal_m)) * p.dt
    tail = float(np.abs(h2[-int(p.settle_window_s / p.dt):]).max()) * 1000
    sw = 0.0
    for i in range(1, len(a2)):
        sw += abs(a2[i] - a2[i - 1])
    return pos, neg, mx, over, tail, sw


def main() -> int:
    p = Params()
    sc = next(s for s in build_scenarios(p) if s.key == "chain_feedforward")

    print("=" * 92)
    print("BM-01 前馈时序扫描　对象：两池串联 + 渠段前馈　阶跃 +%.0f m3/s" % sc.step_dq)
    print("tau_chain_s = %.0f s　ff_gain = %.2f　tau_ff_s = %.0f s　"
          "正常条带 ±%.2f m" % (p.tau_chain_s, p.ff_gain, p.tau_ff_s, p.band_normal_m))
    print("=" * 92)

    for src, title in [("gate1", "源 = 上游闸门出流 Q_G1（现状）"),
                       ("inflow", "源 = 上游来水 q_up（传感器一直放到源头）")]:
        print()
        print("[%s]" % title)
        print("%10s %11s %11s %11s %11s %11s" %
              ("前馈延迟", "正向峰值", "负向峰值", "|最大偏差|", "越带时长", "闸门总行程"))
        print("%10s %11s %11s %11s %11s %11s" %
              ("(s)", "(mm)", "(mm)", "(mm)", "(s)", "(开度累计)"))
        print("-" * 92)
        best = None
        for d in [0, 60, 120, 180, 240, 300, 360, 420, 480, 600, 900, 1200]:
            h2, a2 = run(p, sc, float(d), src)
            pos, neg, mx, over, tail, sw = stats(p, h2, a2)
            mark = ""
            if best is None or mx < best[1]:
                best = (d, mx)
                mark = "  ←"
            print("%10d %11.1f %11.1f %11.1f %11.1f %11.3f%s" %
                  (d, pos, neg, mx, over, sw, mark))
        print("-" * 92)
        print("最优前馈延迟 ≈ %.0f s（|最大偏差| %.1f mm）" % (best[0], best[1]))

    # 基准：完全关掉前馈（对应 chain_local_only 的时序）
    sc_off = next(s for s in build_scenarios(p) if s.key == "chain_local_only")
    h2o, a2o = run(p, sc_off, 0.0)
    pos, neg, mx, over, tail, sw = stats(p, h2o, a2o)
    print()
    print("-" * 92)
    print("无前馈基准   %11.1f %11.1f %11.1f %11.1f %11.3f" % (pos, neg, mx, over, sw))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
