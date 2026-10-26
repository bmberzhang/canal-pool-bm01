# -*- coding: utf-8 -*-
"""
BM-01 对照组：提前知道（前馈）vs 水到了才动（纯本地反馈）
================================================================================
问题
----
前面几轮扫的都是「提前知道」这一侧的参数（信号源、时机、增益）。本脚本补上对照组，
并回答一个更要紧的问题：**前馈的收益里，有多少是「信息提前」，有多少只是「不受
300 s 控制周期限制」？**

四组工况同对象、同阶跃（+72 m3/s）、同一条判定码：
  L1  纯本地反馈，控制周期 300 s（真实值）—— 水到了才动
  F1  上游信号实时前馈，gain=1.00（现状）
  F2  上游信号实时前馈，gain=0.82（推荐）
  L2  纯本地反馈，控制周期 5 s —— 同样是纯反馈，但不受慢控制周期限制
      → F1 与 L2 的差，才是「信息提前」真正值多少

时序口径（重要）
----------------
`Pool.control()` 里死区 `deadband_m` **只门控积分项**；比例项 `kp*e` 在每
一个控制周期都直接作用，只受测点分辨率 `meas_resolution_m = 1 mm` 限制。
所以本脚本报两个水位阈值：
  ≥1 mm  → 比例通道可以动作（量化分辨率）
  ≥10 mm → 死区，积分通道开始累积

口径
----
不改 BM-01 参数与判据；调用现成的 simulate()/judge()，只做对比与计时。

用法
----
  python bm01_ff_vs_local.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, build_scenarios, simulate, judge


def first_cross(t, y, thr, t0=0.0):
    """自 t0 起 |y| 首次超过 thr 的时刻(s)；从未超过返回 None。"""
    m = (np.abs(y) > thr) & (t >= t0)
    return float(t[int(np.argmax(m))]) if m.any() else None


def pool2_inflow(p: Params, q1_out, n):
    """复算池2 入流（串联一阶滞后状态），给「水什么时候到」一个物理数。

    ponytail: 与 simulate() 里的递推式重复；若串联逻辑改动需同步两处。
    """
    lag = np.empty(n)
    lag[0] = p.q0
    beta = p.dt / max(p.tau_chain_s, p.dt)
    for k in range(1, n):
        lag[k] = lag[k - 1] + beta * (q1_out[k - 1] - lag[k - 1])
    return lag


def main() -> int:
    p0 = Params()
    sc_ff = next(s for s in build_scenarios(p0) if s.key == "chain_feedforward")
    sc_lo = next(s for s in build_scenarios(p0) if s.key == "chain_local_only")

    cases = [
        ("L1 水到了才动（纯反馈 300s）", sc_lo, 1.00, 300.0),
        ("F1 提前知道 gain=1.00（现状）", sc_ff, 1.00, 300.0),
        ("F2 提前知道 gain=0.82（推荐）", sc_ff, 0.82, 300.0),
        ("L2 纯反馈但控制周期 5s", sc_lo, 1.00, 5.0),
    ]

    p = Params()
    print("=" * 108)
    print("BM-01 对照组　阶跃 +%.0f m3/s @ t=%.0f s　池2 正常带 ±%.0f mm　"
          "死区 %.0f mm　tau_chain %.0f s　闸门速率上限 %.5f /s"
          % (sc_ff.step_dq, p.t_step_s, p.band_normal_m * 1000, p.deadband_m * 1000,
             p.tau_chain_s, p.rate_limit_per_s))
    print("=" * 108)

    rows = []
    for name, sc, g, tc in cases:
        p.ff_gain = g
        p.ctrl_period_s = tc
        r = simulate(p, sc)
        j = judge(p, sc, r)
        s = r["series"]
        t, h2, a2, n = s["t"], s["h2_dev"], s["a2"], len(s["t"])
        lag = pool2_inflow(p, s["q1_out"], n)
        m50 = lag >= p.q0 + 0.5 * sc.step_dq
        rows.append(dict(
            name=name, sc=sc, r=r, j=j, t=t, h2=h2, tc=tc,
            t_arrive=first_cross(t, lag - p.q0, 1e-9),
            t_half=float(t[int(np.argmax(m50))]) if m50.any() else None,
            t_p=first_cross(t, h2, p.meas_resolution_m, p.t_step_s),
            t_i=first_cross(t, h2, p.deadband_m, p.t_step_s),
            t_act=first_cross(t, a2 - p.a0, 1e-6, p.t_step_s),
            travel=float(np.abs(np.diff(a2)).sum()),
        ))
    p.ctrl_period_s = p0.ctrl_period_s
    p.ff_gain = 1.0

    # ---- 时序 ----
    ts = p0.t_step_s
    print("\n[时序] 相对阶跃施加时刻 t=%.0f s 的延迟（s）" % ts)
    print("%-30s %11s %11s %11s %11s %11s" %
          ("工况", "水到(入流变)", "入流到50%", "水位≥1mm", "水位≥10mm", "闸门动"))
    print("-" * 108)
    for d in rows:
        rel = lambda x: "—" if x is None else "%+.0f" % (x - ts)
        print("%-30s %11s %11s %11s %11s %11s" %
              (d["name"], rel(d["t_arrive"]), rel(d["t_half"]),
               rel(d["t_p"]), rel(d["t_i"]), rel(d["t_act"])))
    print("-" * 108)

    # ---- 结果 ----
    print("\n[结果] 下游池 pool2")
    print("%-30s %9s %9s %9s %10s %10s %6s" %
          ("工况", "|峰|mm", "峰时刻s", "越带s", "闸门行程", "末段mm", "判定"))
    print("-" * 108)
    for d in rows:
        m2 = d["r"]["metrics"]["pool2"]
        pk = float(np.abs(d["h2"]).max()) * 1000
        tpk = float(d["t"][int(np.argmax(np.abs(d["h2"])))])
        print("%-30s %9.1f %9.0f %9.0f %10.4f %10.1f %6s" %
              (d["name"], pk, tpk, m2["time_over_normal_band_s"], d["travel"],
               m2["final_window_max_dev_m"] * 1000, d["j"]["code"]))
    print("-" * 108)

    print("\n[对照] 上游池 pool1——前馈只作用在下游闸门，本列应完全一致")
    for d in rows:
        m1 = d["r"]["metrics"]["pool1"]
        print("  %-30s |峰| %.1f mm　越带 %.0f s" %
              (d["name"], m1["peak_dev_m"] * 1000, m1["time_over_normal_band_s"]))

    # ---- 收益分解 ----
    pk = lambda d: float(np.abs(d["h2"]).max()) * 1000
    L1, F1, F2, L2 = rows
    ov = lambda d: d["r"]["metrics"]["pool2"]["time_over_normal_band_s"]
    print()
    print("[收益分解]")
    print("  干净单变量对照（唯一变量 = 有无前馈；两侧控制周期均 300 s）：")
    print("    L1 水到了才动 → F1 提前知道：%.1f → %.1f mm，降 %.0f%%；越带 %.0f → %.0f s；"
          "判定 %s → %s" % (pk(L1), pk(F1), 100 * (1 - pk(F1) / pk(L1)),
                            ov(L1), ov(F1), L1["j"]["code"], F1["j"]["code"]))
    print("    再整定增益（F2, gain=0.82）：%.1f mm，较纯反馈降 %.0f%%"
          % (pk(F2), 100 * (1 - pk(F2) / pk(L1))))
    print("  前馈的动作量反而更小：闸门行程 %.4f → %.4f（F1）→ %.4f（F2）"
          % (L1["travel"], F1["travel"], F2["travel"]))
    print("  前馈不靠「调得更猛」赢，靠「调得更准」。")
    print()
    print("  非对照参考（L2 纯反馈 + 控制周期 5 s）：%.1f mm" % pk(L2))
    print("    ⚠ 不可归因：ctrl_period_s 同时作用于两个池，L2 里上游池 pool1 峰值也从 "
          "%.1f 降到 %.1f mm，" % (L1["r"]["metrics"]["pool1"]["peak_dev_m"] * 1000,
                                   L2["r"]["metrics"]["pool1"]["peak_dev_m"] * 1000))
    print("      即「传到池2的扰动本身」已经变了，因此不能把收益拆成「快」与「信息」两块。")
    print("      要干净拆分需让 pool1 保持 300 s、只加快 pool2 —— 当前 simulate() 无此参数。")
    print()
    print("  闸门首次动作时刻：纯反馈 %+.0f s　前馈 %+.0f s　（早 %.0f s，约一个控制周期）"
          % (L1["t_act"] - ts, F1["t_act"] - ts, L1["t_act"] - F1["t_act"]))

    # ---- 2×2 析因：把「信息提前」与「控制快」真正拆开 ----
    # 关键：只给 pool2 换控制周期，上游池固定 300 s。直接改 p.ctrl_period_s 会同时
    #       加快两个池，传到下游的扰动本身也变了 —— 那就不是单变量对照了。
    from bm01_ff_delay_sweep import run

    tc_fast = 5.0
    fac = [
        ("A 无前馈 · pool2 慢 300s", sc_lo, False, 1.00, 300.0),
        ("B 无前馈 · pool2 快 5s", sc_lo, False, 1.00, tc_fast),
        ("C 前馈 · pool2 慢 300s", sc_ff, True, 1.00, 300.0),
        ("D 前馈 · pool2 快 5s", sc_ff, True, 1.00, tc_fast),
        ("E 前馈0.82 · pool2 快 5s", sc_ff, True, 0.82, tc_fast),
    ]
    print()
    print("=" * 108)
    print("[2×2 析因] 同一对象、同一阶跃；只改 pool2 控制周期，上游池固定 300 s")
    print("=" * 108)
    print("%-26s %9s %9s %11s %11s" % ("工况", "|峰|mm", "越带s", "闸门行程", "pool1|峰|"))
    print("-" * 108)
    out, up_ref = [], None
    for name, s_, use_ff, g, tc2 in fac:
        p.ff_gain = g if use_ff else 1.0
        tr = {}
        h2_, a2_ = run(p, s_, 0.0, "gate1", ctrl_period_p2=tc2, trace=tr)
        pk_ = float(np.abs(h2_).max()) * 1000
        ov_ = float(np.sum(np.abs(h2_) > p.band_normal_m)) * p.dt
        trv = float(np.abs(np.diff(a2_)).sum())
        up = float(np.abs(tr["h1"]).max()) * 1000
        # 自检：前馈只作用在下游闸门，上游池时序必须逐点相同
        if up_ref is None:
            up_ref = tr["h1"].copy()
        else:
            assert float(np.abs(up_ref - tr["h1"]).max()) == 0.0, \
                "上游池被改动，析因不成立：%s" % name
        out.append(pk_)
        print("%-26s %9.1f %9.0f %11.4f %11.1f" % (name, pk_, ov_, trv, up))
    p.ff_gain = 1.0
    print("-" * 108)
    print("自检 OK：五行 pool1 |峰| 完全一致 + pool1 时序逐点相同 → 上游扰动没变，析因成立。")

    A, B, C, D, E = out
    print()
    print("[收益归因] 以 A（水到了才动）%.1f mm 为基线" % A)
    print("  只把 pool2 控制周期加快（A→B）：%.1f mm，降 %.0f%%　←「快」的贡献"
          % (B, 100 * (1 - B / A)))
    print("  只加前馈、控制周期不动（A→C）：%.1f mm，降 %.0f%%　←「信息」的贡献"
          % (C, 100 * (1 - C / A)))
    print("  两者都上（A→D）：%.1f mm，降 %.0f%%" % (D, 100 * (1 - D / A)))
    print("  再加增益整定（A→E）：%.1f mm，降 %.0f%%" % (E, 100 * (1 - E / A)))
    print("  若两者独立可加，D 应 ≈ %.1f mm；实测 %.1f mm → 交互项 %+.1f mm"
          % (B + C - A, D, D - (B + C - A)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
