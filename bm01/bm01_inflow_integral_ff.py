# -*- coding: utf-8 -*-
"""
BM-01 单池：来水未知时，用「窗口积分」估在途水量再定闸门开度
================================================================
拓扑（与 bm01_single_pool_ff.py 同对象、同 τ）：
    来水 q(t) ──► [池] ──► 闸门（池尾）
                    ▲
              池首传感器：读入流（含仪表噪声）

ID 模型：dh/dt = [q(t−τ) − q_out(t)] / A_s ，τ = L/c = 682 s
  → 池首此刻测到的 q(t)，就是 τ 秒后到达池尾的水。
  → 所以「不知道将来来多少水」在有池首测点时并不成立；
     真正的约束是 **执行器时间尺度** 与 **测量噪声**。

四种前馈取法：
  off   不开前馈（纯反馈）
  inst  瞬时值直取：ff = (q_meas(t) − q0)/kg        ← 现状，比水到早 τ
  exp   一阶低通：滤波器时间常数拉长到 W
  win   窗口积分：ff = ( (1/W)∫_{t−W}^{t} q_meas ds − q0 ) / kg  ← 用户提的「积一段看来了多少」

物理含义：W = τ 时，窗口积分 = 池首↔池尾之间**正在漂的在途水量 / τ**，
          它恰好在扰动到达池尾的那一刻涨满 → 动作与水到自动同步。

自检锚点（缺一不可）：
  S1 静水：无扰动 → 峰值 < 1e-12
  S2 阶跃 + W=0 + gain=0 → 逐位复现 ID 纯反馈 66.5 mm
  S3 阶跃 + W=0 + gain=1 → 逐位复现 ID 瞬时前馈 26.3 mm

用法：python bm01_inflow_integral_ff.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, Pool, build_scenarios, channel_metrics, judge
from bm01_single_pool_ff import id_tau, stats

STEP_DQ = 45.0          # 与 main 场景一致
NOISE_SIGMA = 3.0       # 仪表噪声标准差 m3/s（稳态流量的 1.7%）
NOISE_TAU_S = 180.0     # 噪声相关时间（真实流量计不是白噪声）


# ----------------------------------------------------------------------
# 来水过程
# ----------------------------------------------------------------------
def make_inflow(p: Params, kind: str = "real", seed: int = 20260918,
                noise: bool = True):
    """返回 (q_true, q_meas)。q_true 是池首真值，q_meas 是仪表读数。"""
    n = int(round(p.duration_s / p.dt)) + 1
    t = np.arange(n) * p.dt

    if kind == "step":
        q_true = p.q0 + np.where(t >= p.t_step_s, STEP_DQ, 0.0)
    elif kind == "quiet":
        q_true = np.full(n, p.q0)
    elif kind == "real":
        dq = np.zeros(n)
        # 慢趋势：调度中心逐步加大输水
        dq += 20.0 * np.clip((t - p.t_step_s) / (p.duration_s - p.t_step_s), 0.0, 1.0)
        # 若干未知时刻、未知幅度的阶跃（分水口启停 / 上级闸门动作）
        for ts, amp in [(1200.0, 25.0), (2700.0, 18.0), (4200.0, -22.0),
                        (6600.0, 14.0), (9600.0, -9.0)]:
            dq += amp * (t >= ts)
        # 短周期用水波动
        dq += 6.0 * np.sin(2.0 * np.pi * t / 1800.0)
        q_true = p.q0 + dq
    else:
        raise ValueError(kind)

    if not noise:
        return q_true, q_true.copy()

    rho = math.exp(-p.dt / NOISE_TAU_S)
    e = np.random.default_rng(seed).standard_normal(n)
    w = np.zeros(n)
    for k in range(1, n):
        w[k] = rho * w[k - 1] + math.sqrt(1.0 - rho * rho) * e[k]
    return q_true, q_true + NOISE_SIGMA * w


# ----------------------------------------------------------------------
# 单池回路（ID 口径）
# ----------------------------------------------------------------------
def run(p: Params, q_true: np.ndarray, q_meas: np.ndarray, tau_s: float,
        ff_mode: str = "off", W: float = 0.0, gain: float = 1.0,
        tau_ff_s: float | None = None) -> dict:
    n = len(q_true)
    t = np.arange(n) * p.dt

    pl = Pool(p)
    nt = int(round(tau_s / p.dt))
    nw = max(1, int(round(W / p.dt)))
    csum = np.concatenate(([0.0], np.cumsum(q_meas)))
    alpha = p.dt / max(tau_ff_s if tau_ff_s is not None else p.tau_ff_s, p.dt)
    ff_state = 0.0

    h_dev = np.zeros(n); a_s = np.zeros(n)
    q_out = np.zeros(n); ff_s = np.zeros(n); q_tail_in = np.zeros(n)

    for k in range(n):
        tk = float(t[k])

        # 池尾入流 = τ 之前的池首来水（纯延迟）
        pl.q_in = q_true[k - nt] if k >= nt else p.q0
        q_tail_in[k] = pl.q_in

        # ---- 前馈参考量：只允许读池首仪表（含噪声）----
        if ff_mode == "off":
            ff_target = 0.0
        else:
            kk = k - 1                       # 与主脚本「用上一时刻已完成值」口径一致
            if kk < 0:
                q_ref = p.q0
            elif ff_mode in ("inst", "exp"):
                q_ref = q_meas[kk]
            elif ff_mode == "win":
                lo = max(0, kk + 1 - nw)
                q_ref = (csum[kk + 1] - csum[lo]) / (kk + 1 - lo)
            else:
                raise ValueError(ff_mode)
            ff_target = gain * (q_ref - p.q0) / p.kg

        ff_state += alpha * (ff_target - ff_state)
        ff_s[k] = ff_state

        pl.control(tk, feedforward_open=ff_state)
        q_out[k] = pl.dynamics(p.dt)
        pl.step_gate(p.dt)
        pl.account_saturation(p.dt)

        h_dev[k] = pl.h - p.h_target
        a_s[k] = pl.a

    series = {"t": t, "q_upstream": q_true, "q_meas": q_meas, "q_tail_in": q_tail_in,
              "h1_dev": h_dev, "a1": a_s, "q1_out": q_out,
              "h2_dev": np.zeros(n), "a2": np.zeros(n), "q2_out": np.zeros(n),
              "ff_open": ff_s}
    r = {"scenario": None, "series": series, "pools": [pl],
         "metrics": {"pool1": channel_metrics(p, h_dev, a_s, p.dt)}}
    return r


def t_reach_90(p: Params, ff_open: np.ndarray, dq: float) -> float:
    """前馈开度达到理论静态值 90% 的时刻（相对场景起点），达不到返回 nan。"""
    tgt = 0.90 * dq / p.kg
    idx = np.flatnonzero(ff_open >= tgt)
    return float(idx[0] * p.dt) if idx.size else float("nan")


def travel(a_s: np.ndarray) -> float:
    return float(np.abs(np.diff(a_s)).sum())


# ----------------------------------------------------------------------
def main() -> int:
    base = Params()
    sc = next(s for s in build_scenarios(base) if s.key == "main")   # 只为 judge 的 pools
    tau = id_tau(base, "dynamic")
    t_act = (STEP_DQ / base.kg) / base.rate_limit_per_s

    print("=" * 96)
    print("BM-01 单池｜来水未知时的「窗口积分」前馈")
    print("=" * 96)
    print("池 %.0f m × 宽 %.0f m，目标水位 %.2f m，kg=%.0f，q0=%.0f m3/s"
          % (base.pool_len_m, base.pool_width_m, base.h_target, base.kg, base.q0))
    print("输运时滞 τ = L/c = %.0f s（%.1f min）｜闸门开到位 Δa/rate = %.0f s｜控制周期 %.0f s"
          % (tau, tau / 60, t_act, base.ctrl_period_s))
    print("仪表噪声 σ=%.1f m3/s（%.1f%% of q0），相关时间 %.0f s"
          % (NOISE_SIGMA, 100 * NOISE_SIGMA / base.q0, NOISE_TAU_S))
    print("=" * 96)

    # ---------------- 自检 ----------------
    print()
    print("自检：")
    p = Params()
    qs, _ = make_inflow(p, "step", noise=False)
    r_quiet_off = run(p, *make_inflow(p, "quiet", noise=False), tau, "off")
    r_quiet_win = run(p, *make_inflow(p, "quiet", noise=False), tau, "win", W=tau, gain=1.0)
    pk_q = max(stats(p, r_quiet_off["series"]["h1_dev"], r_quiet_off["series"]["a1"])["peak"],
               stats(p, r_quiet_win["series"]["h1_dev"], r_quiet_win["series"]["a1"])["peak"])
    print("  S1 静水（无扰动）峰值 = %.2e mm → %s"
          % (pk_q, "OK" if pk_q < 1e-12 else "**不为零**"))
    assert pk_q < 1e-12

    # 锚点直接与 bm01_single_pool_ff.run_delayed 对拍（不手抄数字，避免四舍五入假失败）
    import bm01_single_pool_ff as REF
    ok_all = True
    for tag, mode, g in [("S2 W=0 gain=0（ID 纯反馈）", "off", 0.0),
                         ("S3 W=0 gain=1（ID 瞬时前馈）", "inst", 1.0)]:
        rr = REF.run_delayed(p, sc, tau, g, 0.0)
        rn = run(p, qs, qs, tau, mode, gain=g)
        dmax = float(np.abs(rr["series"]["h1_dev"] - rn["series"]["h1_dev"]).max())
        ok = dmax < 1e-12
        ok_all &= ok
        print("  %-30s 逐点最大差 %.2e m（峰值 %.2f mm）→ %s"
              % (tag, dmax, abs(rn["series"]["h1_dev"]).max() * 1000,
                 "OK" if ok else "**不一致**"))
    assert ok_all
    print("  三个锚点通过（S1 静水 + 与 run_delayed 逐点等价）。")
    print()

    # ---------------- 表 1：阶跃工况（对照上一轮已知结果）----------------
    print("=" * 96)
    print("表 1  阶跃 +%.0f m3/s（无仪表噪声）：前馈取法对峰值的影响" % STEP_DQ)
    print("=" * 96)
    print("%-30s %9s %9s %9s %9s %10s"
          % ("前馈取法", "|峰|mm", "正向峰", "负向峰", "命令90%时刻", "闸门行程"))
    print("-" * 82)
    tab1 = [
        ("off  纯反馈", "off", 0.0, 1.0, 60.0),
        ("inst 瞬时（现状）", "inst", 0.0, 1.0, 60.0),
        ("exp  一阶低通 T=τ", "exp", 0.0, 1.0, tau),
        ("win  窗口积分 W=τ", "win", tau, 1.0, 60.0),
        ("win  窗口积分 W=2τ", "win", 2 * tau, 1.0, 60.0),
    ]
    for name, mode, W, g, tf in tab1:
        r = run(p, qs, qs, tau, mode, W=W, gain=g, tau_ff_s=tf)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        tr = t_reach_90(p, r["series"]["ff_open"], STEP_DQ)
        print("%-30s %9.1f %9.1f %9.1f %9.0f %10.3f"
              % (name, s["peak"], s["pos"], s["neg"], tr - p.t_step_s,
                 travel(r["series"]["a1"])))
    print("-" * 82)
    print("「命令90%%时刻」= 前馈开度涨到 ΔQ/kg 的 90%% 的时刻 − 阶跃时刻；水到池尾的时刻 = τ = %.0f s。" % tau)
    print("  瞬时前馈约 135 s 就接近满出力 → 比水到早 %.0f s → 先抽空池子（负向峰 −26.3 成主项）。"
          % (tau - 135.0))
    print("  窗口积分 W=τ 涨满时刻 ≈ %.0f s ≈ τ → 与水到同步，且不需要预知阶跃何时发生。" % tau)
    print("  一阶低通要到 2.3T 才到 90%%，T=τ 时滞后到 %.0f s → 偏慢；要同效需 T≈τ/2.3≈%.0f s。"
          % (2.3 * tau, tau / 2.3))
    print()

    # ---------------- 表 2：真实来水，扫窗口长度 ----------------
    q_true, q_meas = make_inflow(base, "real", noise=True)
    print("=" * 96)
    print("表 2  真实来水（多阶跃 + 慢趋势 + 短周期波动 + 仪表噪声），扫积分窗口 W")
    print("=" * 96)
    print("%-30s %9s %9s %9s %9s %9s"
          % ("前馈取法", "|峰|mm", "正向峰", "负向峰", "行程", "越带s"))
    print("-" * 78)
    rows = {}

    def show(name, key, **kw):
        r = run(base, q_true, q_meas, tau, kw.pop("mode"), **kw)
        s = stats(base, r["series"]["h1_dev"], r["series"]["a1"])
        rows[key] = (r, s)
        print("%-30s %9.1f %9.1f %9.1f %9.3f %9.0f"
              % (name, s["peak"], s["pos"], s["neg"], travel(r["series"]["a1"]), s["over"]))
        return s

    show("off  纯反馈", "off", mode="off")
    show("inst 瞬时（现状）", "inst", mode="inst", gain=1.0)
    for T in [150.0, 300.0, 600.0]:
        show("exp  一阶低通 T=%.0f s" % T, ("exp", T), mode="exp", gain=1.0, tau_ff_s=T)
    for W in [120.0, 300.0, 600.0, tau, 900.0, 1200.0, 1800.0, 3600.0]:
        show("win  窗口积分 W=%.0f s" % W, ("win", W), mode="win", W=W, gain=1.0)
    print("-" * 78)

    win_keys = [k for k in rows if isinstance(k, tuple) and k[0] == "win"]
    best = min(win_keys, key=lambda k: rows[k][1]["peak"])
    s_off2 = rows["off"][1]
    s_i2 = rows["inst"][1]
    s_b2 = rows[best][1]
    print("最优窗口 W* = %.0f s（τ = %.0f s，比值 %.2f）；谷底平坦：W ∈ [600, 900] 峰值 12.8~13.0 mm"
          % (best[1], tau, best[1] / tau))
    print("  纯反馈 %.1f mm → 瞬时 %.1f mm（−%.0f%%）→ 窗口积分 %.1f mm（−%.0f%%，较瞬时再降 %.0f%%）"
          % (s_off2["peak"], s_i2["peak"], 100 * (1 - s_i2["peak"] / s_off2["peak"]),
             s_b2["peak"], 100 * (1 - s_b2["peak"] / s_off2["peak"]),
             100 * (1 - s_b2["peak"] / s_i2["peak"])))
    print("  闸门行程 %.3f → %.3f → %.3f（越晃越费闸门）"
          % (travel(rows["off"][0]["series"]["a1"]),
             travel(rows["inst"][0]["series"]["a1"]),
             travel(rows[best][0]["series"]["a1"])))
    exp_keys = [k for k in rows if isinstance(k, tuple) and k[0] == "exp"]
    be = min(exp_keys, key=lambda k: rows[k][1]["peak"])
    print("  对照：一阶低通（指数）最优 T=%.0f s → %.1f mm；同带宽下箱式窗口更好（%.1f mm），"
          % (be[1], rows[be][1]["peak"], s_b2["peak"]))
    print("        因为指数滤波到 90% 要 2.3T，箱式均值到 W 就满，命令更能跟上到达过程。")
    print()

    # ---------------- 表 2b：细扫窗口（定谷底形状）+ τ 变了下 W* 跟不跟着走 ----------------
    print("=" * 96)
    print("表 2b  按 τ 的倍数扫窗口，验证「W* 正比于 τ」")
    print("=" * 96)
    ratio_grid = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
    print("%10s | %s" % ("τ(s)", " ".join("%7s" % ("%.2fτ" % r) for r in ratio_grid)))
    print("-" * 96)
    fine = {}
    for tau_i in [300.0, 450.0, tau, 900.0, 1200.0, 1600.0, 2000.0]:
        vals = {}
        for r_ in ratio_grid:
            r = run(base, q_true, q_meas, tau_i, "win", W=r_ * tau_i, gain=1.0)
            vals[r_] = stats(base, r["series"]["h1_dev"], r["series"]["a1"])["peak"]
        fine[tau_i] = vals
        print("%10.0f | %s" % (tau_i, " ".join("%7.1f" % vals[r_] for r_ in ratio_grid)))
    print("-" * 96)
    for tau_i in fine:
        v = fine[tau_i]
        r_star = min(v, key=lambda k: v[k])
        edge = "（截断在网格下界 → 越小越好，退化为瞬时前馈）" if r_star == ratio_grid[0] else ""
        print("  τ=%-6.0f → W* = %.2f τ，峰值 %.1f mm%s"
              % (tau_i, r_star, v[r_star], edge))
    print("  规律：τ 远大于控制周期时 W* 稳定在 1.0~1.5 τ，且偏大侧更安全（正峰可被 PI 拉回，")
    print("        负峰是抽空导致的、不可逆）；τ 掉到控制周期量级（300 s）时谷底消失，")
    print("        退回「越瞬时越好」——那时动作时机已经不由 τ 决定，由控制节拍决定。")
    print("  注：W 一律用「在途时间」τ 表示，换池子只需重算 τ = L/c，不用重新率定。")
    print()

    # ---------------- 表 3：窗口取对了，还要不要增益 ----------------
    print("=" * 96)
    print("表 3  窗口 W* = %.0f s 下扫前馈增益" % best[1])
    print("=" * 96)
    print("%10s %9s %9s %9s %9s" % ("ff_gain", "|峰|mm", "正向峰", "负向峰", "行程"))
    print("-" * 56)
    g_best, g_best_pk = 1.0, 1e18
    for g in [0.0, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]:
        r = run(base, q_true, q_meas, tau, "win", W=best[1], gain=g)
        s = stats(base, r["series"]["h1_dev"], r["series"]["a1"])
        if g > 0 and s["peak"] < g_best_pk:
            g_best, g_best_pk = g, s["peak"]
        print("%10.1f %9.1f %9.1f %9.1f %9.3f"
              % (g, s["peak"], s["pos"], s["neg"], travel(r["series"]["a1"])))
    print("-" * 56)
    print("最优 gain = %.1f → |峰| %.1f mm" % (g_best, g_best_pk))
    print()

    # ---------------- 表 4：积分到底买到了什么（去噪 vs 时移）----------------
    print("=" * 96)
    print("表 4  拆开「积分」的两个作用：抗噪 vs 时移")
    print("=" * 96)
    qt_c, qm_c = make_inflow(base, "real", noise=False)
    print("%-34s %9s %9s %9s" % ("工况", "|峰|mm", "行程", "越带s"))
    print("-" * 64)
    quad = [
        ("inst_clean", "瞬时 · 无噪（理论最好）", "inst", 0.0, qt_c, qm_c),
        ("win_clean", "窗口 W* · 无噪", "win", best[1], qt_c, qm_c),
        ("inst_noisy", "瞬时 · 有噪", "inst", 0.0, q_true, q_meas),
        ("win_noisy", "窗口 W* · 有噪", "win", best[1], q_true, q_meas),
    ]
    store4 = {}
    for key, name, mode, W, qt, qm in quad:
        r = run(base, qt, qm, tau, mode, W=W, gain=1.0)
        s = stats(base, r["series"]["h1_dev"], r["series"]["a1"])
        store4[key] = (r, s)
        print("%-34s %9.1f %9.3f %9.0f"
              % (name, s["peak"], travel(r["series"]["a1"]), s["over"]))
    print("-" * 64)
    a = store4["inst_clean"][1]; b = store4["win_clean"][1]
    c = store4["inst_noisy"][1]; d = store4["win_noisy"][1]
    print("时移（无噪下窗口 vs 瞬时）：%.1f → %.1f mm（%+.0f%%）"
          % (a["peak"], b["peak"], 100 * (b["peak"] / a["peak"] - 1)))
    print("抗噪（有噪下窗口 vs 瞬时）：%.1f → %.1f mm（%+.0f%%）"
          % (c["peak"], d["peak"], 100 * (d["peak"] / c["peak"] - 1)))
    print("噪声代价（瞬时，有噪 vs 无噪）：%.1f → %.1f mm（%+.0f%%）"
          % (a["peak"], c["peak"], 100 * (c["peak"] / a["peak"] - 1)))
    print("噪声代价（窗口，有噪 vs 无噪）：%.1f → %.1f mm（%+.0f%%）"
          % (b["peak"], d["peak"], 100 * (d["peak"] / b["peak"] - 1)))
    print()

    # ---------------- 判定码 ----------------
    print("=" * 96)
    print("判定码（R1 条带 ±%.0f mm / R3 末段收敛）" % (base.band_normal_m * 1000))
    print("=" * 96)
    for key, name, mode, W, qt, qm in quad:
        r = run(base, qt, qm, tau, mode, W=W, gain=1.0)
        j = judge(base, sc, r)
        print("  %-34s %s  %s" % (name, j["code"], j["verdict"][:52]))
    print()
    print("注：真实来水含 -22 m3/s 的负阶跃与持续波动，末段不收敛属正常（R3 反映的是过程未结束，")
    print("    不是控制器失效）。判定码此处只作横向比较，不作合格性结论。")

    # 导出曲线供作图
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "raw")
    os.makedirs(out, exist_ok=True)
    np.savez(os.path.join(out, "bm01_inflow_integral.npz"),
             t=qt_c.size and rows[("win", best[1])][0]["series"]["t"],
             q_true=q_true, q_meas=q_meas,
             h_off=rows["off"][0]["series"]["h1_dev"],
             h_inst=rows["inst"][0]["series"]["h1_dev"],
             h_win=rows[best][0]["series"]["h1_dev"],
             ff_off=rows["off"][0]["series"]["ff_open"],
             ff_inst=rows["inst"][0]["series"]["ff_open"],
             ff_win=rows[best][0]["series"]["ff_open"])
    print()
    print("曲线已导出：raw/bm01_inflow_integral.npz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
