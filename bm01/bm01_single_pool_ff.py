# -*- coding: utf-8 -*-
"""
BM-01 单池单闸门：池首传感器测入流 → 池尾闸门提前动作
================================================================
拓扑（与 bm01_pool_step 的 main 场景同对象）：
    来水 q_up ──► [池] ──► 闸门（池尾）
                    ▲
                池首传感器：直接读 q_up，零信息延迟

对照：
  A 纯反馈        闸门只看本池水位（= bm01_pool_step 的 main，ff_gain=0）
  B 前馈 + 反馈   闸门开度额外加 ff_gain * (q_up - q0) / kg

自检锚点：ff_gain=0 / delay=0 必须逐位复现 main 场景的池1 指标。

用法：
  python bm01_single_pool_ff.py              # 全表
  python bm01_single_pool_ff.py --no-limit   # 顺带看「执行器不限速」时前馈能到多少
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import (Params, Pool, build_scenarios, channel_metrics, judge,
                            gate_flow)

G = 9.81


def id_tau(p: Params, mode: str = "dynamic") -> float:
    """ID 模型的输运时滞 τ = L/c。

    dynamic : c = v + √(g·A/B)（动力波速，只需要几何，水位扰动的传播速度）
    gate    : c = (1/B)·dQ/dh，dQ/dh 取闸门率定曲线
              —— BM-01 没有渠道率定曲线（无坡度/糙率/边坡），只能用闸门这条，
                 结果是「蓄量时间常数 A_s/(dQ/dh)」而不是波速，量纲对但物理含义不同。
    """
    h, B = p.h_target, p.pool_width_m
    A = B * h
    v = p.q0 / A
    if mode == "dynamic":
        c = v + math.sqrt(G * A / B)
    elif mode == "gate":
        dqdh = (gate_flow(p, p.a0, h + 1e-4) - gate_flow(p, p.a0, h - 1e-4)) / 2e-4
        c = dqdh / B
    else:
        raise ValueError(mode)
    return p.pool_len_m / c


def run(p: Params, sc, ff_gain: float = 0.0, ff_delay_s: float = 0.0) -> dict:
    """单池回路。ff_gain / ff_delay_s 只影响池尾闸门的前馈开度。

    前馈参考量取 q_up（池首来水），nd=0 时表示"传感器就在池首、信息瞬时到位"。
    kd = k-1-nd：与 bm01_pool_step 里 q1_out[k-1] 的口径一致（用上一时刻的已完成值）。
    """
    n = int(round(p.duration_s / p.dt)) + 1
    t = np.arange(n) * p.dt

    pool = Pool(p)
    ff_state = 0.0
    nd = int(round(ff_delay_s / p.dt))
    alpha = p.dt / max(p.tau_ff_s, p.dt)

    q_up_s = np.zeros(n)
    h_dev = np.zeros(n)
    a_s = np.zeros(n)
    q_out = np.zeros(n)
    ff_s = np.zeros(n)
    q_up = p.q0

    for k in range(n):
        tk = float(t[k])
        q_up = p.q0 + (sc.step_dq if tk >= p.t_step_s else 0.0)
        pool.q_in = q_up
        q_up_s[k] = q_up

        kd = k - 1 - nd
        q_ref = q_up_s[kd] if kd >= 0 else p.q0
        ff_target = ff_gain * (q_ref - p.q0) / p.kg
        ff_state += alpha * (ff_target - ff_state)
        ff_s[k] = ff_state

        pool.control(tk, feedforward_open=ff_state)
        q_out[k] = pool.dynamics(p.dt)
        pool.step_gate(p.dt)
        pool.account_saturation(p.dt)

        h_dev[k] = pool.h - p.h_target
        a_s[k] = pool.a

    series = {"t": t, "q_upstream": q_up_s, "h1_dev": h_dev, "a1": a_s,
              "q1_out": q_out, "h2_dev": np.zeros(n), "a2": np.zeros(n),
              "q2_out": np.zeros(n), "ff_open": ff_s}
    r = {"scenario": sc, "series": series, "pools": [pool],
         "metrics": {"pool1": channel_metrics(p, h_dev, a_s, p.dt)}}
    r["judgement"] = judge(p, sc, r)
    return r


def stats(p: Params, h_dev: np.ndarray, a_s: np.ndarray) -> dict:
    dev = np.abs(h_dev)
    pos, neg = float(h_dev.max()) * 1000, float(h_dev.min()) * 1000
    travel = float(np.abs(np.diff(a_s)).sum())
    k = int(np.argmax(dev))
    return {
        "pos": pos, "neg": neg, "peak": float(dev.max()) * 1000,
        "t_peak": float(k * p.dt), "travel": travel,
        "over": float(np.sum(dev > p.band_normal_m) * p.dt),
        "final": float(dev[int(round((p.duration_s - p.settle_window_s) / p.dt)):].max()) * 1000,
    }


# ----------------------------------------------------------------------
# 带输运时滞的单池（= ID 模型）：池首测到来水后，池尾要等 τ 才感觉到
#     dh/dt = [ q_in(t-τ) - q_out(t) ] / A
# 集总模型里 τ=0，全池只有一个水位，"池首放传感器"和"池尾放传感器"没区别。
# 只有 τ>0 时，「池首传感器 + 池尾闸门」才是一个真实的提前量。
# ----------------------------------------------------------------------
def run_delayed(p: Params, sc, tau_s: float,
                ff_gain: float = 0.0, ff_delay_s: float = 0.0) -> dict:
    n = int(round(p.duration_s / p.dt)) + 1
    t = np.arange(n) * p.dt

    pl = Pool(p)                       # 复用控制律 / 限速 / 饱和统计 / 蓄量积分
    nt = int(round(tau_s / p.dt))      # 输运时滞步数
    nd = int(round(ff_delay_s / p.dt))
    alpha = p.dt / max(p.tau_ff_s, p.dt)
    ff_state, q_up = 0.0, p.q0

    q_up_s = np.zeros(n); h_dev = np.zeros(n); a_s = np.zeros(n)
    q_out = np.zeros(n); ff_s = np.zeros(n); q_in_s = np.zeros(n)

    for k in range(n):
        tk = float(t[k])
        q_up = p.q0 + (sc.step_dq if tk >= p.t_step_s else 0.0)
        q_up_s[k] = q_up

        # 前馈：池首传感器零延迟读到 q_up（这才是"传感器在池上游"）
        kd = k - 1 - nd
        q_ref = q_up_s[kd] if kd >= 0 else p.q0
        ff_state += alpha * (ff_gain * (q_ref - p.q0) / p.kg - ff_state)
        ff_s[k] = ff_state

        # 池尾入流 = τ 之前的来水
        pl.q_in = q_up_s[k - nt] if k >= nt else p.q0
        q_in_s[k] = pl.q_in

        pl.control(tk, feedforward_open=ff_state)
        q_out[k] = pl.dynamics(p.dt)
        pl.step_gate(p.dt)
        pl.account_saturation(p.dt)

        h_dev[k] = pl.h - p.h_target
        a_s[k] = pl.a

    series = {"t": t, "q_upstream": q_up_s, "q_tail_in": q_in_s, "h1_dev": h_dev,
              "a1": a_s, "q1_out": q_out, "h2_dev": np.zeros(n), "a2": np.zeros(n),
              "q2_out": np.zeros(n), "ff_open": ff_s}
    r = {"scenario": sc, "series": series, "pools": [pl],
         "metrics": {"pool1": channel_metrics(p, h_dev, a_s, p.dt)}}
    r["judgement"] = judge(p, sc, r)
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-limit", action="store_true",
                    help="另跑一组「闸门不限速」，用于判断瓶颈是否在执行器")
    args = ap.parse_args()

    base = Params()
    sc = next(s for s in build_scenarios(base) if s.key == "main")
    print("=" * 90)
    print("BM-01 单池单闸门｜池首传感器 + 池尾闸门｜阶跃 +%.0f m3/s" % sc.step_dq)
    print("池长 %.0f m × 水面宽 %.0f m = %.2e m²｜目标水位 %.2f m｜初始开度 %.2f"
          % (base.pool_len_m, base.pool_width_m, base.area_m2, base.h_target, base.a0))
    print("闸门系数 kg=%.0f｜向上余量 %.0f m3/s｜开度限速 %.4f /s（满行程 %.0f min）"
          % (base.kg, base.gate_up_capacity, base.rate_limit_per_s, 1.0 / base.rate_limit_per_s / 60))
    print("理想前馈开度 = ΔQ/kg = %.0f/%.0f = %.4f（静态等价）"
          % (sc.step_dq, base.kg, sc.step_dq / base.kg))
    print("=" * 90)

    # ---------- 自检：ff_gain=0 必须等于 main ----------
    p0 = Params()
    r_ref = run(p0, sc, 0.0, 0.0)
    from bm01_pool_step import simulate as sim_official
    r_off = sim_official(p0, sc)
    same = (abs(r_ref["metrics"]["pool1"]["peak_dev_m"]
                - r_off["metrics"]["pool1"]["peak_dev_m"]) < 1e-12
            and abs(r_ref["metrics"]["pool1"]["final_window_max_dev_m"]
                    - r_off["metrics"]["pool1"]["final_window_max_dev_m"]) < 1e-12)
    print("自检 ff_gain=0 vs main 场景：peak %.5f m / %.5f m，末段 %.5f / %.5f → %s"
          % (r_ref["metrics"]["pool1"]["peak_dev_m"], r_off["metrics"]["pool1"]["peak_dev_m"],
             r_ref["metrics"]["pool1"]["final_window_max_dev_m"],
             r_off["metrics"]["pool1"]["final_window_max_dev_m"],
             "OK" if same else "**不一致**"))
    assert same, "单池回路与原 simulate() 不等价，先修脚本"
    print()

    # ---------- 增益扫描 ----------
    print("池首传感器（零信息延迟），扫前馈增益：")
    print("%8s %10s %10s %10s %10s %10s %10s %8s"
          % ("ff_gain", "正向峰", "负向峰", "|峰|mm", "越带s", "末段mm", "行程", "判定"))
    print("-" * 90)
    rows = {}
    for g in [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]:
        p = Params()
        p.ff_gain = g
        r = run(p, sc, g, 0.0)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        rows[g] = (r, s)
        mark = ""
        print("%8.2f %10.1f %10.1f %10.1f %10.0f %10.1f %10.3f %8s%s"
              % (g, s["pos"], s["neg"], s["peak"], s["over"], s["final"],
                 s["travel"], r["judgement"]["code"], mark))
    print("-" * 90)
    g_best = min(rows, key=lambda k: rows[k][1]["peak"])
    print("峰值最优 ff_gain = %.2f → |峰| %.1f mm（纯反馈 %.1f mm，降 %.0f%%），越带 %.0f s，判定 %s"
          % (g_best, rows[g_best][1]["peak"], rows[0.0][1]["peak"],
             100 * (1 - rows[g_best][1]["peak"] / rows[0.0][1]["peak"]),
             rows[g_best][1]["over"], rows[g_best][0]["judgement"]["code"]))

    # ---------- 信息延迟敏感性（传感器挪远多少才失效）----------
    print()
    print("传感器位置 = 信息延迟；固定 ff_gain=%.2f：" % g_best)
    print("%10s %10s %10s %10s" % ("延迟s", "|峰|mm", "越带s", "行程"))
    print("-" * 46)
    for d in [0, 60, 120, 180, 300, 420, 600, 900]:
        p = Params()
        p.ff_gain = g_best
        r = run(p, sc, g_best, float(d))
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        print("%10d %10.1f %10.0f %10.3f" % (d, s["peak"], s["over"], s["travel"]))
    print("-" * 46)

    # ---------- 瓶颈判定：执行器限速 ----------
    if args.no_limit:
        print()
        print("把闸门限速放开（rate_limit_per_s=1.0，等于瞬时到位）：")
        print("%8s %10s %10s %10s" % ("ff_gain", "|峰|mm", "越带s", "行程"))
        print("-" * 46)
        for g in [0.0, 0.6, 0.8, 1.0, 1.2]:
            p = replace(Params(), ff_gain=g, rate_limit_per_s=1.0)
            r = run(p, sc, g, 0.0)
            s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
            print("%8.2f %10.1f %10.0f %10.3f" % (g, s["peak"], s["over"], s["travel"]))
        print("-" * 46)

    # ---------- 带输运时滞的单池：池首/池尾真正分离 ----------
    print()
    print("=" * 92)
    print("单池 + 输运时滞 τ：池首传感器在 τ 秒前就知道要来水（τ=0 即集总，无提前量）")
    print("=" * 92)

    p = Params()
    r0 = run_delayed(p, sc, 0.0, 0.0, 0.0)
    s0 = stats(p, r0["series"]["h1_dev"], r0["series"]["a1"])
    lump = rows[0.0][1]["peak"]
    ok = abs(s0["peak"] - lump) < 1e-9
    print("自检 τ=0 必须退化回集总单池：|峰| %.2f mm vs %.2f mm → %s"
          % (s0["peak"], lump, "OK" if ok else "**不一致**"))
    assert ok, "τ=0 未退化回集总模型"
    print()

    print("%8s | %10s %10s %10s | %10s %9s %8s"
          % ("τ(s)", "纯反馈峰", "前馈d=0", "前馈时机匹配", "最优delay", "行程", "判定"))
    print("-" * 92)
    for tau in [0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2400.0, 3000.0, 3600.0]:
        p = Params()
        r_off = run_delayed(p, sc, tau, 0.0, 0.0)
        s_off = stats(p, r_off["series"]["h1_dev"], r_off["series"]["a1"])
        best = (None, 1e18, None, None)
        for d in np.arange(0.0, tau + 900.0 + 1.0, 300.0):
            r = run_delayed(p, sc, tau, 1.0, float(d))
            s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
            if s["peak"] < best[1]:
                best = (float(d), s["peak"], s["travel"], r)
        r_zero = run_delayed(p, sc, tau, 1.0, 0.0)
        s_zero = stats(p, r_zero["series"]["h1_dev"], r_zero["series"]["a1"])
        print("%8.0f | %10.1f %10.1f %10.1f | %10.0f %9.3f %8s"
              % (tau, s_off["peak"], s_zero["peak"], best[1], best[0],
                 best[2], best[3]["judgement"]["code"]))
    print("-" * 92)
    print("注意「纯反馈峰」那一列**不是**常数：τ 取 300 的整数倍时全踩在控制节拍上（有效滞后 = 300 s）")
    print("  → 恒 72.7 mm；τ 落在节拍中间时反而更低（有效滞后 = 300 − 相位，最小 65.1 mm）。")
    print("  上一轮写「纯反馈峰值与 τ 无关」是采样错觉，已在此更正。")
    print("  「前馈时机匹配」地板与 τ 无关 = 只要把动作按 τ 平移就够。")

    # ---------- ID 模型：τ 由波速算出，不是扫的 ----------
    print()
    print("=" * 92)
    print("ID 模型（积分时滞）：dh/dt = [q_in(t−τ) − q_out(t)] / A_s，τ = L/c")
    print("=" * 92)
    h, B = base.h_target, base.pool_width_m
    area = B * h
    v = base.q0 / area
    c_d = v + math.sqrt(G * area / B)
    tau_dyn = id_tau(base, "dynamic")
    tau_gate = id_tau(base, "gate")
    print("几何：池长 L=%.0f m，水面宽 B=%.0f m，水深 h=%.2f m → A=%.0f m²，v=Q/A=%.3f m/s"
          % (base.pool_len_m, B, h, area, v))
    print("  动力波速 c_d = v + √(gA/B) = %.3f + %.3f = %.3f m/s  → τ = %.0f s = %.1f min"
          % (v, math.sqrt(G * area / B), c_d, tau_dyn, tau_dyn / 60))
    dqdh = (gate_flow(base, base.a0, h + 1e-4) - gate_flow(base, base.a0, h - 1e-4)) / 2e-4
    print("  闸门率定口径 dQ/dh = %.2f m²/s → c = %.4f m/s → τ = %.0f s = %.1f h"
          % (dqdh, dqdh / B, tau_gate, tau_gate / 3600))
    print("  ↑ BM-01 无渠道率定曲线（无坡度/糙率/边坡），这条算出来是「蓄量时间常数」不是波速；"
          "且 %.1f h > 仿真时长 %.1f h，池尾全程收不到扰动 → 不可用。**取动力波速口径。**"
          % (tau_gate / 3600, base.duration_s / 3600))
    print()

    p = Params()
    r_chk = run_delayed(p, sc, tau_dyn, 0.0, 0.0)
    k_chk = int((base.t_step_s + tau_dyn - 60.0) / base.dt)
    d_chk = abs(r_chk["series"]["h1_dev"][k_chk])
    ok = d_chk < 1e-12
    print("自检 纯时滞（不是一阶滞后）：t = t_step+τ−60 s 处池尾偏差 = %.2e m → %s"
          % (d_chk, "OK" if ok else "**非零，退化成了一阶滞后**"))
    assert ok, "时滞不是纯延迟"
    print()

    print("对照组（同一对象、同一阶跃 +%.0f m³/s、同一判定码）：" % sc.step_dq)
    print("%-26s %10s %10s %10s %8s %9s"
          % ("工况", "|峰|mm", "越带s", "末段mm", "判定", "闸门行程"))
    print("-" * 80)
    cases = [
        ("A 集总·纯反馈", 0.0, 0.0, 0.0),
        ("B 集总·前馈", 0.0, 1.0, 0.0),
        ("C ID·纯反馈", tau_dyn, 0.0, 0.0),
        ("D ID·前馈（信号一到就动）", tau_dyn, 1.0, 0.0),
    ]
    store = {}
    for name, tau, g, d in cases:
        p = Params()
        r = run_delayed(p, sc, tau, g, d)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        store[name] = (r, s)
        print("%-26s %10.1f %10.0f %10.1f %8s %9.3f"
              % (name, s["peak"], s["over"], s["final"],
                 r["judgement"]["code"], s["travel"]))
    print("-" * 80)

    # 时机（ID 口径下重新寻优）
    print()
    print("ID 模型下扫前馈时机（ff_gain=1.0，τ=%.0f s）：" % tau_dyn)
    print("%12s %11s %11s %11s %11s" % ("前馈延迟s", "正向峰", "负向峰", "|峰|mm", "行程"))
    print("-" * 62)
    best_d, best_pk = 0.0, 1e18
    for d in [0.0, 150.0, 300.0, 450.0, 600.0, 750.0, 900.0, 1200.0]:
        p = Params()
        r = run_delayed(p, sc, tau_dyn, 1.0, d)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        if s["peak"] < best_pk:
            best_d, best_pk = d, s["peak"]
        print("%12.0f %11.1f %11.1f %11.1f %11.3f" % (d, s["pos"], s["neg"], s["peak"], s["travel"]))
    print("-" * 62)
    print("最优时机 %.0f s（|峰| %.1f mm）；对照 τ − 闸门动作滞后 = %.0f − 600 = %.0f s"
          % (best_d, best_pk, tau_dyn, max(0.0, tau_dyn - 600.0)))

    # 增益（在最优时机上）
    print()
    print("ID 模型下扫前馈增益（delay=%.0f s）：" % best_d)
    print("%10s %11s %11s %11s %11s %8s" % ("ff_gain", "正向峰", "负向峰", "|峰|mm", "行程", "判定"))
    print("-" * 70)
    for g in [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5]:
        p = Params()
        p.ff_gain = g
        r = run_delayed(p, sc, tau_dyn, g, best_d)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        print("%10.2f %11.1f %11.1f %11.1f %11.3f %8s"
              % (g, s["pos"], s["neg"], s["peak"], s["travel"], r["judgement"]["code"]))
    print("-" * 70)

    # 未到最优时机，前馈提前量还能再挖多少 —— 拆控制节拍与执行器
    print()
    print("ID 模型（τ=%.0f s）下拆瓶颈：" % tau_dyn)
    print("%-34s %10s %10s %8s %9s"
          % ("工况", "|峰|mm", "越带s", "判定", "闸门行程"))
    print("-" * 78)
    dec = [
        ("ID 纯反馈（现口径）", tau_dyn, 0.0, 0.0, None, None),
        ("ID 纯反馈 · 控制周期 5 s", tau_dyn, 0.0, 0.0, 5.0, None),
        ("ID 纯反馈 · 放开限速", tau_dyn, 0.0, 0.0, None, 1.0),
        ("ID 前馈 d=%.0f · 限速 50 min" % best_d, tau_dyn, 1.0, best_d, None, None),
        ("ID 前馈 d=%.0f · 放开限速" % best_d, tau_dyn, 1.0, best_d, None, 1.0),
    ]
    for name, tau, g, d, per, rate in dec:
        p = Params()
        if per is not None:
            p = replace(p, ctrl_period_s=per)
        if rate is not None:
            p = replace(p, rate_limit_per_s=rate)
        p.ff_gain = g
        r = run_delayed(p, sc, tau, g, d)
        s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
        print("%-34s %10.1f %10.0f %8s %9.3f"
              % (name, s["peak"], s["over"], r["judgement"]["code"], s["travel"]))
    print("-" * 78)
    print()
    print("注：τ 由波速公式给出（动力口径），不再是自由参数；上面的 τ 扫描表仅作敏感性。")
    print("    前馈命令在 Pool.control() 内下发，与水位测点同 300 s 节拍 → 有效提前量只能取"
          " {%.0f, %.0f, %.0f} s 三档（对应前馈延迟 0/300/600 s），最优提前量 %.0f s。"
          % (tau_dyn, tau_dyn - 300.0, max(0.0, tau_dyn - 600.0), tau_dyn - best_d))
    print("    该最优提前量 ≈ 闸门开到位所需时间 Δa/rate = %.0f s：过早 → 抽空池子（负向峰），"
          "过晚 → 水到了闸门还没开完。" % ((sc.step_dq / base.kg) / base.rate_limit_per_s))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
