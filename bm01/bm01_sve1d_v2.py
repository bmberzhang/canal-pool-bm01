# -*- coding: utf-8 -*-
"""
BM-01 一维模型 v2：把"更贴近实际"变成可量化的开关
================================================================================
v1 是什么
---------
v1（`bm01_sve1d_compare.py`）已经把池内 ∂/∂x 项拿回来了（有限体积 + Rusanov +
显式欧拉），结论是：集总 0D 把峰值低估约 2.2×，差异全部来自"闸门量到的是局部水位"
而非池均水位（160.8 mm vs 72.7 mm）。v1 用 120/240/480 三级网格，逐级增量
13.7 → 2.1 mm（正在收敛），故其 Richardson 外推 160.8 mm 是可信的。

v2 要回答的问题
---------------
v1 还剩哪些"不贴近实际"的地方？每一处值多少 mm？本文件把它们做成开关，逐项量化。

开关清单（默认全关 = 与 v1 逐行等价，用作自检基准）
-----------------------------------------------------
S2  断面改梯形（边坡 m）—— **有效**
    v1 用矩形 A=B·h、水面宽 B 恒为池宽。真实渠池多为梯形：
      A(h) = h(b + mh)，水面宽 B(h) = b + 2mh，静水压力项 g·I1(h)（I1=对水面静矩）。
    边坡把"水位涨、水面宽跟着涨"带进来 → 同样来水，蓄量增加更快、水位涨得更少
    → 峰值天然被压小。本开关量出的就是这个效应。
    注：S2 同时启用完整水力半径 R=A/P（梯形下才有区别）。

关掉的两项（已实测确认无物理效果，不列入成果）
-------------------------------------------------
· 完整水力半径 R=A/P（矩形断面）：m=0 时 P=b+2h 与 v1 的 (B+2h) **恒等**，
  增量必然为 0。不是"影响小"，是同一个式子，不该拿来充数。
· 闸门改内部建筑物：本包闸室单元占全池水面 1/120，把出口从边界搬进单元后
  影响在本参数下不可辨（<0.5 mm）。要让它有实质区别，需要闸门**局部收缩**
  （侧收缩、胸墙、淹没出流）这类真实结构——留待有断面资料时再做。

输出
----
raw/bm01_sve1d_v2.json —— 逐开关峰值、末段、越带、相对基准的增量与收敛值。
raw/bm01_sve1d_v2_sensor.json —— 传感器通道网格无关性诊断（单点 vs 闸前固定长度平均）。

**2026-10-01 修正（三处，均为口径错误，结论方向不变）**
1. **纵坡**：梯形档此前沿用 v1 的**矩形反解**纵坡（`V1.normal_slope`，P=b+2h 低估了梯形湿周
   b+2h√(1+m²)），反解出的 S0 偏大 1.38/1.52/1.79 倍（m=1/1.5/2.5）⇒ 初值 (h_target, q0)
   不是均匀流、净驱动力 S0−Sf≠0。现改用 `section_normal_slope()`，旧值在 `legacy_rect_slope` 留档。
2. **峰值指标**：旧表用 120/240 **两级**网格做 Richardson 外推（2×156.61−142.88=170.3 mm），
   等于丢掉了 v1 已有的 480 格那一级——**用最粗两级外推，必然高估**。现改为网格序列 + 最细网格
   收敛值，不再做两点外推。
3. **传感器通道**：单点 h[-1] 的回路增益 ∝ 1/dx，网格加密到 960 格时越过稳定边界，
   出现**晚期增幅极限环**（全程峰伪高到 212 mm、闸前 527 mm；同一跑的前 2h 峰仍为 158.7 mm）。
   现改用闸前固定物理长度平均（`gate_avg`，默认 100 m）——响应速率分母为 B·L_s，与 dx 无关。

**因此峰值指标取 `early_peak_mm`（前 2h 窗）**：物理峰出现在 t≈36 min，而上述晚期病态出现在
t>150 min，取前 2h 窗即可与晚期伪峰解耦。实测该指标对网格与传感器类型都不敏感（漂移 ≤1.3%）。

口径（不得外推）
----------------
证据等级 `synthetic-1d-sve`，**不改变 BM-01 的判定码与证据等级**，
参数仍为量级示例（assume）、非中线现场数据、不构成工程结论。
本文件只为把"每个物理简化值多少"变成一个具体数字（论文"方法验证/敏感性"用）。

用法
----
  python bm01_sve1d_v2.py            # 跑全部开关
  python bm01_sve1d_v2.py --selftest # 只跑自检（开关全关 vs v1）
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import bm01_sve1d_compare as V1

G = V1.G
MANNING_N = V1.MANNING_N
CFL = V1.CFL


# ======================================================================
# 断面几何：矩形（m=0）/ 梯形（m>0）
# ======================================================================
class Section:
    """A(h)、B(h)、P(h)、静水压力积分 I1(h)=∫0^h (h−η)B(η)dη。"""

    def __init__(self, b: float, m: float = 0.0):
        self.b = b        # 底宽
        self.m = m        # 边坡（水平:垂直），m=0 → 矩形

    def A(self, h):
        return h * (self.b + self.m * h)

    def B(self, h):
        return self.b + 2.0 * self.m * h

    def P(self, h):
        # 湿周：底宽 + 两侧坡面
        return self.b + 2.0 * h * math.sqrt(1.0 + self.m * self.m)

    def I1(self, h):
        # ∫0^h (h−η)·B(η)dη = b h²/2 + m h³/3
        return 0.5 * self.b * h * h + (self.m / 3.0) * h * h * h

    def is_rect(self):
        return abs(self.m) < 1e-12


def section_normal_slope(sec: "Section", h: float, Q: float, n: float) -> float:
    """按**实际断面**反解纵坡，使 (h, Q) 恰为均匀流（正常水深）。

    矩形下与 V1.normal_slope 恒等（P=b+2h、A=b·h），故不改变任何既有结论；
    梯形下两者差 2 倍以上，必须用本式，否则初值根本不是均匀流。
    """
    A = sec.A(h)
    P = sec.P(h)
    R = A / P if P > 0 else h
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


# ======================================================================
# 核心求解器 v2
# ======================================================================
def sve_1d_v2(p: dict, ncell: int, step_dq: float, *,
              ctrl_on: str = "gate",
              m_slope: float = 0.0,
              n: float = MANNING_N,
              cfl: float = CFL,
              s0: float = None,
              sensor_len_m: float = 100.0,
              early_window_s: float = 7200.0,
              estimator: str = None) -> dict:
    """
    一维圣维南（守恒形式 + Rusanov 通量 + 显式欧拉）。
    m_slope=0 时与 v1 逐行等价（自检保证）。

    s0：渠底纵坡。None（默认）= 沿用 V1.normal_slope 的矩形反解，行为与既有版本完全一致；
        传入数值则使用给定坡度——梯形断面下必须传 section_normal_slope 的结果，
        否则 (h_target, q0) 不是均匀流，初值会带一个假瞬变。

    ctrl_on：控制通道。
        "gate"     = 单点传感器（最末单元 h[-1]，历史行为）——**回路增益 ∝ 1/dx**，
                     网格越细越容易在 300 s 采样下失稳（实测 960 格出现增幅极限环）。
        "gate_avg" = 闸前固定物理长度 sensor_len_m 内取平均（**回路增益与 dx 无关**）。
        "mean"     = 全池平均（对照组）。
        "off"      = 不开控制器。

    新增（2026-10-12，E-OBS 可观测性实验）：ctrl_on="est" —— 用**工程上真能装**的估计量
        替代理想全池平均，检验判定码 A1 是否保得住：
        estimator="midpoint" = 单点池中点水位 h[ncell//2]（1 支水位计）
        estimator="mp3"      = 池内 3 点等距平均（0.25L/0.5L/0.75L，3 支水位计）
        estimator="mp5"      = 池内 5 点等距平均（5 支水位计）
        estimator="volbal"   = 容积平衡估计 ĥ = h0 + ∫(q_in−q_gate)dt /(B0·L)（**不需水位计**，
                               只需上下游流量计；与 0D 积分蓄量模型同构）。
    estimator 仅在 ctrl_on="est" 时生效；其它通道行为与本文件历史版本逐位一致。
    sensor_len_m：仅 ctrl_on="gate_avg" 生效，默认 100 m。

    2026-10-01 修正：
      (1) 梯形档必须传 s0=section_normal_slope(...)，否则初值不是均匀流；
      (2) 单点传感器在细网格上会失稳，"gate_avg" 是网格无关的替代；
      (3) 峰值指标增加 early_peak（前 2 h 内），把"物理瞬变峰"与"晚期数值伪峰"分开。
    """
    B0, L, ht = p["pool_width_m"], p["pool_len_m"], p["h_target"]
    sec = Section(B0, m_slope)
    dx = L / ncell
    S0 = V1.normal_slope(p, n) if s0 is None else s0
    h = [ht] * ncell
    Q = [p["q0"]] * ncell
    # 固定物理长度的闸前平均传感器（网格无关）
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))

    wmax = p["q0"] / (B0 * ht) + math.sqrt(G * ht)
    lo = cfl * dx / wmax
    dt = V1.pick_dt(dx, wmax) if abs(cfl - CFL) < 1e-15 else None
    if dt is None:
        dt = next(c for c in V1.DT_CANDIDATES if c <= lo)
    nt = int(round(p["duration_s"] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = V1.Ctrl(p)

    out = {"t": [], "dev_mean": [], "dev_gate": [], "dev_head": [], "gate": [], "q_out": []}
    peak_mean = peak_gate = 0.0
    early_peak = 0.0
    tail_mean = 0.0
    est_track_err = 0.0
    t_tail = p["duration_s"] - p["settle_window_s"]
    over_band_dt = 0
    vol_net = 0.0
    c1 = dt / dx
    sqrt = math.sqrt
    hmin = 1e-3
    # 完整水力半径开关：仅当梯形时才有别于 v1 的宽浅近似
    friction_correct = not sec.is_rect()

    for k in range(nt):
        tk = k * dt
        q_up = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)
        q_gate = p["kg"] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)

        if ctrl_on != "off":
            if ctrl_on == "gate":
                meas = h[-1]                       # 单点（历史行为，回路增益 ∝ 1/dx）
            elif ctrl_on == "gate_avg":
                meas = sum(h[ncell - nsens:]) / nsens   # 固定物理长度，dx 无关
            elif ctrl_on == "est":
                if estimator == "midpoint":
                    meas = h[ncell // 2]                       # 1 支水位计（池中点）
                elif estimator == "mp3":
                    idx = [ncell // 4, ncell // 2, (3 * ncell) // 4]
                    meas = (h[idx[0]] + h[idx[1]] + h[idx[2]]) / 3.0
                elif estimator == "mp5":
                    idx = [(j * ncell) // 6 for j in (1, 2, 3, 4, 5)]
                    meas = sum(h[i] for i in idx) / 5.0
                elif estimator == "volbal":
                    meas = ht + vol_net / (B0 * L)              # 容积平衡（只用流量计）
                else:
                    raise ValueError("未知 estimator: %r" % (estimator,))
                est_track_err = max(est_track_err, abs(meas - sum(h) / ncell))
            else:
                meas = sum(h) / ncell              # 全池平均（对照）
            ctrl.control(tk, meas)
            ctrl.step_gate(dt)

        # ---------- 数值通量 ----------
        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif
        F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = hx[j], hx[j + 1]
            if h1 < hmin: h1 = hmin
            if h2 < hmin: h2 = hmin
            A1, A2 = sec.A(h1), sec.A(h2)
            if A1 < 1e-6: A1 = 1e-6
            if A2 < 1e-6: A2 = 1e-6
            Bc1, Bc2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / Bc1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / Bc2)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            # 静水压力项：矩形 g·A²/(2B)（=g·I1）；梯形必须用 I1
            P1 = 0.5 * G * A1 * A1 / Bc1
            P2 = 0.5 * G * A2 * A2 / Bc2
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2 + P1 + P2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate

        # ---------- 单元推进 ----------
        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            A = sec.A(hi)
            if A < 1e-6: A = 1e-6
            Qi = Q[i]
            if friction_correct:
                R = A / sec.P(hi) if sec.P(hi) > 0 else hi
            else:
                R = A / (B0 + 2.0 * hi)                     # v1 的宽浅近似
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        vol_net += (q_up - q_gate) * dt

        dm = sum(h) / ncell - ht
        dg = h[-1] - ht
        if abs(dm) > peak_mean: peak_mean = abs(dm)
        if abs(dg) > peak_gate: peak_gate = abs(dg)
        if tk <= early_window_s and abs(dm) > early_peak: early_peak = abs(dm)
        if tk >= t_tail and abs(dm) > tail_mean: tail_mean = abs(dm)
        if abs(dm) > p["band_normal_m"]: over_band_dt += 1
        if k % rec_every == 0:
            out["t"].append(tk); out["dev_mean"].append(dm); out["dev_gate"].append(dg)
            out["dev_head"].append(h[0] - ht); out["gate"].append(ctrl.a)
            out["q_out"].append(q_gate)

    # 质量守恒：净入流 vs 池内蓄量增量。
    # 精确式 ∫[A(x)−A(ht)]dx = Δh·[b + m(x+ht)]（从断面几何严格推出），
    # 矩形下退化为 B0·Δh，与既有版本恒等；梯形下不要用 B(x)·Δh ——
    # 那只在 Δh→0 时成立，Δh≈335 mm 时会把守恒误差抬到 1.8%（是公式误差，不是求解器发散）。
    dv = 0.0
    for x in h:
        dv += (x - ht) * (sec.b + sec.m * (x + ht)) * dx
    return {"ncell": ncell, "dx_m": dx, "dt_s": dt, "steps": nt, "S0": S0,
            "ctrl_on": ctrl_on, "m_slope": m_slope, "sensor_len_m": sensor_len_m,
            "estimator": estimator, "est_track_err_mm": est_track_err * 1000,
            "nsens_cells": nsens, "early_window_s": early_window_s,
            "series": out,
            "peak_mean_mm": peak_mean * 1000, "peak_gate_mm": peak_gate * 1000,
            "early_peak_mm": early_peak * 1000,
            "final_mean_mm": tail_mean * 1000,
            "over_band_s": over_band_dt * dt,
            "vol_net_m3": vol_net, "vol_storage_m3": dv,
            "mass_balance_rel_err": abs(dv - vol_net) / max(abs(vol_net), 1.0)}


# ======================================================================
# 开关清单：只列**实测有效**的
# ======================================================================
SWITCHES = [
    ("S2a", "梯形断面 边坡 m=1.0", dict(m_slope=1.0)),
    ("S2b", "梯形断面 边坡 m=1.5", dict(m_slope=1.5)),
    ("S2c", "梯形断面 边坡 m=2.5", dict(m_slope=2.5)),
]

# 闸前传感器长度（m）：固定物理长度 ⇒ 回路增益与 dx 无关。
# 出处见 bm01_sve1d_v2_sensor.py：单点传感器 h[-1] 在 960 格上产生增幅极限环（峰值 212 mm 伪峰）。
SENSOR_LEN_M = 100.0


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 一维 v2：断面形状开关量化")
    ap.add_argument("--case", default=".")
    ap.add_argument("--grids", default="240,480")
    ap.add_argument("--only", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    case = os.path.abspath(args.case)
    with open(os.path.join(case, "raw", "bm01_main.json"), encoding="utf-8") as fh:
        raw = json.load(fh)
    p = raw["params"]
    step_dq = float(raw["scenario"]["step_dq"])
    rec = raw["metrics"]["pool1"]
    grids = [int(x) for x in args.grids.split(",")]
    t0 = time.time()

    c = math.sqrt(G * p["h_target"])
    print("BM-01 一维 v2　参数来自 raw/bm01_main.json（单一事实来源）")
    print(f"池长 {p['pool_len_m']:.0f} m × 底宽 {p['pool_width_m']:.0f} m，h_target={p['h_target']:.2f} m，"
          f"Q0={p['q0']:.0f} m³/s，阶跃 +{step_dq:.0f} m³/s，条带 ±{p['band_normal_m']*1000:.0f} mm")
    print(f"波速 c={c:.2f} m/s，L/c={p['pool_len_m']/c/60:.1f} min，2L/c={2*p['pool_len_m']/c/60:.1f} min\n")

    # ---- 自检：矩形 + 单点传感器 应逐行复现 v1 ----
    if args.selftest:
        v1ref = {"120": 142.88143604489179, "240": 156.6079038724011}
        print("[自检] 与 v1 记录值比对（矩形断面 + 单点传感器 = 历史行为）")
        ok = True
        for nc in (120, 240):
            r = sve_1d_v2(p, nc, step_dq, ctrl_on="gate", m_slope=0.0)
            d = abs(r["peak_mean_mm"] - v1ref[str(nc)])
            good = d < 0.5
            ok = ok and good
            print(f"   {nc:>4} 格  本文件 {r['peak_mean_mm']:.2f} mm vs "
                  f"v1 {v1ref[str(nc)]:.2f} mm　差 {d:.3f} mm → {'一致 ✓' if good else '不一致 ✗'}")
        print(f"   结论：{'矩形 + 单点传感器下本文件与 v1 逐行等价 ✓' if ok else '存在差异 ✗'}")
        return 0 if ok else 1

    # ---- 基准：矩形（结构与 v1 相同，传感器改为闸前固定长度平均）----
    print(f"[基准] 矩形断面 m=0；传感器 = 闸前 {SENSOR_LEN_M:.0f} m 平均（网格无关）")
    base = {}
    for nc in grids:
        base[nc] = sve_1d_v2(p, nc, step_dq, ctrl_on="gate_avg", m_slope=0.0,
                             sensor_len_m=SENSOR_LEN_M)
        r = base[nc]
        print(f"   {nc:>4} 格  dt={r['dt_s']:>6} s｜池均峰值 {r['peak_mean_mm']:7.1f} mm｜"
              f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰值 {r['peak_gate_mm']:7.1f} mm｜"
              f"末段 {r['final_mean_mm']:7.1f} mm｜守恒 {r['mass_balance_rel_err']:.1e}")
    base_conv = base[grids[-1]]["early_peak_mm"]
    spread = abs(base[grids[-1]]["early_peak_mm"] - base[grids[-2]]["early_peak_mm"])
    base_conv_full = base[grids[-1]]["peak_mean_mm"]
    print(f"   收敛峰值（**前 2h 窗**，取最细网格 {grids[-1]}）= {base_conv:.1f} mm"
          f"（与 {grids[-2]} 格差 {spread:.1f} mm，{spread/base_conv*100:.1f}%）"
          f"｜全程峰 {base_conv_full:.1f} mm（含晚期项，仅作参考）")

    # ---- 「旧表」口径复现留档（2026-10-01 补）----
    # 文档 README §9 / §4.2 与 summary.md 里有一列"旧表"值（170.3 / 116.2 / 102.4 / 82.0 mm）。
    # 那是**更正前三项口径叠加**的结果：单点传感器 h[-1] + 矩形反解纵坡 + 全程峰 + 120/240 两级
    # Richardson 外推。此前只有 legacy_rect_slope 留了其中"错纵坡"一项（且是最细网格收敛口径），
    # 读者拿证据包核不出那张变化量表。这里把旧口径整体跑出来存档，使旧值在包内可复现。
    def superseded_two_level(m_slope: float) -> dict:
        a = sve_1d_v2(p, 120, step_dq, ctrl_on="gate", m_slope=m_slope)["peak_mean_mm"]
        b = sve_1d_v2(p, 240, step_dq, ctrl_on="gate", m_slope=m_slope)["peak_mean_mm"]
        return {"peak_120_mm": round(a, 2), "peak_240_mm": round(b, 2),
                "extrapolated_2b_minus_a_mm": round(2.0 * b - a, 2),
                "note": "更正前口径复现（单点传感器 + 矩形反解纵坡 + 全程峰 + 两级外推）"}

    base_superseded = superseded_two_level(0.0)

    # ---- 逐开关 ----
    # 2026-10-01 修正（三项）：
    #  (1) 纵坡：梯形档此前沿用 v1 的**矩形反解**纵坡（V1.normal_slope）。矩形湿周 P=b+2h
    #      低估了梯形湿周 b+2h√(1+m²)，反解出的 S0 偏大（m=1/1.5/2.5 下偏大 1.38/1.52/1.79 倍）
    #      ⇒ 初值 (h_target, q0) 不是均匀流，净驱动力 S0−Sf≠0，模型自带假加速瞬变。
    #      现改用 section_normal_slope（断面感知），旧值在 legacy_rect_slope 留档。
    #  (2) 传感器：单点 h[-1] 的回路增益 ∝ 1/dx，960 格上出现增幅极限环（伪峰 212 mm）。
    #      现改用闸前固定长度平均（gate_avg），dx 无关；点传感器留档见 sensor_channel_study。
    #  (3) 峰值指标：弃用 Richardson 外推（序列非单调，外推无效），改报**最细网格收敛值**，
    #      并同时给出前 2 h 的 early_peak（物理瞬变窗）作交叉核对。
    rows = []
    only = args.only.split(",") if args.only else None
    for tag, name, kw in SWITCHES:
        if only and tag not in only:
            continue
        print(f"\n[{tag}] {name}")
        sec = Section(p["pool_width_m"], kw.get("m_slope", 0.0))
        s0_sec = section_normal_slope(sec, p["h_target"], p["q0"], MANNING_N)
        s0_leg = V1.normal_slope(p)
        print(f"   纵坡：断面感知 {s0_sec:.6e}（1/{1/s0_sec:,.0f}）｜"
              f"矩形反解 {s0_leg:.6e}（1/{1/s0_leg:,.0f}）"
              f"　矩形式偏大 {s0_leg/s0_sec:.2f} 倍")
        res, res_leg = {}, {}
        for nc in grids:
            r = sve_1d_v2(p, nc, step_dq, ctrl_on="gate_avg", s0=s0_sec,
                          sensor_len_m=SENSOR_LEN_M, **kw)
            res[nc] = r
            print(f"   {nc:>4} 格  dt={r['dt_s']:>6} s｜池均峰值 {r['peak_mean_mm']:7.1f} mm｜"
                  f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰值 {r['peak_gate_mm']:7.1f} mm｜"
                  f"末段 {r['final_mean_mm']:6.1f} mm｜越带 {r['over_band_s']/60:5.2f} min｜"
                  f"守恒 {r['mass_balance_rel_err']:.1e}")
            res_leg[nc] = sve_1d_v2(p, nc, step_dq, ctrl_on="gate_avg",
                                    sensor_len_m=SENSOR_LEN_M, **kw)
        conv = res[grids[-1]]["early_peak_mm"]
        conv_leg = res_leg[grids[-1]]["early_peak_mm"]
        conv_full = res[grids[-1]]["peak_mean_mm"]
        d_conv = conv - base_conv
        sp = abs(res[grids[-1]]["early_peak_mm"] - res[grids[-2]]["early_peak_mm"])
        print(f"   收敛峰值（前2h窗，{grids[-1]} 格）= {conv:.1f} mm　对基准 {base_conv:.1f} mm 的增量 "
              f"{d_conv:+.1f} mm（{d_conv/base_conv*100:+.1f}%）"
              f"　[网格差 {sp:.1f} mm / {sp/conv*100:.1f}%]")
        print(f"   全程峰 {conv_full:.1f} mm（含晚期项，仅作参考）")
        print(f"   [留档] 同一档改用矩形反解纵坡：收敛峰值 {conv_leg:.1f} mm"
              f"　→ 纵坡修正带来 {conv - conv_leg:+.1f} mm")
        st = superseded_two_level(kw.get("m_slope", 0.0))
        print(f"   [留档-旧表口径] 120格 {st['peak_120_mm']:.2f} / 240格 {st['peak_240_mm']:.2f} "
              f"→ 两级外推 {st['extrapolated_2b_minus_a_mm']:.1f} mm（= README §9「旧表」列）")
        rows.append({"tag": tag, "name": name, "switches": {k: v for k, v in kw.items()},
                     "slope_S0_section_aware": s0_sec,
                     "superseded_table": st,
                     "peak_mean_mm": {str(nc): res[nc]["peak_mean_mm"] for nc in grids},
                     "early_peak_mm": {str(nc): res[nc]["early_peak_mm"] for nc in grids},
                     "peak_gate_mm": {str(nc): res[nc]["peak_gate_mm"] for nc in grids},
                     "final_mean_mm": {str(nc): res[nc]["final_mean_mm"] for nc in grids},
                     "over_band_s": {str(nc): res[nc]["over_band_s"] for nc in grids},
                     "mass_balance_rel_err": {str(nc): res[nc]["mass_balance_rel_err"] for nc in grids},
                     "converged_mm": conv, "converged_grid": grids[-1],
                     "converged_metric": "early_peak_mm（前 2h 窗）",
                     "converged_fullrun_mm": conv_full,
                     "grid_spread_mm": sp, "delta_vs_base_mm": d_conv,
                     "delta_pct": d_conv / base_conv * 100.0,
                     "legacy_rect_slope": {
                         "slope_S0": s0_leg, "ratio_over_section_aware": s0_leg / s0_sec,
                         "peak_mean_mm": {str(nc): res_leg[nc]["peak_mean_mm"] for nc in grids},
                         "early_peak_mm": {str(nc): res_leg[nc]["early_peak_mm"] for nc in grids},
                         "peak_gate_mm": {str(nc): res_leg[nc]["peak_gate_mm"] for nc in grids},
                         "final_mean_mm": {str(nc): res_leg[nc]["final_mean_mm"] for nc in grids},
                         "over_band_s": {str(nc): res_leg[nc]["over_band_s"] for nc in grids},
                         "converged_mm": conv_leg,
                         "note": "更正前的记录（矩形反解纵坡）；保留用于说明错纵坡的影响，不作结论依据"}})

    out_raw = os.path.join(case, "raw")
    payload = {
        "case_id": "BM-01", "purpose": "1d_model_fidelity_upgrade",
        "evidence_grade": "synthetic-1d-sve",
        "note": "断面形状开关逐项量化；不改变 BM-01 证据等级；量级示例参数，非现场结论。"
                "2026-10-01 修正：梯形档此前沿用矩形反解纵坡（偏大 1.38~1.79 倍），"
                "已改用断面感知纵坡（section_normal_slope），旧值在 legacy_rect_slope 留档。"
                "传感器由单点 h[-1] 改为闸前固定长度平均（gate_avg），旧行为见 sensor_channel_study。"
                "峰值指标弃用 Richardson 外推（序列非单调），改报最细网格收敛值 converged_mm。",
        "solver": {"base": "v1 conservative Rusanov + explicit Euler",
                   "deps": "python stdlib only", "manning_n": MANNING_N, "cfl": CFL,
                   "slope_S0_rect_legacy": V1.normal_slope(p),
                   "slope_S0_note": "矩形基准档（m=0）两种反解恒等，故基准不受纵坡修正影响；"
                                    "梯形档实际使用的纵坡见各 switch_results 的 slope_S0_section_aware",
                   "sensor_channel": "gate_avg",
                   "sensor_len_m": SENSOR_LEN_M,
                   "sensor_note": "闸前固定物理长度内取平均水位。单点 h[-1] 的回路增益 ∝ 1/dx，"
                                  "960 格上产生增幅极限环（伪峰 212 mm、闸前 527 mm）；"
                                  "详见 bm01_sve1d_v2_sensor.py 与 sensor_channel_study。",
                   "peak_metric": "converged_mm = 最细网格峰均；early_peak_mm = 前 2 h 内峰均（物理瞬变窗）",
                   "wave_celevity_mps": c,
                   "travel_time_min": p["pool_len_m"] / c / 60},
        "params_source": "raw/bm01_main.json",
        "grids": grids,
        "baseline_rect": {"peak_mean_mm": {str(nc): base[nc]["peak_mean_mm"] for nc in grids},
                          "early_peak_mm": {str(nc): base[nc]["early_peak_mm"] for nc in grids},
                          "peak_gate_mm": {str(nc): base[nc]["peak_gate_mm"] for nc in grids},
                          "converged_mm": base_conv, "converged_grid": grids[-1],
                          "converged_metric": "early_peak_mm（前 2h 窗）",
                          "converged_fullrun_mm": base_conv_full,
                          "grid_spread_mm": spread,
                          "superseded_table": base_superseded,
                          "extrapolated_mm_legacy": 2.0 * base[grids[-1]]["peak_mean_mm"]
                                                    - base[grids[-2]]["peak_mean_mm"] if len(grids) >= 2 else None},
        "sensor_channel_study": {
            "note": "单点传感器 h[-1] 的回路增益 ∝ 1/dx，细网格上出现晚期增幅极限环"
                    "（矩形 960 格全程峰 212 mm、闸前 527 mm；早峰仍为 158.7 mm 稳定）。"
                    "本包已改用闸前固定长度平均（gate_avg）。完整网格扫描见 "
                    "raw/bm01_sve1d_v2_sensor.json（bm01_sve1d_v2_sensor.py）。",
            "root_cause": "控制器测点/动作点都在最末单元；网格越细该单元蓄量越小，"
                          "回路增益上升，300 s 采样周期下 PI 回路转为增幅极限环。",
            "impact": "只影响**晚期**响应（t>150 min）；物理峰出现在 t≈36 min，"
                      "用前 2h 窗（early_peak_mm）作收敛指标即不受影响。",
        },
        "switch_results": rows,
        "disabled_switches": {
            "R_hydraulic_full_rect": {
                "switch": "完整水力半径 R=A/P（矩形）",
                "reason": "m=0 时 P=b+2h 与 v1 的 (B+2h) 恒等，增量恒为 0——同一个式子，不列入成果"},
            "gate_internal": {
                "switch": "闸门改内部建筑物",
                "reason": "闸室单元占全池水面 1/120，本参数下影响 <0.5 mm 不可辨；"
                          "要有实质区别需闸门局部收缩（侧收缩/胸墙/淹没出流）等真实结构，留待有断面资料时再做"},
        },
        "lumped_recorded_mm": rec["peak_dev_m"] * 1000,
        "runtime_s": round(time.time() - t0, 2),
    }
    jp = os.path.join(out_raw, "bm01_sve1d_v2.json")
    with open(jp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    print(f"\n总用时 {time.time()-t0:.1f} s")
    print(f"输出：{jp}")
    print("证据等级：synthetic-1d-sve（不改变 BM-01 判定码，非现场结论）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
