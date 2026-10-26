# -*- coding: utf-8 -*-
"""生成《BM-01 渠池前馈控制 阶段汇报》用的插图和锚点数字。

只读既有模型 + 既有 raw 数据，不改动任何案例结论。
图宽统一按 15 cm 出（figsize 宽 = 5.9055 in），字号直接按最终尺寸设，
这样 docx 转换器按「像素 ÷ dpi」插入时不会把字缩到看不清。

运行：
  PYTHONPATH=<venv>/pylibs-extra python _make_report_figs.py
"""
from __future__ import annotations

import json
import math
import os
import sys
from dataclasses import replace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Patch
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from bm01_pool_step import Params, build_scenarios          # noqa: E402
from bm01_single_pool_ff import id_tau, run, run_delayed     # noqa: E402
from bm01_single_pool_ff import stats                        # noqa: E402

OUT = os.path.join(HERE, "figures", "report")
os.makedirs(OUT, exist_ok=True)

W_IN = 15.0 / 2.54                      # 15 cm 画布宽
TARGET_CM = 15.0

_avail = {f.name for f in font_manager.fontManager.ttflist}
for _c in ["Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC"]:
    if _c in _avail:
        plt.rcParams["font.sans-serif"] = [_c]
        break
plt.rcParams["axes.unicode_minus"] = False

INK, MUTED, HAIR = "#2c2c2a", "#5f5e5a", "#d3d1c7"
BLUE, AMBER, RED, GREEN = "#185fa5", "#ba7517", "#a32d2d", "#1a7f37"
PANEL, WATER = "#f4f4f2", "#b5d4f4"


def save(fig, name: str) -> str:
    """存 PNG 并把 dpi 元数据改写成「自然宽度 = 15 cm」。"""
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=200, facecolor="white")
    plt.close(fig)
    from PIL import Image
    im = Image.open(p)
    w, _ = im.size
    im.save(p, dpi=(w * 2.54 / TARGET_CM,) * 2)
    return p


# ======================================================================
# 锚点数字
# ======================================================================
def collect() -> dict:
    base = Params()
    sc = next(s for s in build_scenarios(base) if s.key == "main")
    tau = id_tau(base, "dynamic")

    A = base.pool_width_m * base.h_target
    v = base.q0 / A
    c_grav = math.sqrt(9.81 * A / base.pool_width_m)
    c_dyn = v + c_grav

    d = {}
    d["geometry"] = {
        "pool_len_m": base.pool_len_m, "pool_width_m": base.pool_width_m,
        "area_m2": base.area_m2, "h_target_m": base.h_target,
        "q0_m3s": base.q0, "a0": base.a0, "kg": base.kg,
        "gate_up_capacity_m3s": base.gate_up_capacity,
        "step_dq_m3s": sc.step_dq,
        "ctrl_period_s": base.ctrl_period_s,
        "rate_limit_per_s": base.rate_limit_per_s,
        "full_stroke_s": 1.0 / base.rate_limit_per_s,
        "band_normal_mm": base.band_normal_m * 1000,
    }
    try:
        tau_kin_h: float | None = id_tau(base, "gate") / 3600.0
    except Exception:                                        # noqa: BLE001
        tau_kin_h = None
    d["tau"] = {"v_flow": v, "c_gravity_wave": c_grav, "c_dynamic": c_dyn,
                "tau_s": tau, "tau_kinematic_h": tau_kin_h,
                "tau_kin_label": (f"{tau_kin_h:.1f} h" if tau_kin_h else "算不出")}

    # 自由上涨 / 开到位时间
    d["scales"] = {
        "free_rise_mm_per_ctrl": sc.step_dq * base.ctrl_period_s / base.area_m2 * 1000,
        "ideal_ff_open": sc.step_dq / base.kg,
        "t_act_s": (sc.step_dq / base.kg) / base.rate_limit_per_s,
        "tau_over_ctrl": tau / base.ctrl_period_s,
        "tau_over_tact": tau / ((sc.step_dq / base.kg) / base.rate_limit_per_s),
    }

    # ---------- 集总 vs ID ----------
    lump_off = stats(base, run(base, sc, 0.0, 0.0)["series"]["h1_dev"],
                     run(base, sc, 0.0, 0.0)["series"]["a1"])
    r_lump_ff = run(base, sc, 1.0, 0.0)
    lump_ff = stats(base, r_lump_ff["series"]["h1_dev"], r_lump_ff["series"]["a1"])

    id_off = run_delayed(base, sc, tau, 0.0, 0.0)
    s_id_off = stats(base, id_off["series"]["h1_dev"], id_off["series"]["a1"])
    id_ff0 = run_delayed(base, sc, tau, 1.0, 0.0)
    s_id_ff0 = stats(base, id_ff0["series"]["h1_dev"], id_ff0["series"]["a1"])
    id_ff382 = run_delayed(base, sc, tau, 1.0, 382.0)
    s_id_ff382 = stats(base, id_ff382["series"]["h1_dev"], id_ff382["series"]["a1"])

    d["lumped"] = {"off_peak_mm": lump_off["peak"], "ff_peak_mm": lump_ff["peak"]}
    d["id_model"] = {
        "off_peak_mm": s_id_off["peak"],
        "ff_delay0_peak_mm": s_id_ff0["peak"],
        "ff_best_peak_mm": s_id_ff382["peak"],
        "ff_gain_pct": 100 * (1 - s_id_ff382["peak"] / s_id_off["peak"]),
        "ff_delay0_gain_pct": 100 * (1 - s_id_ff0["peak"] / s_id_off["peak"]),
        "margin_before_mm": base.band_normal_m * 1000 - s_id_off["peak"],
        "margin_after_mm": base.band_normal_m * 1000 - s_id_ff382["peak"],
    }

    # ---------- 瓶颈分解（都在 ID 上做） ----------
    br = {}
    r = run_delayed(replace(base, ctrl_period_s=5.0), sc, tau, 0.0, 0.0)
    br["fast_ctrl_only"] = stats(base, r["series"]["h1_dev"], r["series"]["a1"])["peak"]
    r = run_delayed(replace(base, rate_limit_per_s=1.0), sc, tau, 0.0, 0.0)
    br["no_rate_limit_only"] = stats(base, r["series"]["h1_dev"], r["series"]["a1"])["peak"]
    r = run_delayed(replace(base, rate_limit_per_s=1.0), sc, tau, 1.0, 382.0)
    br["ff_plus_no_limit"] = stats(base, r["series"]["h1_dev"], r["series"]["a1"])["peak"]
    d["bottleneck"] = br

    # ---------- 曲线：ID 机理三联图 ----------
    curves = {
        "t_min": id_off["series"]["t"] / 60.0,
        "q_up": id_off["series"]["q_upstream"],
        "q_tail": id_off["series"]["q_tail_in"],
        "h_lumped_mm": run(base, sc, 0.0, 0.0)["series"]["h1_dev"] * 1000,
        "h_id_off_mm": id_off["series"]["h1_dev"] * 1000,
        "h_id_ff_mm": id_ff382["series"]["h1_dev"] * 1000,
        "a_id_off": id_off["series"]["a1"],
        "a_id_ff": id_ff382["series"]["a1"],
        "t_step_min": base.t_step_s / 60.0,
        "t_tail_start_min": (base.t_step_s + tau) / 60.0,
    }
    np.savez(os.path.join(HERE, "raw", "bm01_report_curves.npz"), **curves)
    d["t_step_min"] = base.t_step_s / 60.0
    d["t_tail_start_min"] = curves["t_tail_start_min"]

    with open(os.path.join(HERE, "raw", "bm01_report_anchors.json"), "w",
              encoding="utf-8") as fh:
        json.dump(d, fh, ensure_ascii=False, indent=2)
    return d


# ======================================================================
# 图 1 场景与信息流
# ======================================================================
def fig_scene(d: dict) -> str:
    fig, a = plt.subplots(figsize=(W_IN, 2.78))
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.axis("off")

    bed, surf = 0.300, 0.640
    x0, x1 = 0.185, 0.760

    # 渠池水体
    a.add_patch(Rectangle((x0, bed), x1 - x0, surf - bed,
                          facecolor=WATER, edgecolor="#7fa9d4", lw=0.8, zorder=1))
    a.plot([x0, x1], [surf, surf], color=BLUE, lw=1.5, zorder=3)
    a.plot([x0, x1], [bed, bed], color="#8a8577", lw=2.0, zorder=2)

    # 闸门
    for xg, lab in ((x0 - 0.012, "上游闸门\n（来水）"), (x1 + 0.012, "池尾闸门\n（被控）")):
        a.add_patch(Rectangle((xg - 0.007, bed - 0.02), 0.014, 0.49,
                              facecolor="#5f5e5a", edgecolor="none", zorder=4))
        a.text(xg, bed - 0.035, lab, ha="center", va="top",
               fontsize=6.0, color=MUTED, linespacing=1.45)

    # 来水箭头（标注移到箭头下方，避开池首流量计）
    a.annotate("", xy=(x0 - 0.020, surf - 0.060), xytext=(x0 - 0.078, surf - 0.060),
               arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.3))
    a.text(x0 - 0.026, surf - 0.018, "来水\n180 → 225 m³/s", ha="right", va="bottom",
           fontsize=6.2, color=BLUE, linespacing=1.45)

    # 池首流量计
    xs = x0 + 0.040
    a.plot([xs], [surf], marker="v", color=AMBER, ms=5.5, zorder=6)
    a.text(xs, surf + 0.026, "池首流量计", ha="center", va="bottom",
           fontsize=6.2, color=AMBER)

    # 池尾水位计
    xh = x1 - 0.032
    a.plot([xh], [surf], marker="v", color=GREEN, ms=5.5, zorder=6)
    a.text(xh, surf + 0.026, "池尾水位计", ha="center", va="bottom",
           fontsize=6.2, color=GREEN)

    # 池长标注（放在渠底与控制器之间的空档）
    a.annotate("", xy=(x0, bed - 0.055), xytext=(x1, bed - 0.055),
               arrowprops=dict(arrowstyle="<|-|>", color=MUTED, lw=0.8))
    a.text((x0 + x1) / 2, bed - 0.072, "池长 6000 m｜水面宽 50 m｜目标水位 7.00 m",
           ha="center", va="top", fontsize=6.2, color=MUTED)

    # 控制器
    a.add_patch(Rectangle((0.300, 0.012), 0.340, 0.098,
                          facecolor=PANEL, edgecolor="#b4b2a9", lw=0.8, zorder=3))
    a.text(0.470, 0.078, "下游闸门控制器   PI 反馈 + 前馈", ha="center", va="center",
           fontsize=6.4, color=INK, zorder=4)
    a.text(0.470, 0.040, "每 300 s 下发一次命令", ha="center", va="center",
           fontsize=5.8, color=MUTED, zorder=4)

    arr = dict(arrowstyle="-|>", lw=1.0, shrinkA=0, shrinkB=0)
    # 池首流量计 -> 控制器（前馈通道）
    a.add_patch(FancyArrowPatch((xs, surf - 0.030), (0.330, 0.108),
                                connectionstyle="arc3,rad=-0.25",
                                color=AMBER, **arr))
    a.text(0.222, 0.415, "前馈通道\n零信息延迟", ha="center", va="bottom",
           fontsize=6.0, color=AMBER, linespacing=1.45)
    # 池尾水位计 -> 控制器（反馈通道）
    a.add_patch(FancyArrowPatch((xh, surf - 0.030), (0.610, 0.108),
                                connectionstyle="arc3,rad=0.25",
                                color=GREEN, **arr))
    a.text(0.700, 0.415, "反馈通道\n水位涨了才动", ha="center", va="bottom",
           fontsize=6.0, color=GREEN, linespacing=1.45)
    # 控制器 -> 闸门
    a.add_patch(FancyArrowPatch((0.640, 0.062), (x1 + 0.012, bed - 0.020),
                                connectionstyle="arc3,rad=-0.25",
                                color=INK, **arr))
    a.text(0.998, 0.106, "开度命令\n限速：满行程 ≈50 min", ha="right", va="top",
           fontsize=6.0, color=INK, linespacing=1.45)

    # 标题行
    a.text(0.0, 1.00, "对象：一个池、一座池尾闸门；传感器在池的首端，闸门在池的尾端",
           ha="left", va="top", fontsize=7.2, color=INK)
    a.text(0.0, 0.917, d["scene_note"], ha="left", va="top",
           fontsize=6.2, color=BLUE)

    fig.tight_layout(pad=0.15)
    return save(fig, "r01_scene.png")


# ======================================================================
# 图 2 降阶阶梯
# ======================================================================
def fig_ladder(d: dict) -> str:
    rows = [
        ("完整圣维南方程", "连续方程 + 动量方程（含惯性项）",
         "什么都不丢。代价：要解偏微分方程，要全渠段地形、糙率、边坡", INK, "white"),
        ("扩散波", "动量方程丢掉惯性项，保留压力项与摩阻项",
         "丢掉加速度。急变过程（闸门瞬间启闭）失真", MUTED, "white"),
        ("运动波", "再丢掉压力项，波速 = 流量对水位的导数 ÷ 水面宽",
         f"丢掉回水影响。BM-01 算出 {d['tau']['tau_kinematic_h']:.1f} h，比仿真时长 4 h 还长 → 排除", MUTED, "white"),
        ("积分时滞（ID）", "纯延迟一个传播时间 + 一次水量平衡积分",
         "丢掉池内水面形态，只留一个滞后时间。本汇报用的就是这一级", BLUE, "#eef5fc"),
        ("集总 0D", "传播时间取 0，全池只有一个水位",
         "丢掉相位。「池首」和「池尾」被压成同一个点", MUTED, "white"),
    ]
    fig, a = plt.subplots(figsize=(W_IN, 3.05))
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.axis("off")

    y = 0.855
    h = 0.145
    for i, (name, keep, lose, col, bg) in enumerate(rows):
        hl = (name.startswith("积分时滞"))
        a.add_patch(Rectangle((0.115, y - h), 0.865, h,
                              facecolor=bg, edgecolor=col, lw=1.1 if hl else 0.7,
                              zorder=2))
        a.text(0.128, y - h / 2, name, ha="left", va="center",
               fontsize=7.0, color=col, fontweight="bold", zorder=3)
        a.text(0.325, y - h * 0.30, keep, ha="left", va="center",
               fontsize=6.2, color=INK, zorder=3)
        a.text(0.325, y - h * 0.72, lose, ha="left", va="center",
               fontsize=5.9, color=MUTED, zorder=3)
        if i < len(rows) - 1:
            a.annotate("", xy=(0.048, y - h - 0.030), xytext=(0.048, y - h - 0.002),
                       arrowprops=dict(arrowstyle="-|>", color=HAIR, lw=1.0))
        y -= h + 0.028

    a.annotate("", xy=(0.048, 0.020), xytext=(0.048, 0.855),
               arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.1))
    a.text(0.030, 0.44, "逐级降阶", ha="right", va="center", rotation=90,
           fontsize=7.0, color=MUTED)
    a.text(0.115, 0.995, "越往下：计算量越小，能回答的问题也越少", ha="left",
           va="top", fontsize=7.2, color=INK)
    a.text(0.980, 0.995, "高亮这一级是本汇报的工作面", ha="right", va="top",
           fontsize=6.3, color=BLUE)

    fig.tight_layout(pad=0.15)
    return save(fig, "r02_ladder.png")


# ======================================================================
# 图 3 ID 机理三联
# ======================================================================
def fig_id_mech(d: dict) -> str:
    z = np.load(os.path.join(HERE, "raw", "bm01_report_curves.npz"))
    t = z["t_min"]
    tau = d["tau"]["tau_s"]
    tstep = d["t_step_min"]
    band = d["geometry"]["band_normal_mm"]
    t_tail = tstep + tau / 60.0

    fig, ax = plt.subplots(3, 1, figsize=(W_IN, 4.75),
                           gridspec_kw=dict(height_ratios=[0.80, 1.20, 0.68]))
    fig.subplots_adjust(hspace=0.60, top=0.895, bottom=0.095,
                        left=0.088, right=0.985)
    fig.suptitle(f"积分时滞模型（ID）怎么算的｜传播时间 = {tau:.0f} s",
                 fontsize=7.6, fontweight="bold", y=0.985)

    # (a) 纯延迟
    a = ax[0]
    a.axvspan(tstep, t_tail, color="#fdf3e3", zorder=0)
    a.plot(t, z["q_up"], color=AMBER, lw=1.4, label="池首流量计此刻读到的来水")
    a.plot(t, z["q_tail"], color=BLUE, lw=1.4, ls="--",
           label="池尾实际收到的入流（同一波形，右移一个传播时间）")
    a.annotate("", xy=(t_tail, 212), xytext=(tstep, 212),
               arrowprops=dict(arrowstyle="<|-|>", color=RED, lw=1.0))
    a.text((tstep + t_tail) / 2, 215, f"纯延迟 {tau:.0f} s", ha="center", va="bottom",
           fontsize=6.2, color=RED)
    a.set_ylabel("流量 (m³/s)", fontsize=6.8, color=MUTED)
    a.set_ylim(168, 242)
    a.set_title("(a) 第一半：入流被原样延迟，波形不变、只往后挪",
                fontsize=6.8, color=INK, loc="left", pad=3)
    a.legend(fontsize=6.0, frameon=False, loc="lower right")

    # (b) 水位
    b = ax[1]
    b.axhspan(-band, band, color="#e8f3ff", zorder=0)
    for yy in (band, -band):
        b.axhline(yy, color=GREEN, lw=0.9)
    b.axvline(tstep, color="#888888", lw=0.8, ls="--")
    b.axvline(t_tail, color=BLUE, lw=0.8, ls=":")
    b.plot(t, z["h_lumped_mm"], color=MUTED, lw=1.0, ls="-.",
           label=f"集总 0D（传播时间=0，池首池尾同一个水位）峰值 {d['lumped']['off_peak_mm']:.1f} mm")
    b.plot(t, z["h_id_off_mm"], color=RED, lw=1.3,
           label=f"ID 纯反馈 峰值 {d['id_model']['off_peak_mm']:.1f} mm")
    b.plot(t, z["h_id_ff_mm"], color=BLUE, lw=1.5,
           label=f"ID 前馈（提前 382 s）峰值 {d['id_model']['ff_best_peak_mm']:.1f} mm")
    b.text(t[-1], band + 5, f"正常调节条带 ±{band:.0f} mm", fontsize=6.0, color=GREEN,
           ha="right", va="bottom")
    b.annotate("阶跃发生", xy=(tstep, 150), xytext=(2, 172),
               fontsize=6.0, color="#777777", ha="left",
               arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.7))
    b.annotate("水到池尾", xy=(t_tail, 150), xytext=(32, 172),
               fontsize=6.0, color=BLUE, ha="left",
               arrowprops=dict(arrowstyle="-", color="#7fa9d4", lw=0.7))
    b.set_ylabel("水位偏差 (mm)", fontsize=6.8, color=MUTED)
    b.set_ylim(-185, 186)
    b.set_title("(b) 第二半：把延迟后的入流减掉出流，再对时间做一次积分，就是水位",
                fontsize=6.8, color=INK, loc="left", pad=3)
    b.legend(fontsize=6.0, frameon=False, loc="lower left", ncol=1)

    # (c) 闸门开度
    c = ax[2]
    c.axvline(tstep, color="#888888", lw=0.8, ls="--")
    c.plot(t, z["a_id_off"], color=RED, lw=1.2, label="纯反馈：水到才动")
    c.plot(t, z["a_id_ff"], color=BLUE, lw=1.3, label="前馈：与水到同步地动")
    c.set_ylabel("闸门开度", fontsize=6.8, color=MUTED)
    c.set_xlabel("时间 (min)", fontsize=6.8, color=MUTED)
    c.set_title("(c) 动作差在哪：前馈把闸门提前到一个斜坡上，不是更早、也不是更晚",
                fontsize=6.8, color=INK, loc="left", pad=3)
    c.legend(fontsize=6.0, frameon=False, loc="lower right")

    for a_ in ax:
        a_.tick_params(labelsize=6.2, colors=MUTED)
        a_.set_xlim(0, t[-1])
        for s in ("top", "right"):
            a_.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a_.spines[s].set_color(HAIR)
        a_.grid(color="#efeeea", lw=0.6)
        a_.set_axisbelow(True)

    return save(fig, "r03_id_mech.png")


# ======================================================================
# 图 4 波速分解
# ======================================================================
def fig_wavespeed(d: dict) -> str:
    f = d["tau"]
    fig, a = plt.subplots(figsize=(W_IN, 1.78))
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.axis("off")

    boxes = [(0.015, f"断面平均流速\n{f['v_flow']:.3f} m/s", BLUE,
              "流量 ÷ 过水断面面积"),
             (0.345, f"重力波项\n{f['c_gravity_wave']:.3f} m/s", "#6f9fd8",
              "水深与断面形状决定"),
             (0.685, f"动力波速\n{f['c_dynamic']:.3f} m/s", INK,
              "两者相加，不需要率定")]
    w = 0.285
    for x, lab, col, sub in boxes:
        a.add_patch(Rectangle((x, 0.56), w, 0.30, facecolor=col,
                              edgecolor="none", zorder=3))
        a.text(x + w / 2, 0.735, lab, ha="center", va="center", fontsize=6.4,
               color="white", linespacing=1.5, zorder=4)
        a.text(x + w / 2, 0.515, sub, ha="center", va="top", fontsize=5.8,
               color=MUTED, zorder=4)
    a.text(0.300 + 0.015, 0.71, "+", ha="center", va="center", fontsize=10,
           color=MUTED)
    a.text(0.630 + 0.015, 0.71, "=", ha="center", va="center", fontsize=10,
           color=MUTED)

    a.text(0.5, 0.30, f"传播时间 = 池长 {d['geometry']['pool_len_m']:.0f} m ÷ "
                      f"{f['c_dynamic']:.3f} m/s = {f['tau_s']:.0f} s"
                      f"（约 {f['tau_s'] / 60:.1f} min）",
           ha="center", va="center", fontsize=7.6, color=RED)
    a.text(0.5, 0.10, "波形不是慢慢淡掉，是原样搬过去——所以模型里用纯延迟，不用扩散方程",
           ha="center", va="center", fontsize=6.2, color=MUTED)

    return save(fig, "r04_wavespeed.png")


# ======================================================================
# 图 5 控制周期
# ======================================================================
def fig_ctrl_period(d: dict) -> str:
    """控制周期的物理意义：控制器每 300 s 只拿到一个读数，两次读数之间水位无人管。

    曲线取真实的 ID 纯反馈时序（阶跃在 900 s 施加），只看前 1300 s。
    """
    g = d["geometry"]
    z = np.load(os.path.join(HERE, "raw", "bm01_report_curves.npz"))
    t_all = z["t_min"] * 60.0                   # s
    h_all = z["h_id_off_mm"]
    sel = t_all <= 4200.0
    t, h = t_all[sel], h_all[sel]

    step = g["ctrl_period_s"]
    idx = [int(round(c / 5.0)) for c in np.arange(0.0, 4200.0, step)]
    pairs = [(h[idx[k + 1]] - h[idx[k]], k) for k in range(len(idx) - 1)]
    jump, k0 = max(pairs)
    k1 = k0 + 1

    fig, a = plt.subplots(figsize=(W_IN, 2.15))
    fig.subplots_adjust(top=0.845, bottom=0.175, left=0.098, right=0.985)
    a.axhline(0, color=HAIR, lw=0.8)
    a.plot(t, h, color=BLUE, lw=1.5, label="池尾水位偏差（真实仿真：ID 纯反馈）")
    for c in np.arange(0.0, 4200.0, step):
        i = int(round(c / 5.0))
        a.axvline(c, color=AMBER, lw=0.8, ls=":", zorder=1)
        a.plot([c], [h[i]], "o", color=AMBER, ms=3.6, zorder=5)
    a.plot([], [], "o--", color=AMBER, ms=3.6, lw=0.8,
           label=f"控制器每 {step:.0f} s 才读一次水位、改一次命令")

    a.annotate("", xy=(t[idx[k1]], h[idx[k1]]), xytext=(t[idx[k0]], h[idx[k0]]),
               arrowprops=dict(arrowstyle="<|-|>", color=RED, lw=1.2))
    a.text((t[idx[k0]] + t[idx[k1]]) / 2, h[idx[k1]] + 3.0,
           f"这一拍 {step:.0f} s 里水位走了 {jump:.0f} mm\n期间没有任何人管它",
           ha="center", va="bottom", fontsize=6.2, color=RED, linespacing=1.5)

    a.set_xlabel("时间 (s)（阶跃在 900 s 施加，水到池尾在 1582 s）",
                 fontsize=6.8, color=MUTED)
    a.set_ylabel("水位偏差 (mm)", fontsize=6.8, color=MUTED)
    a.set_xlim(0, 4200)
    a.set_ylim(-8, max(float(h.max()) * 1.55, 12))
    a.set_title("控制周期的物理意义：两次看水位、改命令之间的间隔；间隔里的水位没人管",
                fontsize=7.0, color=INK, loc="left", pad=6)
    a.legend(fontsize=6.0, frameon=False, loc="upper left")
    a.tick_params(labelsize=6.2, colors=MUTED)
    for s_ in ("top", "right"):
        a.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        a.spines[s_].set_color(HAIR)
    a.grid(color="#efeeea", lw=0.6)
    a.set_axisbelow(True)
    return save(fig, "r05_ctrl_period.png")


# ======================================================================
# 图 6 真实来水下三种取法
# ======================================================================
def fig_real_inflow(d: dict) -> str:
    z = np.load(os.path.join(HERE, "raw", "bm01_inflow_integral.npz"))
    t = z["t"] / 60.0
    band = d["geometry"]["band_normal_mm"]
    off = np.abs(z["h_off"]).max() * 1000
    inst = np.abs(z["h_inst"]).max() * 1000
    win = np.abs(z["h_win"]).max() * 1000

    fig, ax = plt.subplots(2, 1, figsize=(W_IN, 3.35),
                           gridspec_kw=dict(height_ratios=[0.60, 1.0]))
    fig.subplots_adjust(hspace=0.50, top=0.875, bottom=0.135,
                        left=0.088, right=0.985)
    fig.suptitle("来水未知（含仪表噪声）时：三种取法的实测对比",
                 fontsize=7.6, fontweight="bold", y=0.985)

    a = ax[0]
    a.plot(t, z["q_true"], color=MUTED, lw=1.0, label="真实来水")
    a.plot(t, z["q_meas"], color=AMBER, lw=0.7, alpha=0.85, label="池首流量计读数（含噪声）")
    a.set_ylabel("流量 (m³/s)", fontsize=6.8, color=MUTED)
    a.set_title("(a) 池首测到的从来不是一条干净的阶跃", fontsize=6.8,
                color=INK, loc="left", pad=3)
    a.legend(fontsize=6.0, frameon=False, loc="upper right", ncol=2)

    b = ax[1]
    b.plot(t, z["h_off"] * 1000, color=RED, lw=1.2, label=f"纯反馈 峰值 {off:.1f} mm")
    b.plot(t, z["h_inst"] * 1000, color=AMBER, lw=1.1, ls="--",
           label=f"瞬时前馈 峰值 {inst:.1f} mm（追噪声，闸门动作最费）")
    b.plot(t, z["h_win"] * 1000, color=BLUE, lw=1.6,
           label=f"窗口积分前馈 峰值 {win:.1f} mm（降 {100 * (1 - win / off):.0f}%）")
    ymin = min(float((z[k] * 1000).min()) for k in ("h_off", "h_inst", "h_win"))
    ymax = max(float((z[k] * 1000).max()) for k in ("h_off", "h_inst", "h_win"))
    pad_y = 0.30 * max(ymax - ymin, 1.0)
    b.set_ylim(ymin - pad_y, ymax + pad_y)
    b.text(0.5, 0.02, f"（条带为 ±{band:.0f} mm，三种取法全程都在带内；此图已放大看差异）",
           transform=b.transAxes, fontsize=6.0, color=MUTED, va="bottom")
    b.set_ylabel("水位偏差 (mm)", fontsize=6.8, color=MUTED)
    b.set_xlabel("时间 (min)", fontsize=6.8, color=MUTED)
    b.set_title("(b) 把过去一个传播时间的入流积出来当指令，比追瞬时值稳得多",
                fontsize=6.8, color=INK, loc="left", pad=3)
    b.legend(fontsize=6.0, frameon=False, loc="upper right")

    for a_ in ax:
        a_.tick_params(labelsize=6.2, colors=MUTED)
        a_.set_xlim(0, t[-1])
        for s_ in ("top", "right"):
            a_.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"):
            a_.spines[s_].set_color(HAIR)
        a_.grid(color="#efeeea", lw=0.6)
        a_.set_axisbelow(True)

    fig.tight_layout(pad=0.15)
    return save(fig, "r06_real_inflow.png")


# ======================================================================
# 图 7 前馈收益 vs 传播时间
# ======================================================================
def fig_tau_peak(d: dict) -> str:
    z = np.load(os.path.join(HERE, "raw", "bm01_ff_tau_map.npz"))
    tc = d["geometry"]["ctrl_period_s"]
    x = z["tau"] / tc
    off = z["off"]
    g_win = 100 * (off - z["win"]) / off
    g_inst = 100 * (off - z["inst"]) / off
    xp = d["scales"]["tau_over_tact"]

    fig, a = plt.subplots(figsize=(W_IN, 3.30))
    fig.subplots_adjust(top=0.695, bottom=0.155, left=0.098, right=0.985)
    a.axhline(0, color=RED, lw=1.0, ls="--")
    a.axvspan(0, 0.4, color="#fbecec", zorder=0)
    a.text(0.06, -92, "负收益区\n装了不如不装", fontsize=6.0, color=RED,
           va="bottom", linespacing=1.4)
    a.axvline(xp, color=GREEN, lw=0.9, ls=":")

    a.plot(x, g_inst, color=AMBER, lw=1.2, marker="o", ms=2.6,
           label="瞬时前馈：来水一变就照办")
    a.plot(x, g_win, color=BLUE, lw=1.6, marker="s", ms=2.6,
           label="窗口积分前馈：积过去一个传播时间的入流")
    a.annotate("峰值落在传播时间 ≈ 闸门开到位时间处\n（本例 1.8 倍控制周期）",
               xy=(float(xp), float(np.interp(xp, x, g_win))),
               xytext=(3.6, 92), fontsize=6.0, color=GREEN,
               va="top", linespacing=1.5,
               arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=0.8))
    a.set_xlabel("传播时间 ÷ 控制周期（倍）", fontsize=6.8, color=MUTED)
    a.set_ylabel("相比纯反馈的峰值降幅 (%)", fontsize=6.8, color=MUTED)
    a.set_ylim(-110, 108)
    a.set_xlim(-0.15, max(x) + 0.2)
    a.set_title("前馈不是越早越好，是一条单峰曲线——峰在闸门来得及动作的地方",
                fontsize=7.0, color=INK, loc="left", pad=40)
    a.legend(fontsize=6.0, frameon=False, loc="lower right",
             bbox_to_anchor=(1.0, 1.0), ncol=2, handlelength=1.6,
             columnspacing=1.4)
    a.tick_params(labelsize=6.2, colors=MUTED)
    for s_ in ("top", "right"):
        a.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        a.spines[s_].set_color(HAIR)
    a.grid(color="#efeeea", lw=0.6)
    a.set_axisbelow(True)
    return save(fig, "r07_tau_peak.png")


def main() -> int:
    d = collect()
    d["scene_note"] = (
        f"传播时间 = 池长 ÷ 波速 = {d['geometry']['pool_len_m']:.0f} ÷ "
        f"{d['tau']['c_dynamic']:.3f} ≈ {d['tau']['tau_s']:.0f} s"
        f"（水从池首爬到池尾要 {d['tau']['tau_s'] / 60:.1f} min）")
    print(json.dumps(d, ensure_ascii=False, indent=2))
    for fn in (fig_scene, fig_ladder, fig_id_mech, fig_wavespeed,
               fig_ctrl_period, fig_real_inflow, fig_tau_peak):
        print("  出图:", os.path.basename(fn(d)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
