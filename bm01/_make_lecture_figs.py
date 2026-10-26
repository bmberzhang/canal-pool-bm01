# -*- coding: utf-8 -*-
"""生成《BM-01 降阶对照·方法原理讲义》用的三张示意图（纯绘图，不改变任何案例结论）。

运行（需 pylibs-extra 提供 numpy/matplotlib）：
  PYTHONPATH=<venv>/pylibs-extra python _make_lecture_figs.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib import font_manager

OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures"
os.makedirs(OUT, exist_ok=True)

_avail = {f.name for f in font_manager.fontManager.ttflist}
for _c in ["Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC"]:
    if _c in _avail:
        plt.rcParams["font.sans-serif"] = [_c]
        break
plt.rcParams["axes.unicode_minus"] = False

INK, MUTED, HAIR = "#2c2c2a", "#5f5e5a", "#d3d1c7"
BLUE, AMBER, RED = "#185fa5", "#ba7517", "#a32d2d"
PANEL, WATER = "#f4f4f2", "#b5d4f4"


def fig_schematic() -> str:
    """0D 水缸 vs 1D 一排格子。"""
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 3.7))
    for a in axs:
        a.set_xlim(0, 1.3)
        a.set_ylim(-0.12, 1.18)
        a.axis("off")

    # ---- 左：0D ----
    a = axs[0]
    a.add_patch(Rectangle((0.02, 0.10), 1.00, 0.80,
                          facecolor=PANEL, edgecolor="#b4b2a9", lw=0.9, zorder=1))
    a.add_patch(Rectangle((0.025, 0.105), 0.99, 0.455,
                          facecolor=WATER, edgecolor="none", zorder=1))
    a.plot([0.02, 1.02], [0.56, 0.56], color=BLUE, lw=1.8, zorder=3)
    a.plot([0.02, 1.02], [0.56, 0.56], color=MUTED, lw=0.9,
           ls=(0, (4, 3)), zorder=5)
    a.add_patch(Rectangle((1.02, 0.35), 0.04, 0.26,
                          facecolor="#5f5e5a", edgecolor="none", zorder=3))
    a.text(0.52, 0.31, "水面处处水平，正好压在目标水位上", ha="center",
           va="center", fontsize=11, color=MUTED)
    a.text(0.52, 1.02, "集总 0D（一个大水缸）", ha="center", va="bottom",
           fontsize=12, color=INK)
    a.text(0.52, -0.02, "闸门量到的是全池水位", ha="center", va="top",
           fontsize=11, color=MUTED)

    # ---- 右：1D ----
    a = axs[1]
    a.add_patch(Rectangle((0.02, 0.10), 1.00, 0.80,
                          facecolor=PANEL, edgecolor="#b4b2a9", lw=0.9, zorder=1))

    def _surf(xq):
        return (0.540
                - 0.045 * np.exp(-((xq - 0.25) / 0.32) ** 2)
                + 0.190 * np.exp(-((xq - 1.06) / 0.52) ** 2))

    ncell = 12
    edg = np.linspace(0.02, 1.02, ncell + 1)
    top = [float(_surf(0.5 * (edg[k] + edg[k + 1]))) for k in range(ncell)]
    for k in range(ncell):
        a.add_patch(Rectangle((edg[k] + 0.003, 0.105),
                              edg[k + 1] - edg[k] - 0.006, top[k] - 0.105,
                              facecolor=WATER, edgecolor="none", zorder=2))
    sx, sy = [], []
    for k in range(ncell):
        sx += [edg[k] + 0.003, edg[k + 1] - 0.003]
        sy += [top[k], top[k]]
    a.plot(sx, sy, color=BLUE, lw=1.8, zorder=4)
    a.plot([0.02, 1.02], [0.56, 0.56], color=MUTED, lw=0.9,
           ls=(0, (4, 3)), zorder=4)
    a.text(0.035, 0.585, "目标水位", ha="left", va="bottom",
           fontsize=10, color=MUTED, zorder=5)
    a.add_patch(Rectangle((1.02, 0.35), 0.04, 0.26,
                          facecolor="#5f5e5a", edgecolor="none", zorder=3))
    a.text(0.52, 1.02, "一维 1D（切成 N 个小格子）", ha="center", va="bottom",
           fontsize=12, color=INK)
    a.text(0.52, -0.02, "闸门只量到下游端最后一格", ha="center", va="top",
           fontsize=11, color=MUTED)
    a.text(0.52, -0.115, "（示意，非按比例）", ha="center", va="top",
           fontsize=9, color="#999999")

    fig.tight_layout()
    p = os.path.join(OUT, "bm01_lecture_schematic_0d_1d.png")
    fig.savefig(p, dpi=170, facecolor="white")
    plt.close(fig)
    return p


def fig_cfl() -> str:
    """CFL：一个时间步内波最多跨一格。"""
    fig, a = plt.subplots(figsize=(11.0, 3.1))
    a.set_xlim(0, 5.6)
    a.set_ylim(-0.05, 2.85)
    a.axis("off")
    x0, W, H = 0.2, 1.0, 0.62
    mid = x0 + 5 * W / 2.0

    # 上排：稳定
    y1 = 1.86
    for k in range(5):
        a.add_patch(Rectangle((x0 + k * W, y1), W, H,
                              facecolor=PANEL, edgecolor=HAIR, lw=0.8))
    cy = y1 + H / 2.0
    a.annotate("", xy=(1.55, cy), xytext=(0.50, cy),
               arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.4,
                               shrinkA=0, shrinkB=0))
    a.plot([0.50], [cy], "o", color=BLUE, ms=6)
    a.plot([1.55], [cy], "o", color=BLUE, ms=6, alpha=0.35)
    a.text(1.02, cy + 0.10, "Δt", ha="center", fontsize=11, color=BLUE)
    a.text(mid, y1 + H + 0.14, "稳定：Δt 内波最多跨一格", ha="center",
           fontsize=12, color=INK)
    a.text(mid, y1 - 0.14, "每一格的信息都被用上", ha="center", va="top",
           fontsize=11, color=MUTED)

    # 下排：不稳定
    y2 = 0.40
    a.add_patch(Rectangle((x0 + W, y2), W, H,
                          facecolor="#e24b4a", alpha=0.16, edgecolor="none"))
    for k in range(5):
        a.add_patch(Rectangle((x0 + k * W, y2), W, H,
                              facecolor="none", edgecolor=HAIR, lw=0.8))
    cy2 = y2 + H / 2.0
    a.annotate("", xy=(2.60, cy2), xytext=(0.50, cy2),
               arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.4,
                               shrinkA=0, shrinkB=0))
    a.plot([0.50], [cy2], "o", color=BLUE, ms=6)
    a.plot([2.60], [cy2], "o", color=BLUE, ms=6, alpha=0.35)
    a.text(1.55, cy2 + 0.10, "Δt", ha="center", fontsize=11, color=BLUE)
    a.text(mid, y2 + H + 0.14, "不稳定：Δt 内跨过了整格", ha="center",
           fontsize=12, color=INK)
    a.text(x0 + 1.5 * W, y2 - 0.14, "中间这一格的信息从头到尾没被用上",
           ha="center", va="top", fontsize=11, color=RED)

    fig.tight_layout()
    p = os.path.join(OUT, "bm01_lecture_cfl.png")
    fig.savefig(p, dpi=170, facecolor="white")
    plt.close(fig)
    return p


def fig_attribution() -> str:
    """归因柱状图：只换反馈量，峰值差一倍。"""
    labels = ["集总 0D\n（记录值）", "1D 池均\n120 格", "1D 池均\n240 格",
              "1D 闸前\n120 格", "1D 闸前\n240 格", "1D 闸前\n480 格",
              "外推\n（网格无关）"]
    vals = [72.7, 74.7, 75.0, 142.9, 156.6, 158.7, 160.8]
    colors = [BLUE, BLUE, BLUE, AMBER, AMBER, AMBER, RED]

    fig, a = plt.subplots(figsize=(9.8, 4.1))
    bars = a.bar(range(len(vals)), vals, color=colors, width=0.62, zorder=3)
    a.axhline(150, color=RED, ls="--", lw=1.2, zorder=4)
    a.text(-0.42, 152.5, "正常条带上限 150 mm",
           ha="left", va="bottom", fontsize=10, color=RED)
    for b, v in zip(bars, vals):
        a.text(b.get_x() + b.get_width() / 2.0, v + 3, f"{v:.1f}",
               ha="center", fontsize=10, color=INK, zorder=5)

    a.set_xticks(range(len(vals)))
    a.set_xticklabels(labels, fontsize=9.5, color=MUTED)
    a.set_ylabel("池水位峰值偏差 (mm)", fontsize=10, color=MUTED)
    a.set_ylim(0, 195)
    a.tick_params(axis="y", labelsize=10, colors=MUTED)
    a.tick_params(axis="x", length=0)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_color(HAIR)
    a.grid(axis="y", color="#e8e6e0", lw=0.7, zorder=0)
    a.set_axisbelow(True)
    a.legend(handles=[Patch(color=BLUE, label="反馈取池均水位"),
                      Patch(color=AMBER, label="反馈取闸前局部水位"),
                      Patch(color=RED, label="外推 / 越出条带")],
             fontsize=9, frameon=False, ncol=3,
             loc="upper center", bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout()
    p = os.path.join(OUT, "bm01_lecture_attribution.png")
    fig.savefig(p, dpi=170, facecolor="white")
    plt.close(fig)
    return p


def fig_transient() -> str:
    """讲义版时间过程图：读 raw/bm01_sve1d_compare.csv，不重算物理。

    与 figures/bm01_lumped_vs_sve1d.png 同数据、同信息，改为「每页 15 cm 宽」
    的竖排三联版式——原始那张是按 16.4 英寸画布出的，缩到 15 cm 后字号只剩
    约 4 pt，无法阅读。
    """
    import csv
    import json

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "raw", "bm01_sve1d_compare.csv"),
              encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    with open(os.path.join(here, "raw", "bm01_sve1d_compare.json"),
              encoding="utf-8") as fh:
        meta = json.load(fh)

    t = [float(r["t_s"]) / 60.0 for r in rows]
    col = lambda k: [float(r[k]) * 1000 for r in rows]      # noqa: E731
    lump = col("lumped_dev_m")
    f_mean = col("g240_gate_dev_mean_m")
    g_gate = col("g120_gate_dev_gate_m")
    m_mean = col("g120_mean_dev_mean_m")

    lp = meta["lumped_recorded"]["peak_mm"]
    mp = meta["attribution"]["sve_mean_feedback_mm"]["120"]
    fm = meta["sve"]["g240_gate"]["peak_mean_mm"]
    gp = meta["sve"]["g120_gate"]["peak_gate_mm"]
    ext = meta["attribution"]["extrapolated_mm"]
    gfin = meta["lumped_recorded"]["gate_final"]
    BAND, STEP = 150.0, 15.0

    fig, ax = plt.subplots(3, 1, figsize=(5.9055, 6.3))
    fig.suptitle("BM-01 降阶损失对照：集总 0D vs 一维圣维南（同参数·同 PI·量级示例）",
                 fontsize=9.5, fontweight="bold")

    a = ax[0]
    a.axhspan(-BAND, BAND, color="#e8f3ff", zorder=0)
    for y in (BAND, -BAND):
        a.axhline(y, color="#1a7f37", lw=0.9)
    a.axvline(STEP, color="#888888", lw=0.8, ls="--")
    a.plot(t, lump, color="#0066cc", lw=1.4, label="集总 0D（主案例）峰值 73 mm → A1")
    a.plot(t, m_mean, color="#8e44ad", lw=1.0, ls="--",
           label=f"1D 但反馈取池均水位 {mp:.0f} mm（≈0D）")
    a.plot(t, f_mean, color="#1a7f37", lw=1.2,
           label=f"1D 池均水位（240 格）{fm:.0f} mm")
    a.annotate("正常条带 ±150 mm", xy=(6, BAND), xytext=(6, BAND + 14),
               fontsize=6.5, color="#1a7f37")
    a.annotate(f"外推 ~{ext:.0f} mm\n越过条带 → 判定翻成 A3", xy=(28, 158),
               xytext=(62, 108), fontsize=6.5,
               arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))
    a.set_title("(a) 池均水位峰值：0D 73 mm ／ 1D 池均 75 mm ／ 1D 闸前 157 mm",
                fontsize=8)
    a.set_xlabel("时间 (min)", fontsize=7.5)
    a.set_ylabel("水位偏差 (mm)", fontsize=7.5)
    a.set_xlim(0, max(t)); a.set_ylim(-60, 205)
    a.legend(fontsize=6, loc="upper right", framealpha=0.9)
    a.grid(alpha=.25); a.tick_params(labelsize=6.5)

    b = ax[1]
    for y in (BAND, -BAND):
        b.axhline(y, color="#1a7f37", lw=0.9)
    b.axvline(STEP, color="#888888", lw=0.8, ls="--")
    b.plot(t, f_mean, color="#1a7f37", lw=1.2, label="池均水位（池子的平均状态）")
    b.plot(t, g_gate, color="#e67e22", lw=1.0,
           label=f"闸前局部水位（闸门／传感器看到的）峰值 {gp:.0f} mm")
    b.set_title("(b) 闸门量到的是局部水位，不是池均水位", fontsize=8)
    b.set_xlabel("时间 (min)", fontsize=7.5)
    b.set_ylabel("水位偏差 (mm)", fontsize=7.5)
    b.set_xlim(0, max(t))
    b.legend(fontsize=6, loc="upper left", framealpha=0.9)
    b.grid(alpha=.25); b.tick_params(labelsize=6.5)

    c = ax[2]
    c.plot(t, [float(r["lumped_gate"]) for r in rows], color="#0066cc", lw=1.3,
           label="集总 0D 闸门开度")
    c.plot(t, [float(r["g120_gate_gate"]) for r in rows], color="#1a7f37", lw=1.3,
           label="1D SVE 闸门开度")
    c.axhline(gfin, color="#888888", lw=0.9, ls="--",
              label=f"终态开度 {gfin:.3f}")
    c.set_title("(c) 闸门动作（同一 PI·采样周期 300 s·速率限）", fontsize=8)
    c.set_xlabel("时间 (min)", fontsize=7.5)
    c.set_ylabel("开度 (0~1)", fontsize=7.5)
    c.set_xlim(0, max(t))
    c.legend(fontsize=6, loc="upper right", framealpha=0.9)
    c.grid(alpha=.25); c.tick_params(labelsize=6.5)

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    p = os.path.join(OUT, "bm01_lecture_transient.png")
    fig.savefig(p, dpi=170, facecolor="white")
    plt.close(fig)
    return p


if __name__ == "__main__":
    for f in (fig_schematic, fig_cfl, fig_attribution, fig_transient):
        print("写出:", f())
