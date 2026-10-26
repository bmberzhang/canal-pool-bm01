# -*- coding: utf-8 -*-
"""图S1 重画：网格敏感性（12 h、W9 三个算例），修正旧版误用 4 h 数据（142.9/156.6/158.7mm）的问题。

三面板：
 (a) L=6000, ΔQ=45, T=300 s —— 翻转存在：池均 A1（74.7/75.0/74.8），单点 A3（813.8/470.8/1057.6）
 (b) L=2000, ΔQ=45, T=60 s  —— 翻转存在：池均 A1（102.5/102.3/108.2），单点 A3（164.8/175.0/184.1）
 (c) L=2000, ΔQ=45, T=300 s —— 翻转不存在：两档全 A3（≈4200 mm）
数据源 raw/bm01_w9_grid_check.json（12 h）。
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 9

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
OUT = os.path.join(HERE, "figures", "pub_bookcn")

d = json.load(open(os.path.join(RAW, "bm01_w9_grid_check.json"), encoding="utf-8"))
BAND = 150.0
GREEN = "#2ca02c"; RED = "#d62728"; BLUE = "#1a4fcc"; GREY = "#555555"

PANELS = [
    ("C_主案例", "(a) L=6000 m, ΔQ=45, T$_{ctrl}$=300 s", True),
    ("B_短池_快周期", "(b) L=2000 m, ΔQ=45, T$_{ctrl}$=60 s", True),
    ("A_短池_原周期", "(c) L=2000 m, ΔQ=45, T$_{ctrl}$=300 s", False),
]

fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.5))
for ax, (cid, title, flip) in zip(axes, PANELS):
    c = d["cases"][cid]
    rows = c["rows"]
    dx = [r["dx_m"] for r in rows]
    mean = [r["arms"]["mean"]["peak_mean_mm"] for r in rows]
    point = [r["arms"]["point"]["peak_mean_mm"] for r in rows]

    ax.axhspan(0, BAND, color=GREEN, alpha=0.07, lw=0)
    ax.axhline(BAND, color=GREEN, ls="-.", lw=1.3, label="判定条带 ±150 mm")
    ax.plot(dx, mean, "o-", color=GREEN, lw=2.0, ms=6, label="全池平均反馈（池均峰值）")
    ax.plot(dx, point, "s-", color=RED, lw=2.0, ms=6, label="闸前单点反馈（池均峰值）")
    for i, (x, y) in enumerate(zip(dx, mean)):
        ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=7.5, color=GREEN)
    for i, (x, y) in enumerate(zip(dx, point)):
        ax.annotate(f"{y:.1f}", (x, y), textcoords="offset points",
                    xytext=(0, -13 if y > 1000 else 9),
                    ha="center", fontsize=7.5, color=RED)

    ax.set_xscale("linear")
    ax.set_xticks(dx)
    ax.set_xticklabels([("%g" % v) for v in dx], fontsize=8)
    ax.minorticks_off()
    ax.set_xlim(min(dx) * 0.55, max(dx) * 1.18)
    ax.set_xlabel("网格间距 dx [m]", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.grid(True, ls="-", lw=0.4, color="0.88")
    lo = min(min(mean), min(point)); hi = max(max(mean), max(point))
    ax.set_ylim(0, hi * 1.34)
    for s in ax.spines.values():
        s.set_linewidth(0.8)

    tag = "翻转成立" if flip else "无翻转（两档全 A3）"
    ax.text(0.5, 0.035, tag, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8.5, fontweight="bold",
            color=(RED if flip else GREY),
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=(RED if flip else GREY), alpha=0.92))

axes[0].set_ylabel("池均水位偏差峰值 [mm]", fontsize=8)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8,
           frameon=False, bbox_to_anchor=(0.5, -0.02))

fig.tight_layout(rect=[0, 0.06, 1, 1])
os.makedirs(OUT, exist_ok=True)
p = os.path.join(OUT, "figS1_grid.png")
fig.savefig(p, dpi=220, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "figS1_grid.pdf"), bbox_inches="tight")
print("SAVED", p)
