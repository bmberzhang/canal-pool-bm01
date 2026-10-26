# -*- coding: utf-8 -*-
"""图1 重画：受控反例三臂（统一 C_主案例 L=6000/ΔQ=45/T=300s/12h，判定量=池均水位偏差）。
三面板：
  (a) 0D 集总模型（反馈=池蓄量等效水位）→ 池均偏差峰值 72.7mm → A1
  (b) 1D 全池平均反馈 → 池均偏差峰值 75.05mm → A1
  (c) 1D 闸前单点反馈 → 池均偏差峰值 470.78mm → A3（晚期极限环）
数据源 bm01_fig1_rebuilt.json（真跑）。"""
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

d = json.load(open(os.path.join(RAW, "bm01_fig1_rebuilt.json"), encoding="utf-8"))
arms = d["arms"]
BAND = d["band_mm"]

BLUE = "#1a4fcc"; GREEN = "#2ca02c"; RED = "#d62728"; ORANGE = "#e68a00"

fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.4), sharex=True)
fig.suptitle("受控反例：同一阶跃、同一控制器、同一渠池，仅反馈取样尺度不同", fontsize=11, y=0.99)

def style(ax):
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.8)
    ax.grid(True, ls="-", lw=0.4, color="0.88")
    ax.tick_params(labelsize=8)

def draw_panel(ax, arm, color, label, verdict, note):
    style(ax)
    t = arms[arm]["t_h"]
    dev = arms[arm]["dev_mean_mm"]
    peak = arms[arm]["peak_mean_mm"]
    ax.axhspan(-BAND, BAND, color=GREEN, alpha=0.06, lw=0)
    ax.fill_between(t, 0, dev, color=color, alpha=0.12, lw=0)
    ax.plot(t, dev, color=color, lw=2.0, label=label)
    ax.axhline(BAND, color=GREEN, ls="-.", lw=1.2)
    ax.axhline(-BAND, color=GREEN, ls="-.", lw=1.2)
    ax.axvline(0.25, color="0.45", ls=":", lw=0.8)
    # 右上角固定格式：Peak=xx mm / 判定
    ax.text(0.985, 0.94, f"Peak = {peak:.1f} mm\n{verdict}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color=color, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=color, alpha=0.9))
    if note:
        ax.text(0.015, 0.06, note, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color="0.35")
    ax.set_ylabel("池均水位偏差 [mm]", fontsize=8)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

# (a) 0D
ax = axes[0]
draw_panel(ax, "0D", BLUE, "0D 集总模型（蓄量等效反馈）", "A1（可吸收）",
           "阶跃 +45 m³/s 于 0.25 h")
ax.set_ylim(-180, 190)
ax.set_title("(a) 0D 集总模型", fontsize=9)

# (b) 1D 池均
ax = axes[1]
draw_panel(ax, "mean", GREEN, "1D·全池平均反馈", "A1（可吸收）", None)
ax.set_ylim(-180, 190)
ax.set_title("(b) 1D·全池平均反馈（与 0D 几乎重合）", fontsize=9)

# (c) 1D 单点（判定量=池均，展示晚期极限环）
ax = axes[2]
draw_panel(ax, "point", RED, "1D·闸前单点反馈（判定量=池均）", "A3（不可吸收）",
           "晚期极限环：峰值在 t>2 h 才长成")
ax.set_ylim(-550, 550)
ax.set_xlabel("时间 [h]", fontsize=8)
ax.set_title("(c) 1D·闸前单点反馈", fontsize=9)

fig.tight_layout(rect=[0, 0, 1, 0.98])
os.makedirs(OUT, exist_ok=True)
outpath = os.path.join(OUT, "fig1_counterexample_v2.png")
fig.savefig(outpath, dpi=220, bbox_inches="tight")
print("SAVED", outpath)
