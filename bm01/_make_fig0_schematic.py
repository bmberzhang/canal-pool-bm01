# -*- coding: utf-8 -*-
"""图 0 v2：修复文字重叠，重新布局"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

for fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures/pub_bookcn"
BLUE = "#1a4fcc"; WATER = "#7fb3e8"; GREEN = "#2ca02c"; RED = "#d62728"
DARK = "#222222"; GRAY = "#666666"; ORANGE = "#e68a00"

fig, ax = plt.subplots(figsize=(10.0, 5.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 60)
ax.axis("off")

def box(x, y, w, h, fc, ec=DARK, lw=1.2, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0,rounding_size={r}",
                 fc=fc, ec=ec, lw=lw, zorder=3))

def arrow(x1, y1, x2, y2, color=DARK, lw=1.6, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=14, color=color, lw=lw, linestyle=ls, zorder=4))

# ---------- 上游水源 ----------
wx, wy, ww, wh = 2, 30, 13, 20
box(wx, wy, ww, wh, "#eef5fd")
ax.add_patch(Rectangle((wx+0.8, wy+0.8), ww-1.6, 7, fc=WATER, ec="none", zorder=4))
ax.text(wx+ww/2, wy+wh-3.5, "上游水源", ha="center", fontsize=10, fontweight="bold", color=BLUE)
ax.text(wx+ww/2, wy+11.5, "入库流量 $Q_{up}$\n阶跃 +45 m³/s\n($q_0$ = 180 m³/s)",
        ha="center", fontsize=7.5, color=DARK)

# ---------- 明渠渠池 ----------
px, py, pw, ph = 21, 28, 45, 22
box(px, py, pw, ph, "#f5fafd")
ax.text(px+pw/2, py+ph-3.2, "明渠渠池（BM-01）", ha="center",
        fontsize=10.5, fontweight="bold", color=BLUE)
ax.text(px+pw/2, py+ph-7.0, "池长 L = 6000 m　池宽 B = 50 m",
        ha="center", fontsize=7.5, color=DARK)
# 水体 + 波浪线
ax.add_patch(Rectangle((px+1, py+1.2), pw-2, 6.5, fc=WATER, ec="none", alpha=0.85, zorder=2))
xs = [px+1+i*1.0 for i in range(int((pw-2)))] + [px+pw-1]
ys = [py+7.7 + (0.45 if i % 2 == 0 else -0.45) for i in range(len(xs))]
ax.plot(xs, ys, color=BLUE, lw=1.0, zorder=3)
ax.text(px+pw/2, py+3.4, "水位 $h(x,t)$", ha="center", fontsize=8, color="#083a75", zorder=4)

# 进流箭头
arrow(wx+ww, wy+8, px-0.4, py+6, color=BLUE, lw=2.2)
ax.text((wx+ww+px)/2-1.5, wy+11, "$Q_{up}$", fontsize=8, color=BLUE)

# ---------- 节制闸 ----------
gx, gy, gw, gh = 72, 28, 9, 22
box(gx, gy, gw, gh, "#fdf3e8")
ax.text(gx+gw/2, gy+gh-3.2, "节制闸", ha="center", fontsize=9.5, fontweight="bold", color=ORANGE)
ax.plot([gx+gw/2, gx+gw/2], [gy+gh-5.5, gy+gh-14], color=ORANGE, lw=3.5, zorder=4)
ax.text(gx+gw/2, gy+gh-16.8, "闸门开度 a", ha="center", fontsize=7, color=ORANGE)
arrow(gx+gw, gy+6, 97, gy+6, color=DARK, lw=2.0)
ax.text(88, gy+8.6, "$Q_{down}$\n至下游渠段", ha="center", fontsize=7.5, color=DARK)

# ---------- 反馈取样（标注放在渠池框下方，避免重叠） ----------
# ① 闸前单点（红）
m1x = px+pw+2.2
ax.plot([m1x], [py+7.7], marker="v", ms=9, color=RED, zorder=5)
ax.text(m1x, py+11, "①", ha="center", fontsize=9, color=RED, fontweight="bold")
# ② 分布测点（绿）
for fx in [px+6, px+14.5, px+23, px+31.5, px+39]:
    ax.plot([fx], [py+7.7], marker="o", ms=4.5, mfc="white", mec=GREEN, mew=1.4, zorder=5)
ax.text(px+23, py+10.2, "② 分布测点 → 全池平均", ha="center",
        fontsize=7.5, color=GREEN, fontweight="bold")
ax.text(m1x-2.2, py+10.2, "① 闸前单点", ha="right", fontsize=7.5,
        color=RED, fontweight="bold")

# ---------- 控制器 ----------
cx, cy, cw, ch = 33, 7, 34, 13
box(cx, cy, cw, ch, "#eef7ee", ec=GREEN, lw=1.5)
ax.text(cx+cw/2, cy+ch-3, "PID 控制器", ha="center", fontsize=10, fontweight="bold", color=GREEN)
ax.text(cx+cw/2, cy+ch-6.4, "采样周期 $T_{ctrl}$ = 300 s｜限幅 + 积分抗饱和",
        ha="center", fontsize=7.2, color=DARK)
ax.text(cx+cw/2, cy+2.8, "反馈量 = ？（唯一变量）", ha="center",
        fontsize=8.5, color=RED, fontweight="bold")

# 反馈路径
arrow(px+23, py-1.2, cx+9, cy+ch+0.3, color=GREEN, lw=1.6, ls="-.")
arrow(m1x, py-1.2, cx+cw-6, cy+ch+0.3, color=RED, lw=1.6, ls="--")
ax.text(30, 24.0, "路径② 全池平均 → A1", fontsize=7.5, color=GREEN, ha="center")
ax.text(63, 24.0, "路径① 闸前单点 → A3", fontsize=7.5, color=RED, ha="center")
# 执行
arrow(cx+cw, cy+ch/2+1, gx+gw/2+0.6, gy-0.6, color=ORANGE, lw=1.8)
ax.text(74.5, 21.5, "执行", fontsize=7.5, color=ORANGE, ha="center")

# ---------- 底部结论条 ----------
box(12, 0.5, 76, 4.6, "#fffbe8", ec=GRAY, lw=0.8, r=0.8)
ax.text(50, 2.8, "同一阶跃、同一控制器、同一渠池 —— 唯一变量是反馈水位取在哪里：判定在 A1 / A3 之间翻转",
        ha="center", fontsize=8.5, color=DARK)

fig.savefig(os.path.join(OUT, "fig0_schematic.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig0_schematic.png"), dpi=220, bbox_inches="tight")
print("SCHEMATIC V2 DONE")
