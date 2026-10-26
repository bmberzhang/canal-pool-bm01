# -*- coding: utf-8 -*-
"""
BM-01 论文图件 v3 —— 对齐《水系统先进控制与网络协同》书图风格：
- 图顶总标题 + (a)(b)(c) 面板标题
- 粗彩线（蓝 #1f4fd8 / 红 #d62728 / 绿 #2ca02c）、关键曲线下浅色填充
- 阈值线点划线 + 图例带框并写物理含义（如 "正常调节条带 (±150 mm)"）
- 箭头标注关键特征（洪峰/极限环/峰值数值）
- 四边框封闭坐标轴、浅色双向网格、y 轴标签带单位 [mm]
输出 figures/pub_bookcn/
"""
import json, math, os
import numpy as np
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

RAW = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/raw"
OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures/pub_bookcn"
os.makedirs(OUT, exist_ok=True)

BLUE = "#1a4fcc"; RED = "#d62728"; GREEN = "#2ca02c"; ORANGE = "#e68a00"; PURPLE = "#7b52ab"

def style_ax(ax):
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.8)
    ax.grid(True, ls="-", lw=0.4, color="0.88")
    ax.tick_params(labelsize=8)

def load(p):
    with open(os.path.join(RAW, p), encoding="utf-8") as f:
        return json.load(f)

# ============================================================
# 图 1：受控反例 —— 三面板 + 阶跃标注 + 极限环标注 + 峰值 inset
# ============================================================
def fig1():
    main = load("bm01_main.json")
    sv = load("bm01_sve1d_compare.json")
    t0 = np.array(main["series_downsampled"]["t_s"]) / 3600.0
    d0 = np.array(main["series_downsampled"]["pool1_dev_m"]) * 1000.0
    t1 = np.array(sv["series"]["t"]) / 3600.0
    d_mean = np.array(sv["series"]["dev_mean"]) * 1000.0
    d_gate = np.array(sv["series"]["dev_gate"]) * 1000.0
    band = 150.0

    fig, axes = plt.subplots(3, 1, figsize=(6.8, 6.2), sharex=True)
    fig.suptitle("受控反例：反馈取样位置决定可吸收判定", fontsize=11, y=0.98)

    # (a) 0D
    ax = axes[0]; style_ax(ax)
    ax.axhspan(-band, band, color=GREEN, alpha=0.06, lw=0)
    ax.fill_between(t0, 0, d0, color=BLUE, alpha=0.12, lw=0)
    ax.plot(t0, d0, color=BLUE, lw=2.0, label="0D 集总模型水位偏差")
    ax.axhline(band, color=GREEN, ls="-.", lw=1.3)
    ax.axhline(-band, color=GREEN, ls="-.", lw=1.3)
    ax.axvline(0.25, color="0.45", ls=":", lw=0.9)
    ax.annotate("阶跃 +45 m³/s（0.25 h）", xy=(0.25, 40), xytext=(1.1, 95),
                fontsize=7.5, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
    ax.annotate(f"峰值 {d0.max():.1f} mm", xy=(t0[np.argmax(d0)], d0.max()),
                xytext=(2.3, 110), fontsize=7.5, color=BLUE,
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
    ax.set_ylabel("水位偏差 [mm]", fontsize=8)
    ax.set_ylim(-180, 190)
    ax.set_title("(a) 0D 集总模型（判定 A1）", fontsize=9)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # (b) 1D 池均
    ax = axes[1]; style_ax(ax)
    ax.axhspan(-band, band, color=GREEN, alpha=0.06, lw=0)
    ax.plot(t1, d_mean, color=GREEN, lw=2.0, label="1D·全池平均反馈")
    ax.axhline(band, color=GREEN, ls="-.", lw=1.3)
    ax.axhline(-band, color=GREEN, ls="-.", lw=1.3)
    ax.annotate(f"峰值 {d_mean.max():.1f} mm", xy=(t1[np.argmax(d_mean)], d_mean.max()),
                xytext=(2.3, 120), fontsize=7.5, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.9))
    ax.set_ylabel("水位偏差 [mm]", fontsize=8)
    ax.set_ylim(-180, 190)
    ax.set_title("(b) 1D·全池平均反馈（判定 A1，与 0D 几乎重合）", fontsize=9)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    # (c) 1D 闸前单点
    ax = axes[2]; style_ax(ax)
    ax.axhspan(-band, band, color=RED, alpha=0.05, lw=0)
    ax.plot(t1, d_gate, color=RED, lw=1.6, label="1D·闸前单点反馈")
    ax.axhline(band, color=GREEN, ls="-.", lw=1.3, label="正常调节条带 (±150 mm)")
    ax.axhline(-band, color=GREEN, ls="-.", lw=1.3)
    ax.axvspan(1.0, 4.5, color=ORANGE, alpha=0.07, lw=0)
    ax.annotate("持续极限环（不收敛）", xy=(3.0, 430), xytext=(0.7, 480),
                fontsize=7.5, color=ORANGE,
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
    ax.annotate(f"峰值 {d_gate.max():.0f} mm", xy=(t1[np.argmax(d_gate)], d_gate.max()),
                xytext=(2.6, 350), fontsize=7.5, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    ax.set_ylabel("水位偏差 [mm]", fontsize=8)
    ax.set_xlabel("时间 [h]", fontsize=8)
    ax.set_ylim(-380, 580)
    ax.set_title("(c) 1D·闸前单点反馈（判定 A3）", fontsize=9)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(OUT, "fig1_counterexample.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig1_counterexample.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig1] done")

# ============================================================
# 图 2：SWMM 双求解器独立复现（分组棒图 + 判定码 + min~max 须）
# ============================================================
def fig2():
    ext = load("bm01_ext_compare.json")
    from collections import defaultdict
    g = defaultdict(list)
    for r in ext["rows"]:
        g[(r["code"], r["ctrl"])].append(r["peak_pool_mm"])
    order = [("SVE-1D", "mean"), ("EPA SWMM", "pool"),
             ("SVE-1D", "gate"), ("EPA SWMM", "gate")]
    labels = ["本仓库 SVE-1D\n全池平均", "EPA SWMM\n全池平均",
              "本仓库 SVE-1D\n闸前单点", "EPA SWMM\n闸前单点"]
    colors = [GREEN, GREEN, RED, RED]
    codes = ["A1", "A1", "A3", "A3"]
    med = [float(np.median(g[k])) for k in order]
    lo = [float(min(g[k])) for k in order]
    hi = [float(max(g[k])) for k in order]

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    style_ax(ax)
    fig.suptitle("双求解器独立复现：判定翻转跨代码成立", fontsize=10.5, y=0.98)
    x = np.arange(len(order))
    for i in range(len(order)):
        ax.bar(x[i], med[i], color=colors[i], width=0.55, alpha=0.55,
               edgecolor=colors[i], linewidth=1.2)
        ax.vlines(x[i], lo[i], hi[i], color="0.25", lw=1.3)
        ax.text(x[i], hi[i]+14, f"{med[i]:.0f}", ha="center", fontsize=7.5)
        ax.text(x[i], -48, codes[i], ha="center", fontsize=9.5,
                fontweight="bold", color=colors[i])
    ax.axvline(1.5, color="0.5", lw=1.0, ls=":")
    ax.axhline(150, color=ORANGE, ls="-.", lw=1.4, label="正常调节条带 (±150 mm)")
    ax.text(3.4, 160, "判定阈值", fontsize=7, color=ORANGE, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("全池平均水位峰值 [mm]", fontsize=8)
    ax.set_ylim(0, 540)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(OUT, "fig2_swmm.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig2_swmm.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig2] done")

# ============================================================
# 图 3：量测可实现性（4 估计量 + 理想参考 + 阈值 + 0.0mm 标注）
# ============================================================
def fig3():
    obs = load("bm01_obs_estimator.json")
    ref = obs["reference_arms"]["mean"]["peak_mean_mm"]
    est = obs["estimated_arms"]
    arms = ["midpoint", "mp3", "mp5", "volbal"]
    labels = ["中点", "3点平均", "5点平均", "容积平衡"]
    vals = [est[a]["peak_mean_mm"] for a in arms]
    colors = [PURPLE, PURPLE, PURPLE, GREEN]

    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    style_ax(ax)
    fig.suptitle("量测可实现性：容积平衡估计量与理想池均完全重合", fontsize=10.5, y=0.98)
    x = np.arange(len(arms))
    for i, v in enumerate(vals):
        ax.bar(x[i], v, color=colors[i], width=0.5, alpha=0.6,
               edgecolor=colors[i], linewidth=1.2)
        ax.text(i, v+3, f"{v:.1f}", ha="center", fontsize=7.5)
        ax.text(i, -13, "A1", ha="center", fontsize=9, fontweight="bold", color=GREEN)
    ax.plot([3], [vals[3]+12], marker="*", ms=11, color=GREEN)
    ax.annotate("偏差 0.0 mm\n（仅需上下游流量计）", xy=(3, vals[3]+14), xytext=(1.9, 135),
                fontsize=7.5, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.9))
    ax.axhline(ref, color=BLUE, ls="--", lw=1.4, label="理想全池平均 (75.0 mm)")
    ax.axhline(150, color=ORANGE, ls="-.", lw=1.3, label="正常调节条带 (±150 mm)")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("全池平均水位峰值 [mm]", fontsize=8)
    ax.set_ylim(0, 175)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(OUT, "fig3_measurability.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig3_measurability.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig3] done")

# ============================================================
# 图 S1（补充）：网格敏感性
# ============================================================
def figS1():
    sv = load("bm01_sve1d_compare.json")
    mean_v = [sv["sve"]["g120_mean"]["peak_mean_mm"], sv["sve"]["g240_mean"]["peak_mean_mm"]]
    gate_v = [sv["sve"]["g120_gate"]["peak_mean_mm"],
              sv["sve"]["g240_gate"]["peak_mean_mm"],
              sv["sve"]["g480_gate"]["peak_mean_mm"]]
    fig, ax = plt.subplots(figsize=(4.2, 3.1))
    style_ax(ax)
    fig.suptitle("网格敏感性：池均档收敛，单点档未收敛", fontsize=10.5, y=0.98)
    ax.plot([50.0, 25.0], mean_v, "-o", color=GREEN, lw=2.0, ms=6, label="全池平均反馈")
    ax.plot([50.0, 25.0, 12.5], gate_v, "--s", color=RED, lw=2.0, ms=6, label="闸前单点反馈")
    ax.axhline(150, color=ORANGE, ls="-.", lw=1.3, label="正常调节条带 (±150 mm)")
    ax.axhline(72.74, color=BLUE, ls="--", lw=1.2, label="0D 集总 (72.7 mm)")
    ax.annotate("随细化未收敛（网格敏感）", xy=(12.5, gate_v[-1]), xytext=(17, 195),
                fontsize=7.5, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    ax.set_xlabel("网格尺寸 dx [m]（向右细化）", fontsize=8)
    ax.set_ylabel("全池平均水位峰值 [mm]", fontsize=8)
    ax.set_xlim(12, 54); ax.invert_xaxis()
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(OUT, "figS1_grid.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS1_grid.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS1] done")

# ============================================================
# 图 S2（补充）：空间低通机制
# ============================================================
def figS2():
    x = np.linspace(0.02, 2.0, 400)
    theory = np.sin(np.pi*x)**2 / (np.pi*x)
    w5 = load("bm01_w5_kc_identity.json")
    xs = [r["L_over_lambda"] for r in w5["rows"] if r["L_m"] == 6000]
    num = [r["numeric"] for r in w5["rows"] if r["L_m"] == 6000]
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    style_ax(ax)
    fig.suptitle("全池平均 = 空间低通滤波器", fontsize=10.5, y=0.98)
    ax.plot(x, theory, color=BLUE, lw=2.0, label="解析：sin²(πx)/(πx)")
    ax.plot(xs, num, "o", color=RED, ms=6, mfc="none", ls="none", label="1D SVE 数值 (dx=25 m)")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axvline(0.5, color=ORANGE, ls=":", lw=1.1)
    ax.axhline(2/math.pi, color=ORANGE, ls=":", lw=1.1)
    ax.plot([0.5], [2/math.pi], "s", color=ORANGE, ms=6)
    ax.annotate("半功率点 x=0.5\n响应 = 2/π ≈ 0.637", xy=(0.5, 2/math.pi), xytext=(0.85, 0.82),
                fontsize=7.5, color=ORANGE,
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
    ax.set_xlabel("x = L/λ（池长 / 扰动波长）", fontsize=8)
    ax.set_ylabel("归一化响应 [–]", fontsize=8)
    ax.set_ylim(-0.15, 1.08)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(OUT, "figS2_lowpass.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS2_lowpass.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS2] done")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); figS1(); figS2()
    print("BOOKCN-STYLE FIGURES DONE ->", OUT)
