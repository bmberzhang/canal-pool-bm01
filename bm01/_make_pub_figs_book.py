# -*- coding: utf-8 -*-
"""
BM-01 论文图件 v3 —— Litrico & Fromion《Modeling and Control of Hydrosystems》
经典控制书风格：全黑线条、线型区分、四边框、细点网格、灰阶。
输出到 figures/pub_book/（彩色版保留在 figures/pub/ 供对比）。
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
OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures/pub_book"
os.makedirs(OUT, exist_ok=True)

BLACK = "0.0"; DARK = "0.25"; MID = "0.45"

def style_ax(ax):
    """控制书风格：四边框封闭 + 双向细点网格"""
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.8)
    ax.grid(True, which="major", ls=":", lw=0.5, color="0.75", alpha=0.9)
    ax.tick_params(labelsize=8, direction="in")

def load(p):
    with open(os.path.join(RAW, p), encoding="utf-8") as f:
        return json.load(f)

def tag(ax, code, pos=(0.86, 0.88)):
    ax.text(*pos, code, transform=ax.transAxes, fontsize=9, fontweight="bold",
            color=BLACK, ha="center", va="center")

# ============================================================
# 图 1（正文）：受控反例三面板（黑线线型区分 + inset）
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

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
    ax = axes[0]; style_ax(ax)
    ax.axhspan(-band, band, color="0.90", lw=0)
    ax.plot(t0, d0, color=BLACK, lw=1.6, ls="-")
    ax.set_ylim(-160, 170)
    ax.set_title("(a) 0D 集总", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8); ax.set_ylabel("水位偏差 (mm)", fontsize=8)
    tag(ax, "A1")
    ax.text(0.25, 165, "+45 m³/s 阶跃", fontsize=7, color=DARK, ha="center")

    ax = axes[1]; style_ax(ax)
    ax.axhspan(-band, band, color="0.90", lw=0)
    ax.plot(t1, d_mean, color=BLACK, lw=1.6, ls="-")
    ax.set_ylim(-160, 170)
    ax.set_title("(b) 1D·全池平均反馈", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8)
    tag(ax, "A1")

    ax = axes[2]; style_ax(ax)
    ax.axhspan(-band, band, color="0.90", lw=0)
    ax.plot(t1, d_gate, color=BLACK, lw=1.6, ls="-")
    ax.set_ylim(-350, 550)
    ax.set_title("(c) 1D·闸前单点反馈", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8)
    tag(ax, "A3", (0.86, 0.90))
    ax.text(2.7, 520, "持续极限环", fontsize=7, color=DARK, ha="center")

    # inset：0D vs 池均 峰值放大（线型区分：0D 实线 + 空心圆，池均 虚线 + 实心点）
    axins = axes[2].inset_axes([0.08, 0.10, 0.42, 0.40])
    axins.plot(t0, d0, color=BLACK, lw=1.3, ls="-", marker="o", ms=2.5, mfc="none", markevery=8)
    axins.plot(t1, d_mean, color=DARK, lw=1.3, ls="--")
    axins.axhline(band, color=MID, ls=":", lw=0.8)
    axins.set_xlim(0, 1.6); axins.set_ylim(0, 90)
    axins.tick_params(labelsize=6, direction="in")
    axins.set_title("峰值放大（0D ≈ 池均）", fontsize=6.5)
    axins.legend(["0D", "1D 池均"], fontsize=6, frameon=False, loc="lower right")
    for s in axins.spines.values():
        s.set_linewidth(0.6)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_counterexample.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig1_counterexample.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig1] done")

# ============================================================
# 图 2（正文）：SWMM 双求解器复现（灰阶棒图 + 纹理 + 判定码）
# ============================================================
def fig2():
    ext = load("bm01_ext_compare.json")
    from collections import defaultdict
    g = defaultdict(list)
    for r in ext["rows"]:
        g[(r["code"], r["ctrl"])].append(r["peak_pool_mm"])
    order = [("SVE-1D", "mean"), ("EPA SWMM", "pool"),
             ("SVE-1D", "gate"), ("EPA SWMM", "gate")]
    labels = ["本仓库\n全池平均", "SWMM\n全池平均", "本仓库\n闸前单点", "SWMM\n闸前单点"]
    codes = ["A1", "A1", "A3", "A3"]
    med = [float(np.median(g[k])) for k in order]
    lo = [float(min(g[k])) for k in order]
    hi = [float(max(g[k])) for k in order]

    fig, ax = plt.subplots(figsize=(4.0, 2.9))
    style_ax(ax)
    x = np.arange(len(order))
    for i in range(len(order)):
        fc = "0.85" if codes[i] == "A1" else "0.55"   # A1 浅灰 / A3 深灰
        hatch = "//" if codes[i] == "A1" else "xx"
        ax.bar(x[i], med[i], color=fc, width=0.55, edgecolor=BLACK,
               linewidth=0.9, hatch=hatch)
        ax.vlines(x[i], lo[i], hi[i], color=BLACK, lw=1.2)
        ax.text(x[i], hi[i]+12, f"{med[i]:.0f}", ha="center", fontsize=7)
        ax.text(x[i], -45, codes[i], ha="center", fontsize=9, fontweight="bold", color=BLACK)
    ax.axvline(1.5, color=MID, lw=1.0, ls=":")
    ax.axhline(150, color=BLACK, ls="--", lw=1.1)
    ax.text(3.35, 158, "±150 mm 条带", fontsize=7, color=DARK, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_ylim(0, 520)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_swmm.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig2_swmm.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig2] done")

# ============================================================
# 图 3（正文）：量测可实现性（灰阶 + 星号 + 阈值线）
# ============================================================
def fig3():
    obs = load("bm01_obs_estimator.json")
    ref = obs["reference_arms"]["mean"]["peak_mean_mm"]
    est = obs["estimated_arms"]
    arms = ["midpoint", "mp3", "mp5", "volbal"]
    labels = ["中点", "3点平均", "5点平均", "容积平衡"]
    vals = [est[a]["peak_mean_mm"] for a in arms]
    fcs = ["0.55", "0.70", "0.80", "0.95"]

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    style_ax(ax)
    x = np.arange(len(arms))
    for i, v in enumerate(vals):
        ax.bar(x[i], v, color=fcs[i], width=0.5, edgecolor=BLACK, linewidth=0.8)
        ax.text(i, v+2.5, f"{v:.0f}", ha="center", fontsize=7)
        ax.text(i, -12, "A1", ha="center", fontsize=8, fontweight="bold", color=BLACK)
    # 容积平衡：星号 + 0.0 mm
    ax.plot([3], [vals[3]+13], marker="*", ms=10, color=BLACK)
    ax.text(3, vals[3]+20, "与理想值偏差 0.0 mm", ha="center", fontsize=7, color=DARK)
    ax.axhline(ref, color=BLACK, ls="--", lw=1.2)
    ax.text(3.42, ref-6, "理想全池平均 75.0 mm", fontsize=7, color=DARK, ha="right")
    ax.axhline(150, color=BLACK, ls=":", lw=1.0)
    ax.text(3.42, 153, "±150 mm 判定阈值", fontsize=7, color=DARK, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_ylim(0, 170)
    ax.set_title("工程可实现估计量（均判 A1）", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_measurability.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig3_measurability.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig3] done")

# ============================================================
# 图 S1（补充）：网格敏感性（黑线 + 线型区分 + 注记）
# ============================================================
def figS1():
    sv = load("bm01_sve1d_compare.json")
    mean_v = [sv["sve"]["g120_mean"]["peak_mean_mm"], sv["sve"]["g240_mean"]["peak_mean_mm"]]
    gate_v = [sv["sve"]["g120_gate"]["peak_mean_mm"],
              sv["sve"]["g240_gate"]["peak_mean_mm"],
              sv["sve"]["g480_gate"]["peak_mean_mm"]]
    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    style_ax(ax)
    ax.plot([50.0, 25.0], mean_v, "-o", color=BLACK, lw=1.6, ms=5, mfc="none", label="全池平均反馈")
    ax.plot([50.0, 25.0, 12.5], gate_v, "--s", color=DARK, lw=1.6, ms=5, label="闸前单点反馈")
    ax.axhline(150, color=BLACK, ls=":", lw=1.0)
    ax.text(14, 158, "±150 mm 条带", fontsize=7, color=DARK)
    ax.axhline(72.74, color=MID, ls="-.", lw=1.0)
    ax.text(14, 62, "0D 集总 = 72.7 mm", fontsize=7, color=MID)
    ax.set_xlabel("网格尺寸 dx (m)　← 细化", fontsize=8)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_xlim(10, 55); ax.invert_xaxis()
    ax.annotate("随细化未收敛（网格敏感）", xy=(12.5, gate_v[-1]), xytext=(19, 190),
                fontsize=7, color=DARK, arrowprops=dict(arrowstyle="->", color=DARK, lw=0.8))
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figS1_grid.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS1_grid.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS1] done")

# ============================================================
# 图 S2（补充）：空间低通（实线解析 / 空心圆数值，控制书 Bode 感）
# ============================================================
def figS2():
    x = np.linspace(0.02, 2.0, 400)
    theory = np.sin(np.pi*x)**2 / (np.pi*x)
    w5 = load("bm01_w5_kc_identity.json")
    xs = [r["L_over_lambda"] for r in w5["rows"] if r["L_m"] == 6000]
    num = [r["numeric"] for r in w5["rows"] if r["L_m"] == 6000]
    fig, ax = plt.subplots(figsize=(3.8, 2.7))
    style_ax(ax)
    ax.plot(x, theory, color=BLACK, lw=1.8, ls="-", label="解析：sin²(πx)/(πx)")
    ax.plot(xs, num, "o", color=BLACK, ms=5.5, mfc="none", ls="none", label="1D SVE 数值（dx=25 m）")
    ax.axhline(0, color=MID, lw=0.8)
    ax.axvline(0.5, color=MID, ls=":", lw=1.0)
    ax.axhline(2/math.pi, color=MID, ls=":", lw=1.0)
    ax.plot([0.5], [2/math.pi], "s", color=BLACK, ms=5)
    ax.annotate("半功率点 x=0.5，响应=2/π≈0.637", xy=(0.5, 2/math.pi), xytext=(0.78, 0.86),
                fontsize=7, color=DARK, arrowprops=dict(arrowstyle="->", color=DARK, lw=0.8))
    ax.set_xlabel("x = L/λ（池长 / 扰动波长）", fontsize=8)
    ax.set_ylabel("归一化响应 |ĥ_mean| / δh_0", fontsize=8)
    ax.set_ylim(-0.15, 1.05)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figS2_lowpass.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS2_lowpass.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS2] done")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); figS1(); figS2()
    print("BOOK-STYLE FIGURES DONE ->", OUT)
