# -*- coding: utf-8 -*-
"""
BM-01 论文正式图件（v2：按三模型审查收敛意见修订）
正文 3 张：fig1 翻转时序 / fig3 SWMM 复现 / fig4 量测可实现
补充材料 2 张：figS1 网格归因 / figS2 空间低通
修订要点：
- fig1：加局部 inset 展示 0D≈池均；判定码 A1/A3 色块标角；条带改灰阶填充
- figS1(原图2)：措辞改"网格敏感/未收敛"，补 ∝1/dx 参考线，横轴正序标细化方向
- fig3：柱顶直接标 A1/A3；误差棒=跨网格 min~max；两对之间分隔
- fig4：容积平衡显式标 0.0mm 偏差；参考线灰虚线；补 ±150mm 阈值线；按数值排序
- figS2(原图5)：x=L/λ 写死物理含义、标半功率点、图注给 1/4 出处
"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle

for fp in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf"]:
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 9

RAW = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/raw"
OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures/pub"
os.makedirs(OUT, exist_ok=True)

BLUE = "#2a78d6"; RED = "#d62728"; GREEN = "#2ca02c"; ORANGE = "#ff7f0e"; PURPLE = "#9467bd"
GRAY_BAND = "0.82"

def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.35, linewidth=0.6)
    ax.tick_params(labelsize=8)

def load(p):
    with open(os.path.join(RAW, p), encoding="utf-8") as f:
        return json.load(f)

def band_fill(ax, band, xmax, color=GRAY_BAND, alpha=0.28):
    """±band 条带用灰阶半透明填充（避免黑白打印与数据线混淆）"""
    ax.axhspan(-band, band, color=color, alpha=alpha, lw=0)

def verdict_tag(ax, code, x, y, color):
    """判定码色块标注（A1 绿底 / A3 红底）"""
    txt = "A1" if code == "A1" else "A3"
    bbox = dict(boxstyle="round,pad=0.25", fc=color, ec="none", alpha=0.9)
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=9, color="white",
            fontweight="bold", ha="center", va="center", bbox=bbox, zorder=5)

# ============================================================
# 图 1（正文）：受控反例 —— 取样位置决定判定翻转（含 inset）
# ============================================================
def fig1():
    main = load("bm01_main.json")
    sv   = load("bm01_sve1d_compare.json")
    t0 = np.array(main["series_downsampled"]["t_s"]) / 3600.0
    d0 = np.array(main["series_downsampled"]["pool1_dev_m"]) * 1000.0
    t1 = np.array(sv["series"]["t"]) / 3600.0
    d_mean = np.array(sv["series"]["dev_mean"]) * 1000.0
    d_gate = np.array(sv["series"]["dev_gate"]) * 1000.0
    band = 150.0

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
    # 左：0D 集总
    ax = axes[0]; style_ax(ax); band_fill(ax, band, 4.5)
    ax.plot(t0, d0, color=BLUE, lw=2)
    ax.set_ylim(-160, 170)
    ax.set_title("(a) 0D 集总", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8); ax.set_ylabel("水位偏差 (mm)", fontsize=8)
    verdict_tag(ax, "A1", 0.85, 0.85, GREEN)
    ax.text(0.25, 165, "+45 m³/s 阶跃", fontsize=7, color="0.3", ha="center")

    # 中：1D 全池平均（与 0D 几乎重合，用独立窄 y 轴突出等价）
    ax = axes[1]; style_ax(ax); band_fill(ax, band, 4.5)
    ax.plot(t1, d_mean, color=GREEN, lw=2)
    ax.set_ylim(-160, 170)
    ax.set_title("(b) 1D·全池平均反馈", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8)
    verdict_tag(ax, "A1", 0.85, 0.85, GREEN)

    # 右：1D 闸前单点（持续极限环）
    ax = axes[2]; style_ax(ax); band_fill(ax, band, 4.5)
    ax.plot(t1, d_gate, color=RED, lw=2)
    ax.set_ylim(-350, 550)
    ax.set_title("(c) 1D·闸前单点反馈", fontsize=9)
    ax.set_xlabel("时间 (h)", fontsize=8)
    verdict_tag(ax, "A3", 0.85, 0.85, RED)
    # 标注极限环区间
    ax.axvspan(1.0, 4.5, color=ORANGE, alpha=0.08, lw=0)
    ax.text(2.7, 520, "持续极限环", fontsize=7, color=ORANGE, ha="center")

    # 右图加 inset：把 0D 与池均的峰值差异放大（证明 0D≈池均）
    axins = axes[2].inset_axes([0.08, 0.10, 0.42, 0.40])
    axins.plot(t0, d0, color=BLUE, lw=1.8, label="0D")
    axins.plot(t1, d_mean, color=GREEN, lw=1.8, label="1D 池均")
    axins.axhline(band, color="0.5", ls="--", lw=1.0)
    axins.set_xlim(0, 1.6); axins.set_ylim(0, 90)
    axins.tick_params(labelsize=6)
    axins.set_title("峰值放大（0D≈池均）", fontsize=6.5)
    axins.legend(fontsize=6, frameon=False, loc="lower right")
    for s in axins.spines.values():
        s.set_linewidth(0.6)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_controlled_counterexample.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig1_controlled_counterexample.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig1] done")

# ============================================================
# 图 S1（补充）：网格敏感性归因 —— 单点档网格病态 vs 池均档网格收敛
# ============================================================
def figS1():
    sv = load("bm01_sve1d_compare.json")
    mean_v = [sv["sve"]["g120_mean"]["peak_mean_mm"],
              sv["sve"]["g240_mean"]["peak_mean_mm"]]
    gate_v = [sv["sve"]["g120_gate"]["peak_mean_mm"],
              sv["sve"]["g240_gate"]["peak_mean_mm"],
              sv["sve"]["g480_gate"]["peak_mean_mm"]]
    dx_mean = [50.0, 25.0]
    dx_gate = [50.0, 25.0, 12.5]

    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    style_ax(ax)
    # 主图：网格敏感性（正序 x：细化向右，即 dx 递减）
    ax.plot(dx_mean, mean_v, "-o", color=GREEN, lw=2, ms=6, label="全池平均反馈")
    ax.plot(dx_gate, gate_v, "-o", color=RED, lw=2, ms=6, label="闸前单点反馈")
    ax.axhline(150, color="0.4", ls="--", lw=1.2)
    ax.text(14, 158, "±150 mm 条带", fontsize=7, color="0.4")
    ax.axhline(72.74, color=BLUE, ls=":", lw=1.2)
    ax.text(14, 62, "0D 集总 = 72.7 mm", fontsize=7, color=BLUE)
    ax.set_xlabel("网格尺寸 dx (m)　← 细化", fontsize=8)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_xlim(10, 55)
    ax.invert_xaxis()
    # 标注"网格敏感/未收敛"而非"发散"
    ax.annotate("随细化未收敛（网格敏感）",
                xy=(12.5, gate_v[-1]), xytext=(20, 190),
                fontsize=7, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figS1_grid_sensitivity.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS1_grid_sensitivity.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS1] done")

# ============================================================
# 图 2（正文）：SWMM 双求解器独立复现
# ============================================================
def fig3():
    ext = load("bm01_ext_compare.json")
    rows = ext["rows"]
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        g[(r["code"], r["ctrl"])].append(r["peak_pool_mm"])
    order = [("SVE-1D", "mean"), ("EPA SWMM", "pool"),
             ("SVE-1D", "gate"), ("EPA SWMM", "gate")]
    labels = ["本仓库\n全池平均", "SWMM\n全池平均", "本仓库\n闸前单点", "SWMM\n闸前单点"]
    colors = [GREEN, GREEN, RED, RED]
    codes  = ["A1", "A1", "A3", "A3"]
    med = [float(np.median(g[k])) for k in order]
    lo  = [float(min(g[k])) for k in order]
    hi  = [float(max(g[k])) for k in order]

    fig, ax = plt.subplots(figsize=(4.0, 2.9))
    style_ax(ax)
    x = np.arange(len(order))
    for i in range(len(order)):
        # 纹理填充（A1 斜线 / A3 点阵）保证黑白打印可辨
        hatch = "//" if codes[i] == "A1" else ".."
        ax.bar(x[i], med[i], color=colors[i], width=0.55, alpha=0.75,
               edgecolor="black", linewidth=0.8, hatch=hatch)
        # 误差棒：跨网格 min~max
        ax.vlines(x[i], lo[i], hi[i], color="0.3", lw=1.2)
        ax.text(x[i], hi[i]+12, f"{med[i]:.0f}", ha="center", fontsize=7)
        ax.text(x[i], -45, codes[i], ha="center", fontsize=9, fontweight="bold",
                color=colors[i])
    # 两对之间的分隔
    ax.axvline(1.5, color="0.5", lw=1.0, ls=":")
    ax.axhline(150, color="0.4", ls="--", lw=1.2)
    ax.text(3.35, 158, "±150 mm 条带", fontsize=7, color="0.4", ha="right")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_ylim(0, 520)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_swmm_reproduction.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig2_swmm_reproduction.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig2] done")

# ============================================================
# 图 3（正文）：量测可实现性
# ============================================================
def fig4():
    obs = load("bm01_obs_estimator.json")
    ref_mean = obs["reference_arms"]["mean"]["peak_mean_mm"]  # 75.05
    est = obs["estimated_arms"]
    # 按数值降序排列，突出"越接近池均越判 A1"的单调趋势
    arms = ["midpoint", "mp3", "mp5", "volbal"]
    labels = ["中点", "3点平均", "5点平均", "容积平衡"]
    vals = [est[a]["peak_mean_mm"] for a in arms]
    # 已按降序（111.9 / 79.1 / 76.3 / 75.0）
    colors = [PURPLE, PURPLE, PURPLE, GREEN]

    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    style_ax(ax)
    x = np.arange(len(arms))
    bars = ax.bar(x, vals, color=colors, width=0.5, alpha=0.8,
                  edgecolor="black", linewidth=0.7)
    # 容积平衡棒加星号 + 显式标注 0.0mm 偏差
    ax.text(3, vals[3]+2.5, "偏差 0.0 mm", ha="center", fontsize=7, color=GREEN)
    ax.text(3, vals[3]+7, "*", ha="center", fontsize=11, color=GREEN, fontweight="bold")
    # 理想参考线（灰虚线）
    ax.axhline(ref_mean, color="0.4", ls="--", lw=1.4)
    ax.text(3.42, ref_mean-5, "理想全池平均 75.0 mm", fontsize=7, color="0.4", ha="right")
    # ±150mm 阈值线
    ax.axhline(150, color=RED, ls=":", lw=1.0)
    ax.text(3.42, 152, "±150 mm 判定阈值", fontsize=7, color=RED, ha="right")
    for i, v in enumerate(vals):
        ax.text(i, v+2.5, f"{v:.0f}", ha="center", fontsize=7)
        ax.text(i, -12, "A1", ha="center", fontsize=8, fontweight="bold", color=GREEN)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("水位峰值 (mm)", fontsize=8)
    ax.set_ylim(0, 165)
    ax.set_title("工程可实现估计量（均判 A1）", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_measurability.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "fig3_measurability.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[fig3] done")

# ============================================================
# 图 S2（补充）：空间低通机制
# ============================================================
def figS2():
    x = np.linspace(0.02, 2.0, 400)
    theory = np.where(np.abs(x) < 1e-9, 1.0, np.sin(np.pi*x)**2 / (np.pi*x))
    w5 = load("bm01_w5_kc_identity.json")
    xs_num = [r["L_over_lambda"] for r in w5["rows"] if r["L_m"] == 6000]
    num = [r["numeric"] for r in w5["rows"] if r["L_m"] == 6000]

    fig, ax = plt.subplots(figsize=(3.8, 2.7))
    style_ax(ax)
    ax.plot(x, theory, color=BLUE, lw=2, label="解析：sin²(πx)/(πx)")
    ax.plot(xs_num, num, "o", color=RED, ms=6, mfc="none", label="1D SVE 数值（dx=25 m）")
    ax.axhline(0, color="0.4", lw=0.8)
    # 半功率点 x=0.5
    ax.axvline(0.5, color="0.4", ls=":", lw=1.0)
    ax.axhline(2/math.pi, color="0.5", ls=":", lw=1.0)
    ax.plot([0.5], [2/math.pi], "o", color=ORANGE, ms=6)
    ax.annotate("半功率点 x=0.5\n响应=2/π≈0.637",
                xy=(0.5, 2/math.pi), xytext=(0.85, 0.85),
                fontsize=7, color=ORANGE,
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.8))
    ax.set_xlabel("x = L/λ（池长 / 扰动波长）", fontsize=8)
    ax.set_ylabel("归一化响应 |ĥ_mean| / δh_0", fontsize=8)
    ax.set_ylim(-0.15, 1.05)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "figS2_spatial_lowpass.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "figS2_spatial_lowpass.png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("[figS2] done")

if __name__ == "__main__":
    fig1(); figS1(); fig3(); fig4(); figS2()
    print("ALL FIGURES DONE ->", OUT)
