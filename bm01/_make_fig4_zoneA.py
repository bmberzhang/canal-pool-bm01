# -*- coding: utf-8 -*-
"""因A 二维可吸收域边界图（正文图4）：起爆ΔQ 随执行器速率，按 rho 分层。
数据源 = C1 审计 bm01_c1_rho_audit.json（真跑，24 点，不重跑）。"""
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

RAW = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/raw"
OUT = r"D:/自主水网测试案例/BM-01_中线渠池_上游阶跃可吸收/figures/pub_bookcn"

with open(os.path.join(RAW, "bm01_c1_rho_audit.json"), encoding="utf-8") as f:
    data = json.load(f)

rows = data["rows"]
rates = [0.5, 1.0, 2.0, 4.0]

# 三档 rho 分组
rhos = [0.44, 0.88, 1.32]
rho_label = {0.44: "rho=0.44", 0.88: "rho=0.88", 1.32: "rho=1.32"}

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))

# 面板(a)：起爆 ΔQ 随执行器速率
ax = axes[0]
colors = {0.44: "#2b6cb0", 0.88: "#c53030", 1.32: "#2f855a"}
markers = {2000.0: "o", 6000.0: "s"}

for rho in rhos:
    xs, ys, ms = [], [], []
    for r in rows:
        if abs(r["rho"] - rho) > 1e-3:
            continue
        xs.append(r["rate_mult"])
        ys.append(r["onset_dq"])
        ms.append(markers[r["L_m"]])
    ax.scatter(xs, ys, color=colors[rho], marker="o", s=34,
               label=f"{rho_label[rho]}", zorder=3)
    # 连线（同 rho 内的 4 个速率档）
    ax.plot(xs, ys, color=colors[rho], lw=1.2, alpha=0.6, zorder=2)
    # 标注 L=2000 与 L=6000 分开的说明：把 L=6000 的点标方框
    for r in rows:
        if abs(r["rho"] - rho) > 1e-3:
            continue
        if r["L_m"] == 6000.0:
            ax.scatter([r["rate_mult"]], [r["onset_dq"]], facecolors="none",
                       edgecolors=colors[rho], marker="s", s=70, linewidths=1.2, zorder=4)

ax.set_xscale("log")
ax.set_xticks(rates)
ax.set_xticklabels(["0.5", "1", "2", "4"])
ax.minorticks_off()
ax.set_xlabel("执行器速率倍数 rate_mult（相对基准 1/3000 s$^{-1}$）")
ax.set_ylabel("起爆阶跃 ΔQ$_\\mathrm{onset}$ (m$^3$/s)")
ax.set_title("(a) 起爆阈值随执行器速率单调升高")
ax.grid(True, which="both", alpha=0.25, lw=0.5)
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.set_xlim(0.35, 5.5)

# 图例说明：圆点 = L=2000，方框 = L=6000
from matplotlib.lines import Line2D
handles, labels = ax.get_legend_handles_labels()
handles += [Line2D([], [], marker="o", ls="", color="gray", label="L=2000"),
            Line2D([], [], marker="s", ls="", mfc="none", mec="gray", label="L=6000")]
ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right")

# 面板(b)：跃升倍数 R 随执行器速率（对数轴，区分边际型/失稳型）
ax2 = axes[1]
for rho in rhos:
    xs, ys = [], []
    for r in rows:
        if abs(r["rho"] - rho) > 1e-3:
            continue
        xs.append(r["rate_mult"])
        ys.append(r["jump_ratio_R"])
    ax2.scatter(xs, ys, color=colors[rho], s=34, zorder=3, label=f"rho={rho}")
    ax2.plot(xs, ys, color=colors[rho], lw=1.2, alpha=0.6, zorder=2)

ax2.axhline(3.0, color="gray", ls="--", lw=0.8, alpha=0.7)
ax2.text(0.42, 4.2, "R=3（边际/失稳分界）", color="gray", fontsize=7.5)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xticks(rates)
ax2.set_xticklabels(["0.5", "1", "2", "4"])
ax2.minorticks_off()
ax2.set_xlabel("执行器速率倍数 rate_mult")
ax2.set_ylabel("跃升倍数 R（对数）")
ax2.set_title("(b) 执行器越快，R 从失稳回落至边际")
ax2.grid(True, which="both", alpha=0.25, lw=0.5)
ax2.set_xlim(0.35, 5.5)
ax2.legend(fontsize=8, frameon=False, loc="upper right")

fig.suptitle("", fontsize=1)
fig.tight_layout()

os.makedirs(OUT, exist_ok=True)
outpath = os.path.join(OUT, "fig4_zoneA_2d_v2.png")
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print("SAVED", outpath)
