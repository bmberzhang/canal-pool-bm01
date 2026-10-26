# -*- coding: utf-8 -*-
"""M3 图：判定阈值敏感性——两个反馈档的峰值分布 vs 阈值，展示翻转稳健性。
数据源 bm01_m3_band_sensitivity.json（零重跑）。"""
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

d = json.load(open(os.path.join(RAW, "bm01_tune_valley.json"), encoding="utf-8"))
bands = [100, 125, 150, 175, 200]

# 收集两档全部峰值
def collect(d, arm):
    ps = []
    for kp, kis in d["maps"][arm]["grid"].items():
        for ki, cell in kis.items():
            pk = cell.get("peak_mean_mm")
            if pk is not None:
                ps.append(pk)
    return ps

mean_pk = collect(d, "mean")
gate_pk = collect(d, "gate_avg")

fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))

# 面板(a)：两档峰值分布直方图 + 阈值线
ax = axes[0]
bins = [50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400]
ax.hist(gate_pk, bins=bins, alpha=0.55, color="#c53030", label="闸前100m平均（gate_avg）", edgecolor="white")
ax.hist(mean_pk, bins=bins, alpha=0.55, color="#2b6cb0", label="全池平均（mean）", edgecolor="white")
ax.axvline(150, color="0.3", ls="--", lw=1.0)
ax.text(151, ax.get_ylim()[1] * 0.95, "150mm", color="0.3", fontsize=8)
ax.axvline(154.23, color="#c53030", ls=":", lw=1.0)
ax.text(155, ax.get_ylim()[1] * 0.82, "154.2mm\n(gate_avg谷底)", color="#c53030", fontsize=7.5, va="top")
ax.set_xlabel("池均峰值 (mm)")
ax.set_ylabel("增益格数")
ax.set_title("(a) 两档峰值分布整体分离")
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.set_xlim(50, 410)

# 面板(b)：A1 占比随阈值
ax2 = axes[1]
a1_mean = [sum(1 for p in mean_pk if p <= b) / len(mean_pk) for b in bands]
a1_gate = [sum(1 for p in gate_pk if p <= b) / len(gate_pk) for b in bands]
ax2.plot(bands, a1_mean, "o-", color="#2b6cb0", lw=1.6, label="mean（全池平均）")
ax2.plot(bands, a1_gate, "s-", color="#c53030", lw=1.6, label="gate_avg（闸前100m平均）")
# 翻转占比差阴影（两线之间）
ax2.fill_between(bands, a1_gate, a1_mean, alpha=0.15, color="0.4")
ax2.axvline(150, color="0.3", ls="--", lw=0.9)
ax2.text(151, 0.92, "150mm", color="0.3", fontsize=8)
ax2.set_xlabel("判定阈值 (mm)")
ax2.set_ylabel("A1 判定占比")
ax2.set_title("(b) 翻转在阈值放宽到谷底之上才收敛")
ax2.set_ylim(0, 1.05)
ax2.set_xticks(bands)
ax2.grid(True, alpha=0.25, lw=0.5)
ax2.legend(fontsize=8, frameon=False, loc="center right")

fig.tight_layout()
os.makedirs(OUT, exist_ok=True)
outpath = os.path.join(OUT, "fig5_band_sensitivity.png")
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print("SAVED", outpath)
