# -*- coding: utf-8 -*-
"""
BM-01 降阶损失对照图（读 raw/bm01_sve1d_compare.csv，不重算物理）
生成 figures/bm01_lumped_vs_sve1d.png，并打印一份压缩序列供内嵌图表使用。
用法：PYTHONPATH=<pylibs-extra> python bm01_sve1d_plot.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "raw", "bm01_sve1d_compare.csv")
JSN = os.path.join(HERE, "raw", "bm01_sve1d_compare.json")
OUT = os.path.join(HERE, "figures", "bm01_lumped_vs_sve1d.png")

with open(CSV, encoding="utf-8-sig") as fh:
    rows = list(csv.DictReader(fh))
with open(JSN, encoding="utf-8") as fh:
    meta = json.load(fh)

t = [float(r["t_s"]) / 60.0 for r in rows]
lump = [float(r["lumped_dev_m"]) * 1000 for r in rows]
c_mean = [float(r["g120_gate_dev_mean_m"]) * 1000 for r in rows]
f_mean = [float(r["g240_gate_dev_mean_m"]) * 1000 for r in rows]
g_gate = [float(r["g120_gate_dev_gate_m"]) * 1000 for r in rows]
m_mean = [float(r["g120_mean_dev_mean_m"]) * 1000 for r in rows]
BAND = 150.0
ALARM = 300.0
STEP = 15.0


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    avail = {f.name for f in font_manager.fontManager.ttflist}
    for c in ["Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC"]:
        if c in avail:
            plt.rcParams["font.sans-serif"] = [c]
            break
    plt.rcParams["axes.unicode_minus"] = False

    lump_peak = meta["lumped_recorded"]["peak_mm"]
    mean_peak = meta["attribution"]["sve_mean_feedback_mm"]["120"]
    ext = meta["attribution"]["extrapolated_mm"]
    gate_peak = meta["sve"]["g120_gate"]["peak_gate_mm"]

    f, ax = plt.subplots(1, 3, figsize=(16.4, 4.9))
    f.suptitle("BM-01 降阶损失对照：集总蓄量模型 vs 一维圣维南"
               "（同参数·同 PI·ΔQ=+45 m³/s·量级示例）", fontsize=13, fontweight="bold")

    # (a) 判定稳健性
    a = ax[0]
    a.axhspan(-BAND, BAND, color="#e8f3ff", zorder=0)
    a.axhline(BAND, color="#1a7f37", lw=1.2)
    a.axhline(-BAND, color="#1a7f37", lw=1.2)
    a.axvline(STEP, color="#888", lw=1.0, ls="--")
    a.plot(t, lump, color="#0066cc", lw=2.0, label=f"集总 0D（主案例）峰值 {lump_peak:.0f} mm → A1")
    a.plot(t, m_mean, color="#8e44ad", lw=1.4, ls="--",
           label=f"1D 但反馈取池均水位 峰值 {mean_peak:.0f} mm（≈0D）")
    a.plot(t, f_mean, color="#1a7f37", lw=1.8,
           label=f"1D 池均水位（240 格）峰值 {meta['sve']['g240_gate']['peak_mean_mm']:.0f} mm")
    a.axhline(0, color="none")
    a.annotate("正常条带 ±150 mm", xy=(8, BAND), xytext=(8, BAND + 12),
               fontsize=9, color="#1a7f37")
    a.annotate(f"外推峰值 ~{ext:.0f} mm\n越过条带 → 判定翻成 A3",
               xy=(28, 158), xytext=(70, 120), fontsize=9,
               arrowprops=dict(arrowstyle="->", color="#555", lw=1.0))
    a.set_title("(a) 判定稳健性：池均水位峰值", fontsize=11)
    a.set_xlabel("时间 (min)"); a.set_ylabel("水位偏差 (mm)")
    a.set_xlim(0, max(t)); a.set_ylim(-60, 205)
    a.legend(fontsize=7.5, loc="upper right"); a.grid(alpha=.25)

    # (b) 传感器看到的是什么
    b = ax[1]
    b.axhline(BAND, color="#1a7f37", lw=1.0)
    b.axhline(-BAND, color="#1a7f37", lw=1.0)
    b.axvline(STEP, color="#888", lw=1.0, ls="--")
    b.plot(t, c_mean, color="#1a7f37", lw=1.7, label="池均水位（池子的平均状态）")
    b.plot(t, g_gate, color="#e67e22", lw=1.3,
           label=f"闸前局部水位（闸门/传感器看到的）峰值 {gate_peak:.0f} mm")
    b.set_title("(b) 闸门量到的是局部水位，不是池均水位", fontsize=11)
    b.set_xlabel("时间 (min)"); b.set_ylabel("水位偏差 (mm)")
    b.set_xlim(0, max(t)); b.legend(fontsize=8, loc="upper left"); b.grid(alpha=.25)

    # (c) 闸门动作
    c = ax[2]
    c.plot(t, [float(r["lumped_gate"]) for r in rows], color="#0066cc", lw=1.8,
           label="集总 0D 闸门开度")
    c.plot(t, [float(r["g120_gate_gate"]) for r in rows], color="#1a7f37", lw=1.8,
           label="1D SVE 闸门开度")
    c.axhline(float(meta["lumped_recorded"]["gate_final"]), color="#888", lw=1.0, ls="--",
              label=f"终态开度 {meta['lumped_recorded']['gate_final']:.3f}")
    c.set_title("(c) 闸门动作（同一 PI、同一采样周期 300 s、同一速率限）", fontsize=11)
    c.set_xlabel("时间 (min)"); c.set_ylabel("开度 (0~1)")
    c.set_xlim(0, max(t)); c.legend(fontsize=8, loc="upper right"); c.grid(alpha=.25)

    f.tight_layout(rect=[0, 0, 1, 0.92])
    f.savefig(OUT, dpi=150)
    print("图已写出：", OUT, os.path.getsize(OUT), "bytes")


def compact():
    """每 10 min 取一点、约 25 点，供内嵌图表使用。"""
    idx = [i for i, v in enumerate(t) if i % 20 == 0]
    ser = lambda arr: [round(arr[i], 1) for i in idx]
    print("LABELS=" + json.dumps([round(t[i]) for i in idx]))
    print("LUMP=" + json.dumps(ser(lump)))
    print("MMEAN=" + json.dumps(ser(m_mean)))
    print("FMEAN=" + json.dumps(ser(f_mean)))
    print("GGATE=" + json.dumps(ser(g_gate)))


if __name__ == "__main__":
    compact()
    try:
        fig()
    except Exception as e:                                  # noqa: BLE001
        print("绘图失败（不影响数据）：", type(e).__name__, e)
        sys.exit(0)
