# -*- coding: utf-8 -*-
"""
BM-01 一维 v2 断面形状敏感性图（读 raw/bm01_sve1d_v2.json，不重算物理）
生成 figures/bm01_sect_shape_sensitivity.png
用法：PYTHONPATH=<pylibs-extra> python bm01_sve1d_v2_plot.py
"""
import json
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
JSN = os.path.join(HERE, "raw", "bm01_sve1d_v2.json")
OUT = os.path.join(HERE, "figures", "bm01_sect_shape_sensitivity.png")

with open(JSN, encoding="utf-8") as fh:
    d = json.load(fh)

BASE = d["baseline_rect"]["converged_mm"]
LUMP = d["lumped_recorded_mm"]
BAND = 150.0
rows = d["switch_results"]
M = [0.0] + [r["switches"]["m_slope"] for r in rows]
PEAK = [BASE] + [r["converged_mm"] for r in rows]
# 更正前（矩形反解纵坡 + Richardson 外推）留档值，用于对照
LEGACY = [d["baseline_rect"].get("extrapolated_mm_legacy")] + \
         [r.get("legacy_rect_slope", {}).get("converged_mm") for r in rows]


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

    f, ax = plt.subplots(1, 2, figsize=(15.0, 4.6))

    # (a) 收敛峰值随边坡变化
    a = ax[0]
    a.axhspan(-10, BAND, color="#eef7ee", zorder=0)
    a.axhline(BAND, color="#1a7f37", lw=1.2, ls="--")
    a.axhline(LUMP, color="#0066cc", lw=1.4, ls=":")
    a.plot(M, LEGACY, "s--", color="#b0b0b0", lw=1.2, ms=5,
           label="更正前（矩形反解纵坡 + Richardson 外推）")
    a.plot(M, PEAK, "o-", color="#e67e22", lw=2.0, ms=7,
           label="1D 圣维南 闸前反馈（最细网格收敛值，前 2h 窗）")
    for m, pk in zip(M, PEAK):
        a.annotate(f"{pk:.0f}", xy=(m, pk), xytext=(0, 9), textcoords="offset points",
                   ha="center", fontsize=9, color="#a0522d")
    a.annotate("正常条带 +150 mm", xy=(0.02, BAND), xytext=(0.02, BAND + 6),
               fontsize=9, color="#1a7f37")
    a.annotate(f"矩形 {BASE:.0f} mm：越过条带 {BASE/BAND-1:+.1%} → A3（裕度已收窄）",
               xy=(0.04, BASE), xytext=(0.30, 196),
               fontsize=9, color="#a32d2d",
               arrowprops=dict(arrowstyle="->", color="#a32d2d", lw=0.9))
    a.annotate(f"集总 0D 记录值 {LUMP:.0f} mm（主案例判定依据）",
               xy=(2.0, LUMP - 2), xytext=(1.18, 26), fontsize=9, color="#0066cc",
               arrowprops=dict(arrowstyle="->", color="#0066cc", lw=0.9))
    a.annotate("梯形 m=2.5 → A1", xy=(2.48, PEAK[-1] + 8), xytext=(2.03, 44),
               fontsize=9, color="#1a7f37",
               arrowprops=dict(arrowstyle="->", color="#1a7f37", lw=0.9))
    a.set_title("(a) 断面形状决定判定：矩形 A3 → 梯形退回 A1", fontsize=11)
    a.set_xlabel("边坡系数 m（0 = 矩形）"); a.set_ylabel("池均水位峰值（mm）")
    a.set_xlim(-0.2, 2.8); a.set_ylim(-10, 215)
    a.legend(fontsize=8, loc="lower left"); a.grid(alpha=.25)

    # (b) 各断面的蓄量系数 dA/dh
    b = ax[1]
    W = 50.0
    HT = 7.0
    for m, col, mk in zip(M, ["#0066cc", "#7f77dd", "#1a7f37", "#e67e22"], ["o", "s", "D", "^"]):
        bB = W + 2 * m * HT
        bb = W
        b.plot([0, 1], [bb, bB], marker=mk, color=col, lw=1.8, ms=6,
               label=f"m={m:g}：水面宽 {bb:.0f} → {bB:.0f} m（+{bB/bb-1:+.0%}）")
    b.set_title("(b) 为什么梯形压峰值：水位越高水面越宽，蓄量涨得更快", fontsize=11)
    b.set_xlabel("相对水位抬升（0 = 目标水位 7.00 m）"); b.set_ylabel("水面宽 B（m）")
    b.set_xlim(-0.05, 1.1); b.set_ylim(44, 100)
    b.legend(fontsize=8, loc="upper left"); b.grid(alpha=.25)

    f.suptitle("BM-01 一维模型 v2：断面形状对峰值的影响"
               "（同参数·同 PI·ΔQ=+45 m³/s·量级示例）", fontsize=12.5, fontweight="bold")
    f.tight_layout(rect=[0, 0, 1, 0.90])
    f.savefig(OUT, dpi=150)

    # 自然宽度 = 15 cm（Word 版心），改写 pHYs 让转换器按 15 cm 插
    with open(OUT, "rb") as fh:
        data = bytearray(fh.read())
    w = struct.unpack(">I", bytes(data[16:20]))[0]
    dpi = int(round(w / (15.0 / 2.54)))
    i = 8
    while i < len(data):
        ln = struct.unpack(">I", bytes(data[i:i + 4]))[0]
        typ = bytes(data[i + 4:i + 8])
        if typ == b"pHYs":
            data[i + 8:i + 17] = struct.pack(">IIB", dpi, dpi, 1)
            break
        i += 12 + ln
        if typ == b"IDAT":
            break
    with open(OUT, "wb") as fh:
        fh.write(bytes(data))
    print("图已写出：", OUT, os.path.getsize(OUT), "bytes",
          f"（{w}px @ {dpi} dpi → {w/dpi*2.54:.1f} cm）")


if __name__ == "__main__":
    try:
        fig()
    except Exception as e:                                  # noqa: BLE001
        print("绘图失败：", type(e).__name__, e)
