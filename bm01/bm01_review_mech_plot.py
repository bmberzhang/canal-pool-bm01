# -*- coding: utf-8 -*-
"""画 A4 机理时程图（Figure S6）：BM-01 主算例 gate 反馈 vs mean 反馈的四组时序对比。
输出 SVG（无 numpy/matplotlib，纯标准库手写），再用 Chrome 渲染 PNG。"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "raw", "bm01_review_mech_ts.json"), encoding="utf-8"))

# ---- 画布 ----
W, H = 900, 720
# 四个子图区域
PANELS = [
    ("pool-mean level (mm)", "dev_mean_mm", 0, 150),
    ("gate-adjacent local level (mm)", "dev_gate_mm", 0, 200),
    ("gate opening a (-)", "gate_a", 0, 1.0),
    ("outflow Q (m³/s)", "q_out", 160, 240),
]
PW, PH = 420, 240
GAPX, GAPY = 30, 30
ML, MT = 70, 40  # margin left/top

# 颜色：gate 红，mean 蓝（中国习惯涨红跌绿，这里用红/蓝区分两档）
C_GATE = "#d62728"
C_MEAN = "#1f77b4"

def series_for(arm):
    return D[arm]

# 时间范围统一 0-12h
TMAX = 12.0

def scale(v, vmin, vmax, lo, hi):
    return lo + (v - vmin) / (vmax - vmin) * (hi - lo)

parts = []
# 标题
parts.append(f'<text x="{W/2}" y="24" text-anchor="middle" font-size="20" font-family="DejaVu Sans" fill="#222">'
             f'Mechanism time series (BM-01, &#916;Q=+45 m&#179;/s, 12 h)</text>')

# 图例
lx = W - ML - 200
parts.append(f'<text x="{lx}" y="50" font-size="13" font-family="DejaVu Sans" fill="{C_GATE}">&#8212; gate feedback (A3)</text>')
parts.append(f'<text x="{lx+170}" y="50" font-size="13" font-family="DejaVu Sans" fill="{C_MEAN}">&#8212; pool-mean feedback (A1)</text>')

for pi, (label, key, vmin, vmax) in enumerate(PANELS):
    col = pi % 2
    row = pi // 2
    x0 = ML + col * (PW + GAPX)
    y0 = MT + 40 + row * (PH + GAPY)
    x1 = x0 + PW
    y1 = y0 + PH
    # 边框
    parts.append(f'<rect x="{x0}" y="{y0}" width="{PW}" height="{PH}" fill="none" stroke="#999" stroke-width="1"/>')
    # 轴
    parts.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#333" stroke-width="1"/>')
    parts.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#333" stroke-width="1"/>')
    # y 标签
    parts.append(f'<text x="{x0-8}" y="{y0+12}" text-anchor="end" font-size="12" font-family="DejaVu Sans" fill="#222">{label}</text>')
    # 刻度
    for tv in range(0, 13, 3):
        tx = scale(tv, 0, TMAX, x0, x1)
        parts.append(f'<text x="{tx}" y="{y1+16}" text-anchor="middle" font-size="11" fill="#555">{tv}h</text>')
    # 150mm 参考线（仅前两个子图）
    if key in ("dev_mean_mm", "dev_gate_mm"):
        yband = scale(150, vmin, vmax, y1, y0)
        parts.append(f'<line x1="{x0}" y1="{yband}" x2="{x1}" y2="{yband}" stroke="#aaa" stroke-dasharray="4,4" stroke-width="1"/>')
        parts.append(f'<text x="{x1-4}" y="{yband-4}" text-anchor="end" font-size="10" fill="#888">±150</text>')
    # 两条曲线
    for arm, color in (("gate", C_GATE), ("mean", C_MEAN)):
        s = D[arm]
        pts = []
        n = len(s["t_h"])
        step = max(1, n // 800)
        for i in range(0, n, step):
            t = s["t_h"][i]
            if t > TMAX:
                break
            v = s[key][i]
            v = max(vmin, min(vmax, v))
            px = scale(t, 0, TMAX, x0, x1)
            py = scale(v, vmin, vmax, y1, y0)
            pts.append(f"{px:.1f},{py:.1f}")
        parts.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="1.5"/>')

svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">' + "".join(parts) + "</svg>"
out = os.path.join(HERE, "figS6_mech_ts.svg")
open(out, "w", encoding="utf-8").write(svg)
print("saved", out)
