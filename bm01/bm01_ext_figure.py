# -*- coding: utf-8 -*-
"""
BM-01 外部门锚：两联对照图（双面板，各自量程）
左 = 池均反馈臂（A1，网格无关）；右 = 闸前单点反馈臂（A3，网格病态）。
读 raw/bm01_ext_compare.json 作图；不需重跑水力。
"""
from __future__ import annotations
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def pts(rows, code, ctrl):
    rs = [r for r in rows if r["code"] == code and r["ctrl"] == ctrl]
    return sorted(rs, key=lambda r: -r["dx_m"])


def draw(rows, path):
    W, H = 940, 430
    mtop, mbot, gap, ml, mr = 56, 62, 60, 66, 26
    panel_w = (W - ml - mr - gap) / 2.0
    ph = H - mtop - mbot
    dx_max = 320.0
    yA, yB = 200.0, 900.0
    s = []
    s.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d"' % (W, H, W, H))
    s.append(' font-family="Segoe UI, Noto Sans SC, sans-serif">')
    s.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    s.append('<text x="%d" y="30" font-size="15.5" font-weight="bold" fill="#111">BM-01 外部门锚：本仓库 SVE-1D vs 外部 EPA SWMM（L=6000 m，ΔQ=+45 m³/s，12 h）</text>' % ml)

    panels = [
        (ml, "池均反馈臂（判定 A1 / 可吸收）", yA, "mean", "pool", "#1f6fd0", "#d97706"),
        (ml + panel_w + gap, "闸前单点反馈臂（判定 A3 / 不可吸收）", yB, "gate", "gate", "#1f6fd0", "#d97706"),
    ]
    for (x0, title, ymax, c_sve, c_swmm, col1, col2) in panels:
        pw = panel_w
        def X(dx):
            return x0 + pw * (1.0 - min(dx, dx_max) / dx_max)
        def Y(v):
            return mtop + ph * (1.0 - min(v, ymax) / ymax)
        # 边框
        s.append('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#fafbfc" stroke="#e3e6ea"/>' % (x0, mtop, pw, ph))
        # y 网格
        nt = 4
        for i in range(nt + 1):
            v = ymax * i / nt
            y = Y(v)
            s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eceff3"/>' % (x0, y, x0 + pw, y))
            s.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#667" text-anchor="end">%.0f</text>' % (x0 - 6, y + 3.5, v))
        # 150 阈值
        if ymax > 150:
            yb = Y(150)
            s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#d33" stroke-dasharray="6,4" stroke-width="1.6"/>' % (x0, yb, x0 + pw, yb))
            s.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#d33">150 mm 条带</text>' % (x0 + 6, yb - 5))
        else:
            s.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#999">（150 mm 条带在此量程之上）</text>' % (x0 + 6, mtop + 14))
        # 标题
        s.append('<text x="%.1f" y="%d" font-size="12.5" font-weight="bold" fill="#222" text-anchor="middle">%s</text>' % (x0 + pw / 2, mtop - 12, title))
        # 数据线
        for code, ctrl, col, lab in (("SVE-1D", c_sve, col1, "本仓库 SVE-1D"), ("EPA SWMM", c_swmm, col2, "外部 EPA SWMM")):
            rs = pts(rows, code, ctrl)
            if not rs:
                continue
            d = " ".join("%.1f,%.1f" % (X(r["dx_m"]), Y(r["peak_pool_mm"])) for r in rs)
            dash = "" if ctrl in ("mean", "pool") else ' stroke-dasharray="0"'
            s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.4"%s/>' % (d, col, dash))
            for r in rs:
                s.append('<circle cx="%.1f" cy="%.1f" r="4.2" fill="#fff" stroke="%s" stroke-width="2.2"/>' % (X(r["dx_m"]), Y(r["peak_pool_mm"]), col))
            # 末端标签
            last = rs[-1]
            s.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>' % (X(last["dx_m"]) - 6, Y(last["peak_pool_mm"]) - 8, col, lab))
        # x 轴标题
        s.append('<text x="%.1f" y="%d" font-size="11.5" fill="#333" text-anchor="middle">网格分辨率 dx（m）· 向右更细</text>' % (x0 + pw / 2, mtop + ph + 42))
        # x 刻度
        for dxv in (300, 200, 100, 50, 25):
            if dxv > dx_max:
                continue
            xx = X(dxv)
            s.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#cfd6dd"/>' % (xx, mtop + ph, xx, mtop + ph + 6))
            s.append('<text x="%.1f" y="%d" font-size="10" fill="#667" text-anchor="middle">%d</text>' % (xx, mtop + ph + 18, dxv))
    # 纵轴标题
    s.append('<text x="16" y="%d" font-size="11.5" fill="#333" transform="rotate(-90 16 %d)" text-anchor="middle">池均峰值（mm）</text>' % (mtop + ph / 2, mtop + ph / 2))
    s.append('</svg>')
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(s))


def main():
    d = json.load(open(os.path.join(HERE, "raw", "bm01_ext_compare.json"), encoding="utf-8"))
    out = os.path.join(HERE, "raw", "bm01_ext_compare2.svg")
    draw(d["rows"], out)
    print("SVG →", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
