# -*- coding: utf-8 -*-
"""
(W8a') 主卖点网格对照汇总：单点 vs 池均离散度，跨 3 个算例
================================================================================
三模型第三轮 C 清单第 3 项：主卖点（"取样位置翻转判定码"）此前只有
L=6000/ΔQ=45 **一个** 算例的三档网格数据 → "n=1 硬伤"。

本脚本**不重跑水力**，而是汇总 W9 已真跑的三档网格结果
（dx = 50 / 25 / 12.5 m），把"池均档"与"单点档"的**网格离散度**并列成表：

  L=6000, ΔQ=45, 300 s  → 主案例（完整翻转）
  L=2000, ΔQ=45,  60 s  → 快周期（完整翻转）
  L=2000, ΔQ=45, 300 s  → 原周期（两臂都崩，对照组）

**预注册判据（先于汇总写死）**：
  P1：出现"池均三档全 A1 且单点三档全 A3"的**翻转算例 ≥ 2 个**；
  P2：每个翻转算例，单点档离散度 / 池均档离散度 ≥ 2。
  通过 = P1 ∧ P2 ⇒ 主卖点从 1 个算例扩到 ≥2 个算例，且"池均网格无关、
  单点网格病态"这一对照在两个独立算例上复现。

离散度定义：spread = (max − min) / min × 100%。
数据源：raw/bm01_w9_grid_check.json（真跑，wall 150.4 s）。
输出：raw/bm01_w8ap_grid_table.json、raw/bm01_w8ap_grid_table.svg、W8ap_网格对照表.md
"""
import json, os

SRC = os.path.join("raw", "bm01_w9_grid_check.json")
OUT_JSON = os.path.join("raw", "bm01_w8ap_grid_table.json")
OUT_SVG = os.path.join("raw", "bm01_w8ap_grid_table.svg")
OUT_MD = "W8ap_网格对照表.md"

ARM_CN = {"mean": "池均（蓄量等效）", "point": "闸前单点"}


def spread(vals):
    lo, hi = min(vals), max(vals)
    return (hi - lo) / lo * 100.0


def main():
    src = json.load(open(SRC, encoding="utf-8"))

    rows = []
    for tag, c in src["cases"].items():
        for arm in ("mean", "point"):
            vals = [r["arms"][arm]["peak_mean_mm"] for r in c["rows"]]
            codes = [r["arms"][arm]["judgement"] for r in c["rows"]]
            rows.append({
                "case": tag, "L_m": c["pool_len_m"], "dq": c["dq"],
                "period_s": c["ctrl_period_s"], "arm": arm,
                "dx50": round(vals[0], 2), "dx25": round(vals[1], 2),
                "dx125": round(vals[2], 2),
                "spread_pct": round(spread(vals), 2),
                "codes": codes,
                "codes_consistent": len(set(codes)) == 1,
            })

    # 翻转算例：池均三档全 A1 且单点三档全 A3
    flips = []
    for c in src["cases"].values():
        m = [r["arms"]["mean"]["judgement"] for r in c["rows"]]
        p = [r["arms"]["point"]["judgement"] for r in c["rows"]]
        if set(m) == {"A1"} and set(p) == {"A3"}:
            flips.append((c["pool_len_m"], c["ctrl_period_s"]))

    # 每个翻转算例的离散度比
    ratios = {}
    for f in flips:
        L, per = f
        key = "L%g_%gs" % (L, per)
        sp = {r["arm"]: r["spread_pct"] for r in rows
              if r["L_m"] == L and r["period_s"] == per}
        ratios[key] = round(sp["point"] / sp["mean"], 2) if sp["mean"] > 0 else None

    P1 = len(flips) >= 2
    P2 = all((v is not None and v >= 2.0) for v in ratios.values())
    passed = bool(P1 and P2)

    out = {
        "purpose": "W8a' 主卖点网格对照汇总（单点 vs 池均，跨 3 算例）",
        "preregistered": {"P1_翻转算例数_ge": 2, "P2_离散度比_ge": 2.0},
        "rows": rows,
        "flip_cases": [{"L_m": L, "ctrl_period_s": p} for L, p in flips],
        "flip_case_count": len(flips),
        "spread_ratio_point_over_mean": ratios,
        "P1_pass": P1, "P2_pass": P2, "W8ap_pass": passed,
    }
    out["verdict"] = (
        "W8a' **通过**：翻转算例 %d 个（≥2），翻转算例的『单点/池均』离散度比 %s（均 ≥2）"
        " ⇒ 主卖点由 1 个算例扩至 %d 个独立算例，且『池均网格无关 vs 单点网格病态』复现。"
        % (len(flips), ratios, len(flips)) if passed else
        "W8a' **不通过**：翻转算例数 %d 或离散度比不足。" % len(flips))

    os.makedirs("raw", exist_ok=True)
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    # ---- Markdown 表 ----
    lines = ["# W8a′ 网格对照表：池均 vs 单点（三档网格 dx=50/25/12.5 m）", "",
             "数据源：`raw/bm01_w9_grid_check.json`（真跑）。spread=(max−min)/min×100%。", "",
             "| 算例 | 反馈量 | dx=50 m | dx=25 m | dx=12.5 m | **spread** | 三档判定码 |",
             "|---|---|---|---|---|---|---|"]
    for tag, c in src["cases"].items():
        for arm in ("mean", "point"):
            r = next(x for x in rows if x["case"] == tag and x["arm"] == arm)
            lines.append("| %s<br>L=%g, ΔQ=%g, T=%g s | %s | %.1f | %.1f | %.1f | **%.2f%%** | %s |"
                         % (tag, r["L_m"], r["dq"], r["period_s"], ARM_CN[arm],
                            r["dx50"], r["dx25"], r["dx125"], r["spread_pct"],
                            "/".join(r["codes"])))
    lines += ["", "## 预注册判据与裁决", "",
              "- P1（翻转算例 ≥2）：**%s** — 实际 %d 个（%s）"
              % ("通过" if P1 else "不通过", len(flips),
                 "、".join("L=%g/T=%gs" % f for f in flips) or "无"),
              "- P2（每个翻转算例 单点/池均 spread ≥2）：**%s** — %s"
              % ("通过" if P2 else "不通过", ratios),
              "", "**裁决：%s**" % ("通过" if passed else "不通过"), "",
              "---", "", out["verdict"]]
    open(OUT_MD, "w", encoding="utf-8").write("\n".join(lines))

    # ---- SVG 条形图（纯标准库，无 numpy）----
    def esc(s):
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    W, H = 760, 430
    ml, mr, mt, mb = 78, 24, 54, 96
    pw, ph = W - ml - mr, H - mt - mb
    ymax = 130.0
    cases_all = list(src["cases"].items())
    ngroup = len(cases_all)
    gw = pw / ngroup
    bw = gw * 0.26
    bar_colors = {"mean": "#2f7ed8", "point": "#d1442f"}
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="Segoe UI, Microsoft YaHei, sans-serif">' % (W, H),
           '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
           '<text x="%d" y="26" font-size="17" font-weight="700" fill="#1a1a1a">三档网格下的峰值离散度：池均 vs 闸前单点</text>' % ml,
           '<text x="%d" y="46" font-size="12" fill="#666666">spread=(max−min)/min×100%%；判定码 A1=可吸收 / A3=不可吸收</text>' % ml]
    # y 轴网格与刻度
    for i in range(0, 6):
        v = ymax * i / 5.0
        y = mt + ph - ph * v / ymax
        svg.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e2e2e2" stroke-width="1"/>' % (ml, y, ml + pw, y))
        svg.append('<text x="%d" y="%.1f" font-size="11" fill="#888888" text-anchor="end">%.0f%%</text>' % (ml - 8, y + 4, v))
    svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#999999" stroke-width="1.2"/>' % (ml, mt, ml, mt + ph))
    svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#999999" stroke-width="1.2"/>' % (ml, mt + ph, ml + pw, mt + ph))
    for gi, (tag, c) in enumerate(cases_all):
        gx = ml + gw * gi
        for ai, arm in enumerate(("mean", "point")):
            r = next(x for x in rows if x["case"] == tag and x["arm"] == arm)
            v = min(r["spread_pct"], ymax)
            bh = ph * v / ymax
            x = gx + gw * 0.22 + ai * (bw + gw * 0.06)
            y = mt + ph - bh
            svg.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="2" fill="%s"/>'
                       % (x, y, bw, bh, bar_colors[arm]))
            lab = ("%.1f%%" % r["spread_pct"]) if r["spread_pct"] < 100 else ("%.0f%%" % r["spread_pct"])
            svg.append('<text x="%.1f" y="%.1f" font-size="11.5" font-weight="700" fill="#1a1a1a" text-anchor="middle">%s</text>'
                       % (x + bw / 2, y - 5, lab))
            svg.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="#444444" text-anchor="middle">%s</text>'
                       % (x + bw / 2, mt + ph + 16, esc("池均" if arm == "mean" else "单点")))
            svg.append('<text x="%.1f" y="%.1f" font-size="10" fill="#777777" text-anchor="middle">%s</text>'
                       % (x + bw / 2, mt + ph + 30, "/".join(r["codes"])))
        svg.append('<text x="%.1f" y="%.1f" font-size="11.5" font-weight="600" fill="#333333" text-anchor="middle">%s</text>'
                   % (gx + gw / 2, mt + ph + 50, esc("L=%g  T=%gs" % (c["pool_len_m"], c["ctrl_period_s"]))))
        svg.append('<text x="%.1f" y="%.1f" font-size="10" fill="#999999" text-anchor="middle">%s</text>'
                   % (gx + gw / 2, mt + ph + 66, esc("ΔQ=%g m³/s" % c["dq"])))
    # 图例
    svg.append('<rect x="%d" y="%d" width="13" height="13" rx="2" fill="%s"/>' % (ml + pw - 268, mt + 6, bar_colors["mean"]))
    svg.append('<text x="%d" y="%d" font-size="12" fill="#333333">池均（蓄量等效）</text>' % (ml + pw - 251, mt + 17))
    svg.append('<rect x="%d" y="%d" width="13" height="13" rx="2" fill="%s"/>' % (ml + pw - 268, mt + 26, bar_colors["point"]))
    svg.append('<text x="%d" y="%d" font-size="12" fill="#333333">闸前单点</text>' % (ml + pw - 251, mt + 37))
    svg.append('</svg>')
    open(OUT_SVG, "w", encoding="utf-8").write("\n".join(svg))

    print(out["verdict"])
    print("\n翻转算例：%s" % ("、".join("L=%g/T=%gs" % f for f in flips) or "无"))
    print("离散度比（单点/池均）：%s" % ratios)
    print("已写入 %s、%s、%s" % (OUT_JSON, OUT_SVG, OUT_MD))


if __name__ == "__main__":
    main()
