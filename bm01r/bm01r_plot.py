# -*- coding: utf-8 -*-
"""
BM-01R 主图：文献参数化把 BM-01 的判定带到了哪
================================================================================
三联图：
  (a) 池均水位偏差时序 —— 0D vs 1D，原参数 vs 真实参数
  (b) 渠道输水能力硬约束 —— 渠底坡度固定时，流量越大所需水深越大
  (c) 各工况峰值对照 —— 条带、Manning 平衡线、判定码

用法：PYTHONPATH=<pylibs-extra> python bm01r_plot.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "BM-01_中线渠池_上游阶跃可吸收")))

import bm01r_lit as R
import bm01r_converge as CV

JSN = os.path.join(HERE, "raw", "bm01r_literature.json")
OUT = os.path.join(HERE, "figures", "bm01r_literature_result.png")
BAND = 150.0

BLUE, RED, AMBER, GREEN, MUTED = "#0066cc", "#a32d2d", "#ba7517", "#1a7f37", "#5f5e5a"


def zero_d_series(cfg):
    import bm01_pool_step as M0
    p = M0.Params(
        pool_len_m=cfg["pool_len_m"], pool_width_m=cfg["water_width_m"],
        q0=cfg["q0"], h_target=cfg["h_target"], a0=cfg["a0"], kg=cfg["kg"],
        duration_s=cfg["duration_s"], t_step_s=cfg["t_step_s"], step_dq=cfg["step_dq"],
    )
    sc = M0.Scenario(key="k", label="l", pools=1, step_dq=cfg["step_dq"],
                     feedforward=False, note="")
    r = M0.simulate(p, sc)
    return r["series"]["t"], r["series"]["h1_dev"]


def main() -> int:
    with open(JSN, encoding="utf-8") as fh:
        d = json.load(fh)

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

    orig, lit = R.make_orig(), R.make_lit()
    sec = R.SVE2.Section(lit["bottom_m"], lit["m_slope"])
    S0 = R.SVE2.section_normal_slope(sec, lit["h_target"], lit["q0"], lit["n"])
    hn_new = CV.normal_depth(lit["q0"] + lit["step_dq"], lit["bottom_m"],
                             lit["m_slope"], lit["n"], S0)
    rise_req = (hn_new - lit["h_target"]) * 1000.0

    f = plt.figure(figsize=(15.0, 8.2))
    gs = f.add_gridspec(2, 2, height_ratios=[1.0, 0.92], hspace=0.38, wspace=0.22)

    # ---------------- (a) 时序 ----------------
    a = f.add_subplot(gs[0, 0])
    a.axhspan(-60, BAND, color="#eef7ee", zorder=0)
    a.axhline(BAND, color=GREEN, lw=1.0, ls="--", zorder=2)
    a.text(0.4, BAND + 8, f"正常调节条带 +{BAND:.0f} mm", fontsize=7.6, color=GREEN)

    for cfg, tag, col, ls, lab in (
        (orig, "A0", BLUE, ":", "0D 集总 · 原参数"),
        (lit, "A1", BLUE, "-", "0D 集总 · 真实参数"),
    ):
        t, dev = zero_d_series(cfg)
        a.plot([x / 3600 for x in t], [x * 1000 for x in dev], color=col, ls=ls,
               lw=1.2, label=lab, zorder=3)

    for key, col, ls, lab in (
        ("B0", RED, ":", "1D 圣维南 · 原参数（矩形）"),
        ("B2", RED, "-", "1D 圣维南 · 真实参数（梯形）"),
    ):
        s = d["runs"][key]["series_sampled"]
        a.plot([x / 3600 for x in s["t_s"]], s["dev_mean_mm"], color=col, ls=ls,
               lw=1.3, label=lab, zorder=4)

    a.axhline(rise_req, color=AMBER, lw=1.1, ls="-.", zorder=2)
    a.text(11.6, rise_req + 14, f"渠道输水能力要求的平衡水深 +{rise_req:.0f} mm",
           fontsize=7.4, color=AMBER)
    a.annotate("0D：峰值 73 / 45 mm\n之后回到零附近",
               xy=(17.0, 6), xytext=(13.4, 168), fontsize=7.4, color=BLUE,
               linespacing=1.4,
               arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8))
    a.annotate("1D 真实参数：稳定在 +330 mm，回不来",
               xy=(19.0, 331), xytext=(9.2, 470), fontsize=7.8, color=RED,
               arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    a.set_title("(a) 池均水位偏差：0D 回到零附近，1D 停在一个回不来的新水位", fontsize=9.6)
    a.set_xlabel("时间 (h)", fontsize=8.4)
    a.set_ylabel("池均水位偏差 (mm)", fontsize=8.4)
    a.set_xlim(0, 24)
    a.set_ylim(-60, 800)
    a.legend(fontsize=7.0, loc="upper left", framealpha=0.92)
    a.grid(alpha=.22)
    a.tick_params(labelsize=7.4)

    # ---------------- (b) 渠道输水能力 ----------------
    b = f.add_subplot(gs[0, 1])
    qs = [lit["q0"] + i for i in range(0, 76, 1)]
    hs = [CV.normal_depth(q, lit["bottom_m"], lit["m_slope"], lit["n"], S0) for q in qs]
    b.plot(qs, [(x - lit["h_target"]) * 1000 for x in hs], color=BLUE, lw=1.6)
    b.axhline(0, color=MUTED, lw=0.6)
    b.axhspan(0, BAND, color="#eef7ee", zorder=0)
    b.axhline(BAND, color=GREEN, lw=1.0, ls="--")
    b.text(qs[1], BAND + 12, f"条带 +{BAND:.0f} mm", fontsize=7.6, color=GREEN)

    q_end = lit["q0"] + lit["step_dq"]
    b.plot([lit["q0"], q_end], [0, rise_req], "o", color=RED, ms=5, zorder=5)
    b.annotate("", xy=(q_end, rise_req), xytext=(q_end, 0),
               arrowprops=dict(arrowstyle="<->", color=RED, lw=1.1))
    b.text(q_end + 2.5, rise_req * 0.52,
           f"+{rise_req:.0f} mm\n= 条带的 {rise_req / BAND:.1f} 倍",
           fontsize=7.8, color=RED, linespacing=1.4)
    b.annotate(f"设计工况\n{lit['q0']:.0f} m³/s", xy=(lit["q0"], 0),
               xytext=(lit["q0"] + 8, 320), fontsize=7.4, color=MUTED,
               linespacing=1.4,
               arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
    b.set_title("(b) 硬约束：渠底坡度不变，通过更大流量需要更大水深", fontsize=9.6)
    b.set_xlabel("流量 Q (m³/s)", fontsize=8.4)
    b.set_ylabel("相对设计水深的抬升 (mm)", fontsize=8.4)
    b.set_ylim(-40, rise_req * 1.16)
    b.grid(alpha=.22)
    b.tick_params(labelsize=7.4)

    # ---------------- (c) 峰值对照 ----------------
    c = f.add_subplot(gs[1, :])
    order = ["A0", "A1", "B0", "B1", "B2", "B3", "B4"]
    names = ["A0\n0D·原", "A1\n0D·真实", "B0\n1D·原(矩形)", "B1\n1D·真实尺度\n(矩形)",
             "B2\n1D·真实(梯形)\nQ=180·超限", "B3\n1D·真实(梯形)\n+25%·超限",
             "B4\n1D·真实(梯形)\nQ=160·合规上限"]
    vals = []
    for k in order:
        v = d["runs"][k].get("extrapolated_peak_mean_mm",
                             d["runs"][k].get("peak_mm"))
        vals.append(v)
    cols = [BLUE, BLUE, RED, RED, RED, RED, AMBER]
    cols[4] = "#7a1f1f"
    bars = c.bar(range(len(order)), vals, width=0.56, color=cols, alpha=0.88, zorder=3)
    for i, (bb, v) in enumerate(zip(bars, vals)):
        c.text(bb.get_x() + bb.get_width() / 2, v + 16, f"{v:.0f}", ha="center",
               fontsize=8.2, color=MUTED)
    c.axhline(BAND, color=GREEN, lw=1.2, ls="--", zorder=4)
    c.text(0.12, BAND - 44, f"条带 +{BAND:.0f} mm", fontsize=8.0, color=GREEN,
           ha="left", zorder=6)
    c.axhline(rise_req, color=AMBER, lw=1.1, ls="-.", zorder=2)
    c.text(6.42, rise_req + 16, f"Manning 平衡 +{rise_req:.0f} mm", fontsize=8.0,
           color=AMBER, ha="right")
    for i, k in enumerate(order):
        r = d["runs"][k]
        code = "A3" if vals[i] > BAND else "A1"
        c.text(i, -0.13 * max(vals) * 0.55 - 30, code, ha="center", fontsize=9.4,
               color=(GREEN if code == "A1" else RED), fontweight="bold")
    c.set_xticks(range(len(order)))
    c.set_xticklabels(names, fontsize=7.6, linespacing=1.5)
    c.set_title("(c) 峰值对照：真实渠段全部越带——连合规上限 Q=160（B4）也越带 29%",
                fontsize=9.6)
    c.set_ylabel("池均水位偏差峰值 (mm)", fontsize=8.4)
    c.set_ylim(-230, max(max(vals), rise_req) * 1.14)
    c.grid(alpha=.22, axis="y")
    c.tick_params(labelsize=7.4)

    f.savefig(OUT, dpi=170, bbox_inches="tight", facecolor="white")
    print(f"图已写出 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
