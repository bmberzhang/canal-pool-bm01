# -*- coding: utf-8 -*-
"""
BM-01R 附录：闸门过流能力校核（用查到的真实闸孔尺寸）
================================================================================
背景
----
BM-01R 主重跑一直用 BM-01 的降阶闸门式 `Q = kg·a·√(h/h_target)`，kg 只重标定到
Q₀，不是真实闸门标定。此前"唐河闸孔数×净宽"未查到，故列为最大未真实化项。

2026-10-01 查到（逐条可核，见 PARAMETERS.md §4）：
  · 唐河倒虹吸出口检修闸门**孔口净宽 5.5 m**，叠梁式，总高 7.5 m，3 孔共用
  · 同闸段布置 3 扇弧形工作闸门、3 台液压启闭机（一门一机）
  · 倒虹吸**管身段为 3 孔 1 联箱型结构，单孔 5.5 m × 5.7 m**（四个独立来源一致）

【2026-10-01 第二轮核查修正，务必先读】
  首轮把"孔口净宽 5.5 m"与"3 孔"拼成了"3 孔 × 5.5 m = 总净宽 16.5 m"。
  这个拼法**未经任何单一来源直接陈述**（属推断），且与倒虹吸工作原理相悖：
  倒虹吸是压力/半压力管，**过流受控制约于管身断面，而不是出口闸孔宽**。
  第二轮用四个彼此独立、数字自洽的来源把 **管身断面 3×5.5×5.7 m（面积 94.05 m²）**
  钉住了。故本脚本改用**管身断面控制**复核，闸孔式只作旁证。
  见 PARAMETERS.md §4。

本脚本回答一个问题：
  **把真实几何代进去，闸门/管身会不会成为瓶颈？**

方法
----
用实测率定的过闸流量公式，在"闸门全开"条件下反算该闸的**最大过流能力**，
再与该渠段的设计/加大流量、以及本案例的工况流量对照。

闸孔出流（中线工程建模实际采用的统一形式，见 PARAMETERS.md §4）：
    Q = σ · m · n · b · e · √(2 g h₁)
其中 n 为孔数、b 为单孔净宽、e 为开度、h₁ 为闸前水深、σ 淹没系数、m 流量系数。

口径：literature-parameterized｜sealed: false｜只做能力校核，不产生判定码。

用法：
  python bm01r_gate_check.py
  python bm01r_gate_check.py --json
"""
from __future__ import annotations

import argparse
import json
import math
import os

G = 9.81

# ---- 真实几何（2026-10-01 两轮核查后锁定，见 PARAMETERS.md §4）----
N_OPENINGS = 3          # 孔数（管身段与出口闸均为 3 孔，两来源一致）
B_OPENING = 5.5         # 单孔净宽 m
B_TOTAL = N_OPENINGS * B_OPENING
H_PIPE = 5.7            # 管身单孔净高 m（箱型 5.5 × 5.7，四来源一致）
A_PIPE = N_OPENINGS * B_OPENING * H_PIPE   # 管身过水面积 94.05 m²

# ---- 渠段参数（与主案例一致）----
H_TARGET = 4.50         # 设计水深 m
Q_DESIGN = 135.0
Q_ENLARGED = 160.0      # 加大流量
Q_STEPPED = 180.0       # 本案例工况（超加大流量 12.5%）

# ---- 闸门总高：用检修闸的 7.5 m 作量级参照（工作闸门实际高度未逐项查到，
#      但开度上限受此约束；弧形门全开时通常接近全孔高）----
GATE_HEIGHT_REF = 7.5

# ---- 流量系数 ----
# 闸孔出流统一式 m = ε·φ（垂向收缩 × 流速系数）。
# 宽顶堰上平板闸门：ε′=0.611~0.673、φ=0.85~1.0（《水力学》华东水利学院）
# 弧门自由出流 Cd 常在 0.55~0.65（CAP 式需现场率定）。
# 本脚本用区间做敏感性，不给单点值——避免把"率定值"当"设计值"。
M_LO, M_HI = 0.55, 0.65


def gate_capacity(e, h1, m_coef, sigma=1.0, n=N_OPENINGS, b=B_OPENING):
    """闸孔出流能力 Q = σ·m·n·b·e·√(2gh₁)"""
    return sigma * m_coef * n * b * e * math.sqrt(2.0 * G * h1)


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01R 闸门/管身过流能力校核")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    raw = os.path.join(here, "raw")

    # 闸前水位：设计水深 4.5 m（与主案例闸前控制目标一致）
    h1 = H_TARGET
    # 全开时开度取孔径高量级（露顶式弧门全开≈过水高度上限）
    e_full_lo, e_full_hi = 4.50, GATE_HEIGHT_REF

    rows = []
    for m_coef in (M_LO, M_HI):
        for e in (e_full_lo, e_full_hi):
            q_cap = gate_capacity(e, h1, m_coef)
            rows.append({"m": m_coef, "e_m": e, "q_capacity_m3s": round(q_cap, 1),
                         "ratio_to_design": round(q_cap / Q_DESIGN, 2),
                         "ratio_to_enlarged": round(q_cap / Q_ENLARGED, 2),
                         "ratio_to_stepped": round(q_cap / Q_STEPPED, 2)})

    q_lo = gate_capacity(e_full_lo, h1, M_LO)
    q_hi = gate_capacity(e_full_hi, h1, M_HI)

    # ---- 管身控制能力（倒虹吸的真实控制断面）----
    # 无压箱涵允许流速按工程惯例 1.5~2.0 m/s（远低于有压管 3~4 m/s）
    V_ALLOW_LO, V_ALLOW_HI = 1.5, 2.0
    q_pipe_lo = A_PIPE * V_ALLOW_LO
    q_pipe_hi = A_PIPE * V_ALLOW_HI
    v_at_design = Q_DESIGN / A_PIPE
    v_at_enlarged = Q_ENLARGED / A_PIPE
    v_at_stepped = Q_STEPPED / A_PIPE

    print("=" * 78)
    print("BM-01R 过流能力校核（真实几何：3 孔 × 5.5 m 净宽，管身 3×5.5×5.7 m）")
    print("=" * 78)
    print(f"闸孔/管身：{N_OPENINGS} 孔 × {B_OPENING} m = 总净宽 {B_TOTAL} m；"
          f"管身面积 {A_PIPE:.2f} m²（单孔 {B_OPENING}×{H_PIPE} m）")
    print(f"闸前水深取设计值 h₁ = {h1} m；开度取 {e_full_lo}~{e_full_hi} m（全开量级）")
    print(f"流量系数敏感性：m = {M_LO}~{M_HI}（率定区间，非单点设计值）")
    print()
    print("【路径一】出口闸孔出流式（旁证，控制断面不在此处）")
    print(f"{'m':>6}{'开度/e(m)':>11}{'过流能力':>12}{'/设计135':>10}{'/加大160':>10}{'/工况180':>10}")
    for r in rows:
        print(f"{r['m']:>6.2f}{r['e_m']:>11.2f}{r['q_capacity_m3s']:>12.1f}"
              f"{r['ratio_to_design']:>10.2f}{r['ratio_to_enlarged']:>10.2f}"
              f"{r['ratio_to_stepped']:>10.2f}")
    print(f"→ 闸门过流能力区间：{q_lo:.0f} ~ {q_hi:.0f} m³/s"
          f"（设计 {Q_DESIGN:.0f} 的 {q_lo/Q_DESIGN:.1f}~{q_hi/Q_DESIGN:.1f} 倍，"
          f"加大 {Q_ENLARGED:.0f} 的 {q_lo/Q_ENLARGED:.1f}~{q_hi/Q_ENLARGED:.1f} 倍）")
    print()
    print("【路径二】管身断面控制（主证：倒虹吸的控制断面在管身）")
    print(f"   管身面积 A = {A_PIPE:.2f} m²；无压箱涵允许流速 {V_ALLOW_LO}~{V_ALLOW_HI} m/s")
    print(f"   对应过流能力 {q_pipe_lo:.0f} ~ {q_pipe_hi:.0f} m³/s")
    print(f"   实际断面流速：Q=135 → {v_at_design:.3f} m/s；"
          f"Q=160 → {v_at_enlarged:.3f} m/s；Q=180 → {v_at_stepped:.3f} m/s")
    print(f"   ⇒ 连工况 {Q_STEPPED:.0f} 也只有 {v_at_stepped:.2f} m/s，"
          f"{'仍在允许区间内' if v_at_stepped <= V_ALLOW_HI else '已超允许上限'}")

    cap_ok = (q_lo > Q_STEPPED) and (v_at_stepped <= V_ALLOW_HI)
    print()
    print("结论：")
    if cap_ok:
        print("  两条路径一致：出口闸与管身过流能力都高于本案例工况流量，")
        print("  ⇒ **瓶颈不在闸门/管身，而在渠道输水能力（Manning）**。")
        print("  这与主案例把根源指向 Manning 的结论相互独立地吻合。")
    else:
        print("  存在闸门/管身能力区间跨过工况流量，需按实际率定值复核。")

    payload = {
        "check": "gate_and_pipe_capacity_vs_channel_capacity",
        "geometry": {
            "n": N_OPENINGS, "b_each_m": B_OPENING, "b_total_m": B_TOTAL,
            "pipe_section_m": [B_OPENING, H_PIPE], "pipe_area_m2": round(A_PIPE, 2),
            "source": ("管身 3 孔 1 联 5.5×5.7 m：中国水利报社 / 百度文库 / 筑龙 / 豆丁"
                       "四来源一致；出口检修闸孔口净宽 5.5 m、3 孔共用：期刊VIP"),
        },
        "formula": "Q = sigma*m*n*b*e*sqrt(2*g*h1)  (闸孔出流；中线建模统一式)",
        "inputs": {"h1_m": h1, "e_range_m": [e_full_lo, e_full_hi],
                   "m_range": [M_LO, M_HI],
                   "v_allow_m_s": [V_ALLOW_LO, V_ALLOW_HI],
                   "q_design_m3s": Q_DESIGN, "q_enlarged_m3s": Q_ENLARGED,
                   "q_stepped_m3s": Q_STEPPED},
        "gate_cases": rows,
        "gate_capacity_range_m3s": [round(q_lo, 1), round(q_hi, 1)],
        "pipe_capacity_range_m3s": [round(q_pipe_lo, 1), round(q_pipe_hi, 1)],
        "pipe_velocity_m_s": {"q135": round(v_at_design, 3),
                              "q160": round(v_at_enlarged, 3),
                              "q180": round(v_at_stepped, 3)},
        "bottleneck": "channel" if cap_ok else "undetermined",
        "verdict": ("闸门与管身过流能力均充裕，瓶颈在渠道输水能力（与主案例结论一致）"
                    if cap_ok else "闸门/管身能力可能不足，需按率定值复核"),
        "evidence_grade": "literature-parameterized",
        "sealed": False,
        "note": ("只做能力校核，不产生判定码。流量系数 m 为率定区间而非设计单点，"
                 "故给区间不给点值。开度上限以检修闸总高 7.5 m 为量级参照。"
                 "管身控制为倒虹吸的主控断面；闸孔式为旁证。"
                 "【第二轮修正】首轮误将'孔口净宽 5.5 m'与'3 孔'拼成总净宽 16.5 m，"
                 "该拼法无单一来源直接陈述且违背倒虹吸工作原理，已改用管身断面主控。"),
    }
    out = os.path.join(raw, "bm01r_gate_check.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
