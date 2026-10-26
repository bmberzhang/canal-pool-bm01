# -*- coding: utf-8 -*-
"""
BM-01R：BM-01 的文献参数化重跑
================================================================================
把 BM-01 的"量级示例参数（assume）"换成南水北调中线**公开文献可核**的
真实渠段参数，回答一个问题：

    判定码 A1（扰动留在控制层、域内收敛）在真实几何下还成立吗？

为什么选这一段
--------------
唐河节制闸 → 放水河节制闸 渠段（桩号 1 046+216 ~ 1 071+783，长 25.567 km）
是石家庄以北**唯一无输水建筑物**的渠段——两水尺断面间无渡槽、倒虹吸、隧洞、
暗涵，也未设分水口门。因此它是本项目里"纯渠池"最干净的真实原型。

参数来源（逐条可核，详见 README.md 的参数来源表）
--------------------------------------------------
· 渠段长 25 567 m、底宽 18.5~23.0 m、边坡 2.0~2.5、设计水深 4.5 m、
  设计流量 135 m³/s、纵坡 1/25000~1/20000、糙率 0.0150
  —— 长江科学院院报《南水北调中线总干渠冰盖糙率观测分析》
· 闸前水位变幅推荐控制在 ±0.15 m
  —— 水利学报《南水北调中线工程冬季输水动态目标水位控制算法研究》
· 运行流速约 1 m/s、中线全线断面七种边坡比 —— 中线公司公开技术方案征集公告

三处自洽性校验（本脚本启动即跑，是"这些参数能不能用"的判据）
------------------------------------------------------------
1. 用 (底宽, 边坡, 水深, 糙率) 反解纵坡，应落在公开区间 1/25000~1/20000；
2. 断面流速 Q/A 应接近公开的"约 1 m/s"；
3. 0D 等效水面面积与原案例差多少（决定 Δh≈ΔQ·Δt/A 的响应快慢）。

析因设计（这是本脚本唯一的方法论要点）
--------------------------------------
一次只动一件事，把"池长/水深/流量变了"与"断面形状变了"分开：

  A0  原参数 · 0D          —— 复现原案例判定（自检）
  A1  真实参数 · 0D        —— 真实参数下 0D 会说什么
  B0  原参数 · 1D 矩形     —— 与 v1 基准对齐
  B1  真实参数 · 1D 矩形   —— 隔离"尺度效应"（池长 4.3×、水面面积 3.5×）
  B2  真实参数 · 1D 梯形   —— 真实几何（含断面形状效应）
  B3  真实参数 · 1D 梯形 +33.75 m³/s —— 同比例扰动（+25%，与主案例同强度）

网格一律两级 + 一阶 Richardson 外推，避免把网格误差当物理效应。

口径（不得外推）
----------------
证据等级 `literature-parameterized`：参数取自**公开文献**，不等于官方活动
R12-standard-benchmarks-r1 的冻结题单值。`sealed: false`，不得声称"通过官方测试"，
不得与官方活动拼接或替代官方结论。

用法
----
  python bm01r_lit.py             # 全量重跑
  python bm01r_lit.py --quick     # 只跑单级网格（快，用于冒烟）
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict

HERE = os.path.dirname(os.path.abspath(__file__))
BM01_DIR = os.path.abspath(os.path.join(HERE, "..", "BM-01_中线渠池_上游阶跃可吸收"))
sys.path.insert(0, BM01_DIR)

import bm01_pool_step as M0           # 0D 集总模型（原案例唯一入口）
import bm01_sve1d_v2 as SVE2          # 1D 圣维南求解器（支持梯形 + 指定坡度）

CASE_ID = "BM-01R"
CASE_TITLE = "唐河—放水河渠段·文献参数化重跑"
OBJECT = "南水北调中线 唐河节制闸→放水河节制闸 渠段"
LOCAL_CAMPAIGN_ID = "local-bm01r-literature-r1"
EVIDENCE_GRADE = "literature-parameterized"
G = 9.81

SOURCE_NOTE = (
    "参数取自公开文献/公开公告（渠段长、底宽、边坡、设计水深、设计流量、纵坡、糙率、"
    "水位变幅推荐值）。**不等于**官方活动 R12-standard-benchmarks-r1 的冻结题单值；"
    "本包 sealed=false，仅为本地可复跑算例，不构成现场结论、不换算 L3/工程授权。"
)


# ======================================================================
# 两套参数
# ======================================================================
def make_orig() -> dict:
    """原 BM-01 参数：全部为量级示例（assume）。取自 raw/bm01_main.json 的字段语义。"""
    return {
        "tag": "A/B0",
        "label": "原 BM-01（量级示例）",
        "pool_len_m": 6000.0,
        "bottom_m": 50.0,
        "m_slope": 0.0,
        "water_width_m": 50.0,          # 矩形：水面宽 = 底宽
        "q0": 180.0,
        "h_target": 7.00,
        "a0": 0.50,
        "kg": 360.0,
        "n": 0.014,
        "duration_s": 86400.0,          # 24 h
        "t_step_s": 900.0,
        "band_normal_m": 0.15,
        "band_alarm_m": 0.30,
        "settle_tol_m": 0.020,
        "settle_window_s": 1800.0,
        "step_dq": 45.0,
    }


def make_lit() -> dict:
    """文献参数化：唐河节制闸→放水河节制闸渠段（公开值取区间中值）。"""
    b, m, h = 20.75, 2.25, 4.50
    return {
        "tag": "A1/B1-B3",
        "label": "唐河—放水河渠段（文献参数化）",
        "pool_len_m": 25567.0,          # 桩号 1046+216 ~ 1071+783
        "bottom_m": b,                  # 18.5~23.0 中值
        "m_slope": m,                   # 2.0~2.5 中值
        "water_width_m": b + 2.0 * m * h,   # = 41.00 m（设计水深处）
        "q0": 135.0,                    # 设计流量
        "h_target": h,                  # 设计水深
        "a0": 0.50,
        "kg": 135.0 / 0.50,             # 重新标定：Q = kg·a·√(h/ht)，a0 时恰为 q0
        "n": 0.0150,                    # 总干渠设计糙率
        "duration_s": 86400.0,          # 24 h：水面面积 3.49× ⇒ 调整时间尺度约 3.5×，
                                        # 按原案例 4 h 收敛推算需 ~14 h，取 24 h 留裕度
        "t_step_s": 900.0,
        "band_normal_m": 0.15,          # ★ 与中线推荐值 ±0.15 m 一致（非自定）
        "band_alarm_m": 0.30,
        "settle_tol_m": 0.020,
        "settle_window_s": 1800.0,
        "step_dq": 45.0,
    }


# ======================================================================
# 自洽性校验
# ======================================================================
def selfcheck(c: dict) -> dict:
    b, m, h, Q, n = c["bottom_m"], c["m_slope"], c["h_target"], c["q0"], c["n"]
    A = h * (b + m * h)
    B = b + 2.0 * m * h
    P = b + 2.0 * h * math.sqrt(1.0 + m * m)
    R = A / P
    S0_section = (n * Q / (A * R ** (2.0 / 3.0))) ** 2
    S0_rectonly = (n * Q / (b * h * (b * h / (b + 2.0 * h)) ** (2.0 / 3.0))) ** 2
    c_ = math.sqrt(G * h)
    L = c["pool_len_m"]
    return {
        "断面": {"底宽_m": b, "边坡_m": m, "过水断面积_m2": round(A, 2),
                 "水面宽_m": round(B, 2), "湿周_m": round(P, 2), "水力半径_m": round(R, 3)},
        "断面流速_m_s": round(Q / A, 3),
        "公开流速口径_m_s": "约 1 m/s",
        "反解纵坡_梯形": round(S0_section, 8),
        "反解纵坡_梯形_倒数": round(1.0 / S0_section),
        "反解纵坡_若误用矩形式": round(S0_rectonly, 8),
        "矩形式偏大倍数": round(S0_rectonly / S0_section, 3),
        "公开纵坡区间": "1/25000 ~ 1/20000 (4.0e-5 ~ 5.0e-5)",
        "纵坡是否落在公开区间": bool(4.0e-5 <= S0_section <= 5.0e-5),
        "波速_m_s": round(c_, 2),
        "L_over_c_min": round(L / c_ / 60, 1),
        "2L_over_c_min": round(2 * L / c_ / 60, 1),
        "水面面积_m2": round(L * B),
    }


# ======================================================================
# 0D
# ======================================================================
def run0d(c: dict) -> dict:
    p = M0.Params(
        pool_len_m=c["pool_len_m"],
        pool_width_m=c["water_width_m"],
        q0=c["q0"], h_target=c["h_target"], a0=c["a0"], kg=c["kg"],
        band_normal_m=c["band_normal_m"], band_alarm_m=c["band_alarm_m"],
        settle_tol_m=c["settle_tol_m"], settle_window_s=c["settle_window_s"],
        dt=5.0, duration_s=c["duration_s"], t_step_s=c["t_step_s"],
        step_dq=c["step_dq"],
    )
    sc = M0.Scenario(key="lit", label=c["label"], pools=1,
                     step_dq=c["step_dq"], feedforward=False, note="")
    r = M0.simulate(p, sc)
    j = M0.judge(p, sc, r)
    cm = r["metrics"]["pool1"]
    return {
        "model": "0d-lumped",
        "water_width_m": c["water_width_m"],
        "area_m2": round(p.area_m2),
        "rise_rate_mm_s": round(c["step_dq"] / p.area_m2 * 1000, 4),
        "peak_mm": round(cm["peak_dev_m"] * 1000, 2),
        "over_normal_band_s": cm["time_over_normal_band_s"],
        "over_alarm_band_s": cm["time_over_alarm_band_s"],
        "final_window_mm": round(cm["final_window_max_dev_m"] * 1000, 3),
        "gate_sat_s": round(max(x.sat_time for x in r["pools"]), 1),
        "gate_up_capacity_m3s": round(p.gate_up_capacity, 1),
        "code": j["code"],
        "rules_pass": {x["rule"]: bool(x["pass"]) for x in j["rules"]},
    }


# ======================================================================
# 1D
# ======================================================================
def run1d(c: dict, ncell: int, *, m_slope: float = 0.0,
          step_dq: float = None, duration_s: float = None) -> dict:
    m = m_slope
    width = c["bottom_m"] if m > 0 else c["water_width_m"]
    sd = c["step_dq"] if step_dq is None else step_dq
    dur = c["duration_s"] if duration_s is None else duration_s
    pd = {
        "pool_width_m": width, "pool_len_m": c["pool_len_m"],
        "h_target": c["h_target"], "q0": c["q0"], "kg": c["kg"],
        "duration_s": dur, "t_step_s": c["t_step_s"],
        "settle_window_s": c["settle_window_s"], "band_normal_m": c["band_normal_m"],
        "a0": c["a0"], "kp": 1.60, "ki": 1.60e-3,
        "ctrl_period_s": 300.0, "deadband_m": 0.010,
        "meas_resolution_m": 0.001, "rate_limit_per_s": 0.02 / 60.0,
    }
    s0 = None
    if m > 0:
        sec = SVE2.Section(width, m)
        s0 = SVE2.section_normal_slope(sec, c["h_target"], c["q0"], c["n"])
    t0 = time.time()
    r = SVE2.sve_1d_v2(pd, ncell, sd, m_slope=m, n=c["n"], s0=s0)
    r["wall_s"] = round(time.time() - t0, 1)
    r["step_dq"] = sd
    r["duration_s"] = dur
    r["section"] = {"bottom_m": width, "m_slope": m,
                    "water_width_target_m": round(width + 2.0 * m * c["h_target"], 2)}
    return r


def pair_extrapolate(a: dict, b: dict) -> float:
    """一阶 Richardson 外推：粗网格为 a、细网格为 b（b 的格数约为 a 的 2 倍）。"""
    return 2.0 * b["peak_mean_mm"] - a["peak_mean_mm"]


# ======================================================================
# 主流程
# ======================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01R 文献参数化重跑")
    ap.add_argument("--quick", action="store_true", help="只跑单级网格")
    ap.add_argument("--out", default=HERE)
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out)
    raw_dir = os.path.join(out_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    orig, lit = make_orig(), make_lit()
    chk = selfcheck(lit)

    print("=" * 78)
    print(f"{CASE_ID}　{CASE_TITLE}")
    print(f"对象：{OBJECT}")
    print(f"活动号 {LOCAL_CAMPAIGN_ID}｜证据等级 {EVIDENCE_GRADE}｜sealed: false")
    print("=" * 78)

    print("\n【一】参数自洽性校验（用公开值交叉验证，不是自证）")
    s = chk["断面"]
    print(f"  断面：底宽 {s['底宽_m']} m，边坡 {s['边坡_m']}，设计水深 {lit['h_target']} m")
    print(f"        过水断面积 {s['过水断面积_m2']} m²，水面宽 {s['水面宽_m']} m，"
          f"湿周 {s['湿周_m']} m，水力半径 {s['水力半径_m']} m")
    print(f"  ① 断面流速 = {chk['断面流速_m_s']} m/s  ← 公开口径 {chk['公开流速口径_m_s']}")
    print(f"  ② 反解纵坡 = {chk['反解纵坡_梯形']:.4e} = 1/{chk['反解纵坡_梯形_倒数']:,}"
          f"  ← 公开区间 {chk['公开纵坡区间']}")
    print(f"     落在公开区间内：{'是 ✓' if chk['纵坡是否落在公开区间'] else '否 ✗'}")
    print(f"  ③ 波速 c = {chk['波速_m_s']} m/s，L/c = {chk['L_over_c_min']} min，"
          f"2L/c = {chk['2L_over_c_min']} min")
    print(f"  附：若误用矩形反解式，纵坡会偏大 {chk['矩形式偏大倍数']} 倍 —— 这是本次"
          f"给一维求解器补 s0 入口的原因")

    print(f"\n【二】尺度对比")
    a_o = orig["pool_len_m"] * orig["water_width_m"]
    a_l = lit["pool_len_m"] * lit["water_width_m"]
    print(f"  水面面积：原 {a_o:,.0f} m² → 真实 {a_l:,.0f} m²（{a_l/a_o:.2f}×）")
    print(f"  同样 +{lit['step_dq']:.0f} m³/s 引起的水位上升速率："
          f"{orig['step_dq']/a_o*1000:.3f} mm/s → {lit['step_dq']/a_l*1000:.3f} mm/s")

    results = {"meta": {
        "case_id": CASE_ID, "title": CASE_TITLE, "object": OBJECT,
        "local_campaign_id": LOCAL_CAMPAIGN_ID, "evidence_grade": EVIDENCE_GRADE,
        "sealed": False, "source_note": SOURCE_NOTE,
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "controller_note": "控制器增益 kp/ki/控制周期/速率限**未重新整定**，"
                           "与原案例相同——差异全部归因于对象参数变化。",
    }, "params": {"orig": orig, "literature": lit},
        "selfcheck": chk, "runs": {}}

    # ---------------- 0D ----------------
    print("\n【三】0D 集总模型（原案例用的就是它）")
    for code, cfg in (("A0", orig), ("A1", lit)):
        r = run0d(cfg)
        results["runs"][code] = r
        print(f"  {code}  {cfg['label']}")
        print(f"      水面宽 {r['water_width_m']:.1f} m，A = {r['area_m2']:,} m²，"
              f"上升速率 {r['rise_rate_mm_s']:.4f} mm/s")
        print(f"      池均峰值 {r['peak_mm']:7.2f} mm｜越正常条带 {r['over_normal_band_s']:7.1f} s"
              f"｜越告警条带 {r['over_alarm_band_s']:7.1f} s")
        print(f"      末段 30 min 最大偏差 {r['final_window_mm']:7.3f} mm｜"
              f"闸门饱和 {r['gate_sat_s']:.1f} s｜向上余量 {r['gate_up_capacity_m3s']:.0f} m³/s")
        print(f"      判定 = {r['code']}　规则 {r['rules_pass']}")

    # ---------------- 1D ----------------
    print("\n【四】1D 圣维南（析因：尺度效应 vs 断面效应）")
    grids = [120, 240]
    lit_grids = [128, 256] if not args.quick else [256]

    plan = [
        ("B0", orig, dict(m_slope=0.0), "原参数 · 矩形（与 v1 基准对齐）"),
        ("B1", lit, dict(m_slope=0.0), "真实尺度 · 矩形（只动尺度）"),
        ("B2", lit, dict(m_slope=2.25), "真实几何 · 梯形 m=2.25（真实断面）"),
        ("B3", lit, dict(m_slope=2.25, step_dq=33.75), "真实几何 · 梯形 · 同比例 +25%"),
        ("B4", lit, dict(m_slope=2.25, step_dq=25.0),
         "真实几何 · 梯形 · 只到加大流量 160 m³/s（工程允许上限）"),
    ]
    for code, cfg, kw, desc in plan:
        gs = grids if code == "B0" else lit_grids
        dur = 14400.0 if code == "B0" else None
        got = {}
        for nc in gs:
            r = run1d(cfg, nc, duration_s=dur, **kw)
            got[nc] = r
            print(f"  {code} [{desc}]  {nc:>4} 格  dx={r['dx_m']:.1f} m  dt={r['dt_s']:>5} s  "
                  f"({r['steps']} 步, {r['wall_s']} s)")
            print(f"        池均峰值 {r['peak_mean_mm']:7.1f} mm｜闸前峰值 {r['peak_gate_mm']:7.1f} mm"
                  f"｜末段 {r['final_mean_mm']:6.2f} mm｜越带 {r['over_band_s']:6.1f} s"
                  f"｜质量守恒 {r['mass_balance_rel_err']:.1e}")
        ext = pair_extrapolate(got[gs[0]], got[gs[-1]]) if len(gs) > 1 else got[gs[0]]["peak_mean_mm"]
        print(f"        → Richardson 外推池均峰值 ≈ {ext:.1f} mm")
        rec = dict(got[gs[-1]])
        s = rec.pop("series", None)
        if s:
            k = max(1, len(s["t"]) // 400)
            rec["series_sampled"] = {
                "t_s": [round(x, 1) for x in s["t"][::k]],
                "dev_mean_mm": [round(x * 1000, 2) for x in s["dev_mean"][::k]],
                "dev_gate_mm": [round(x * 1000, 2) for x in s["dev_gate"][::k]],
                "every_s": round(k * got[gs[-1]]["dt_s"], 1),
            }
        rec["grids"] = gs
        rec["extrapolated_peak_mean_mm"] = round(ext, 1)
        rec["desc"] = desc
        results["runs"][code] = rec

    # ---------------- 汇总 ----------------
    print("\n【五】汇总：真实参数把判定带到了哪")
    print(f"  {'工况':<6}{'模型':<14}{'池均峰值':>10}{'外推':>9}{'越带/s':>9}{'判定':>6}")
    rows = []
    for code in ("A0", "A1", "B0", "B1", "B2", "B3", "B4"):
        r = results["runs"][code]
        peak = r.get("peak_mm", r.get("peak_mean_mm"))
        ext = r.get("extrapolated_peak_mean_mm", peak)
        over = r.get("over_normal_band_s", r.get("over_band_s"))
        model = "0D 集总" if "peak_mm" in r else ("1D 矩形" if r.get("section", {}).get("m_slope", 0) == 0 else "1D 梯形")
        code_j = r.get("code") or ("A1" if ext <= 150.0 else "A3")
        print(f"  {code:<6}{model:<14}{peak:>10.2f}{ext:>9.1f}{over:>9.1f}{str(code_j):>6}")
        rows.append({"id": code, "model": model, "peak_mm": peak,
                     "extrapolated_mm": ext, "over_band_s": over, "code": code_j,
                     "desc": r.get("desc", "")})
    results["summary"] = rows

    out_json = os.path.join(raw_dir, "bm01r_literature.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
