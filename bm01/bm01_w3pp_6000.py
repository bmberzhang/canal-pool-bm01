# -*- coding: utf-8 -*-
"""
(W3'') 把"控制周期 vs 池体波行时间"从单案例升级为可证伪通则
================================================================================
三模型第三轮 B 指定，且**要求预先写死失败分支**：
  - 若 L=6000 在 ΔQ ≤ 80 m³/s 内**找不到**"三档全 A3"起爆点 → 
    判据降格为「经验性设计提示」，**不得写入摘要**；
  - 若能找到 → 起爆点作为第二证据点，检验 T_ctrl 判据是否可写成可证伪通则。

三段：
  第 1 段：L=6000（ncell=240, dx=25 m）扫 ΔQ ∈ {55,60,65,70,75,80}，
          找"三档全 A3"起爆区间（≤50 已由 W0-b 证无起爆）。
  第 2 段：固定 ΔQ（对第 1 段找到的起爆点及其下邻），
          扫 ctrl_period ∈ {60,150,300} s，检验"缩短周期能否救回"是否复现。
  第 3 段：起爆点附近 ΔQ ± 1, ± 2 扰动，检验是否为刀锋（不可复现的跳跃）。
  前置：W9 已通过（三档网格最大偏差 5.48%、判定码一致）。

输出：raw/bm01_w3pp_6000.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

L = 6000.0
NCELL = 240
BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
DQS_SCAN = [55.0, 60.0, 65.0, 70.0, 75.0, 80.0]
PERIODS = [60.0, 150.0, 300.0]
ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def base_p(period=300.0):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": period, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def run(dq, period=300.0):
    arms = {}
    for name, kw in ARMS:
        r = V2.sve_1d_v2(base_p(period), NCELL, dq, m_slope=0.0, **kw)
        arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                      "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3"}
    n_a3 = sum(1 for a in ARMS if arms[a[0]]["judgement"] == "A3")
    return {"dq": dq, "ctrl_period_s": period, "arms": arms, "n_A3_arms": n_a3,
            "ALL_A3": n_a3 == 3}


def main():
    t0 = time.time()
    out = {"purpose": "W3''：把 T_ctrl 机制从单案例升级为可证伪通则（失败分支已预先写死）",
           "config": {"L_m": L, "ncell": NCELL, "dx_m": L / NCELL, "band_mm": BAND_MM,
                      "dqs_scan": DQS_SCAN, "periods_s": PERIODS,
                      "W9_precondition_passed": True},
           "part1_onset_scan": [], "part2_period_compare": [], "part3_boundary_perturb": [],
           "wall_total_s": None}

    print("=== 第 1 段：L=6000 起爆扫描（ΔQ ≤ 80）===", flush=True)
    for dq in DQS_SCAN:
        row = run(dq)
        out["part1_onset_scan"].append(row)
        print("[W3''] ΔQ=%4.1f | point=%8.1f %s | avg100=%8.1f %s | mean=%8.1f %s | ALL_A3=%s"
              % (dq, row["arms"]["point"]["peak_mean_mm"], row["arms"]["point"]["judgement"],
                 row["arms"]["avg100"]["peak_mean_mm"], row["arms"]["avg100"]["judgement"],
                 row["arms"]["mean"]["peak_mean_mm"], row["arms"]["mean"]["judgement"],
                 row["ALL_A3"]), flush=True)

    a3 = [r["dq"] for r in out["part1_onset_scan"] if r["ALL_A3"]]
    ok = [r["dq"] for r in out["part1_onset_scan"] if not r["ALL_A3"]]
    if a3 and ok:
        bracket = [max(ok), min(a3)]
        focus_dq = sum(bracket) / 2.0
    elif a3:
        bracket = [None, min(a3)]
        focus_dq = min(a3)
    else:
        bracket = None
        focus_dq = 45.0        # 无起爆时，退而检验周期对照是否仍成立
    out["part1_onset_bracket_m3s"] = bracket
    out["part1_has_onset"] = bool(a3)
    print("\n>>> L=6000 起爆区间 = %s（起爆 ΔQ 列表 %s）\n" % (bracket, a3), flush=True)

    print("=== 第 2 段：周期对照（ΔQ=%.1f）===" % focus_dq, flush=True)
    p60 = p150 = p300 = None
    for period in PERIODS:
        row = run(float(focus_dq), period)
        out["part2_period_compare"].append(row)
        if period == 60.0: p60 = row
        if period == 150.0: p150 = row
        if period == 300.0: p300 = row
        print("[W3''] period=%5.0fs ΔQ=%4.1f | mean=%8.1f %s | point=%8.1f %s | ALL_A3=%s"
              % (period, focus_dq, row["arms"]["mean"]["peak_mean_mm"],
                 row["arms"]["mean"]["judgement"], row["arms"]["point"]["peak_mean_mm"],
                 row["arms"]["point"]["judgement"], row["ALL_A3"]), flush=True)

    print("\n=== 第 3 段：起爆点附近 ΔQ 扰动 ===", flush=True)
    if bracket and bracket[0] is not None:
        for dq in [bracket[0] - 2, bracket[0] - 1, bracket[0], bracket[1], bracket[1] + 1]:
            row = run(float(dq))
            out["part3_boundary_perturb"].append(row)
            print("[W3''] ΔQ=%5.1f | mean=%8.1f %s | ALL_A3=%s"
                  % (dq, row["arms"]["mean"]["peak_mean_mm"],
                     row["arms"]["mean"]["judgement"], row["ALL_A3"]), flush=True)

    # ---- 预先写死的判定树 ----
    rescued = (p300 and p60 and (not p60["ALL_A3"]) and p300["ALL_A3"])
    p300_mean = p300["arms"]["mean"]["peak_mean_mm"] if p300 else None
    p60_mean = p60["arms"]["mean"]["peak_mean_mm"] if p60 else None

    if out["part1_has_onset"] and rescued:
        out["grade"] = "可证伪通则（tier-1）"
        out["verdict"] = (
            "W3'' **{0}**：L=6000 存在起爆区间 {1}（ΔQ ≤ 80）；在同一 ΔQ={2:.1f} 上，"
            "控制周期 300 s→60 s 使『三档全 A3』被救回（全池平均峰值 {3:.1f} → {4:.1f} mm）"
            "⇒ T_ctrl 机制在两个池长上复现，可写成**可证伪通则**："
            "『当控制周期不小于池体波行时间时，判定码与取样位置无关地整体失败；缩短周期可消除』。"
            .format("通过", bracket, focus_dq, p300_mean or -1, p60_mean or -1))
    elif rescued and not out["part1_has_onset"]:
        out["grade"] = "经验性设计提示（tier-2）"
        out["verdict"] = (
            "W3'' **部分通过**：L=6000 在 ΔQ ≤ 80 内**无**『三档全 A3』起爆点，"
            "但周期对照在 ΔQ={0:.1f} 上仍显示 300 s→60 s 使全池平均峰值 {1:.1f} → {2:.1f} mm。"
            "按预先写死的失败分支，判据**降格为经验性设计提示**，"
            "**不得写入摘要**，只能出现在讨论与局限中。".format(focus_dq, p300_mean or -1, p60_mean or -1))
    else:
        out["grade"] = "未通过（tier-3）"
        out["verdict"] = ("W3'' **未通过**：既未在 L=6000 找到起爆点，周期对照也未显示救回效果 "
                          "⇒ 判据降格为经验性设计提示，不得写入摘要；且须复核 W0-c 的可复现性。")
    out["wall_total_s"] = round(time.time() - t0, 1)

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_w3pp_6000.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("档位:", out["grade"], "｜已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
