# -*- coding: utf-8 -*-
"""
(W3''-fix) 修正 W3'' 第 2 段的取值错误，在**起爆点之上**重做周期对照
================================================================================
W3'' 主脚本用"起爆区间中点" (70+75)/2 = 72.5 做周期对照，而 72.5 位于起爆点**之下**
（该处三档本来就不同时判 A3）⇒ 判据 `p300 全A3 且 p60 不全A3` 必然为假，
输出"未通过"。这是取值 bug，不是结论。

本脚本在 **ΔQ = 76 / 80**（均在 L=6000 的起爆区内，见 W3'' 第 1、3 段）上，
扫 ctrl_period ∈ {60, 150, 300} s，检验"缩短周期能否救回"。

同时补一条关键对照：把 L=2000（失稳型起爆）与 L=6000（边际型起爆）的
**起爆形态**并排写出——前者 49.1 → 3820.7 mm（78 倍，分岔），
后者 137.7 → 152.2 mm（1.10 倍，边际越带）。

输出：raw/bm01_w3pp_fix.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
L = 6000.0
NCELL = 240
DQS = [76.0, 80.0]
PERIODS = [60.0, 150.0, 300.0]
ARMS = (("point", dict(ctrl_on="gate")),
        ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
        ("mean", dict(ctrl_on="mean")))


def base_p(period):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": period, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def main():
    t0 = time.time()
    out = {"purpose": "W3''-fix：在起爆点之上重做周期对照（修正取值 bug）",
           "config": {"L_m": L, "ncell": NCELL, "dqs": DQS, "periods_s": PERIODS},
           "rows": [], "wall_total_s": None}
    for dq in DQS:
        for period in PERIODS:
            arms = {}
            for name, kw in ARMS:
                r = V2.sve_1d_v2(base_p(period), NCELL, dq, m_slope=0.0, **kw)
                arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                              "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3"}
            n_a3 = sum(1 for a in ARMS if arms[a[0]]["judgement"] == "A3")
            row = {"dq": dq, "ctrl_period_s": period, "arms": arms,
                   "n_A3_arms": n_a3, "ALL_A3": n_a3 == 3}
            out["rows"].append(row)
            print("[FIX] L=6000 ΔQ=%5.1f period=%5.0fs | mean=%8.1f %s | point=%8.1f %s | ALL_A3=%s"
                  % (dq, period, arms["mean"]["peak_mean_mm"], arms["mean"]["judgement"],
                     arms["point"]["peak_mean_mm"], arms["point"]["judgement"], row["ALL_A3"]),
                  flush=True)

    rescued = []
    for dq in DQS:
        r300 = next(r for r in out["rows"] if r["dq"] == dq and r["ctrl_period_s"] == 300.0)
        r60 = next(r for r in out["rows"] if r["dq"] == dq and r["ctrl_period_s"] == 60.0)
        r150 = next(r for r in out["rows"] if r["dq"] == dq and r["ctrl_period_s"] == 150.0)
        if r300["ALL_A3"] and not r60["ALL_A3"]:
            rescued.append({"dq": dq, "p300_all_A3": True, "p60_all_A3": False,
                            "p150_all_A3": r150["ALL_A3"],
                            "mean_p300_mm": r300["arms"]["mean"]["peak_mean_mm"],
                            "mean_p60_mm": r60["arms"]["mean"]["peak_mean_mm"],
                            "mean_p150_mm": r150["arms"]["mean"]["peak_mean_mm"]})
    out["rescued"] = rescued
    out["wall_total_s"] = round(time.time() - t0, 1)

    out["verdict"] = (
        ("W3'' **第 2 段通过（修正后）**：在 L=6000 的起爆区内（ΔQ=%s），"
         "控制周期 300 s→60 s 使『三档全 A3』被救回："
         + "；".join("ΔQ=%.0f：全池平均 %.1f → %.1f mm" % (x["dq"], x["mean_p300_mm"], x["mean_p60_mm"])
                     for x in rescued)
         + "。150 s 的效果："
         + "、".join("ΔQ=%.0f %s" % (x["dq"], "也已救回" if not x["p150_all_A3"] else "仍全 A3")
                     for x in rescued)
         + "。⇒ T_ctrl 机制在**两个池长上复现**，可写成可证伪通则。")
        if rescued else
        "W3'' **第 2 段仍未通过**：在 ΔQ=76/80 上缩短周期未改变判定 ⇒ 须复核 W0-c 的可复现性。")

    # 起爆形态并排（不重新计算，直接引用两个脚本的实测值）
    out["onset_signature_contrast"] = {
        "L2000_bifurcation": {"dq_prev": 14, "mean_prev_mm": 49.1,
                              "dq_next": 15, "mean_next_mm": 3820.7,
                              "ratio": 77.8, "type": "失稳型（准则性分岔）",
                              "onset_m3s": [14, 15]},
        "L6000_marginal": {"dq_prev": 70, "mean_prev_mm": 137.7,
                           "dq_next": 75, "mean_next_mm": 152.2,
                           "ratio": 1.10, "type": "边际型（池均稳态响应越带）",
                           "onset_m3s": [70, 75]},
        "note": "两者形态相差 71 倍（77.8 / 1.10）⇒ 是**两个不同机制**，不可混为一谈："
                "短池是控制采样过慢导致的分岔，长池是池均响应平稳越带。"
                "C3 的闭式阈值（ΔQ*=73.5）在长池上吻合到 2%（实测 72~75），"
                "在短池上高估 2.9 倍（实测 14.5）——**C3 只适用于边际型，不适用于失稳型**。"}

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_w3pp_fix.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("\n起爆形态对照:", json.dumps(out["onset_signature_contrast"], ensure_ascii=False)[:400])
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
