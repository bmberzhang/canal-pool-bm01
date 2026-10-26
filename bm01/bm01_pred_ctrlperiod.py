# -*- coding: utf-8 -*-
"""
(W0-c) 因一的第二个子机制：控制周期 vs 池体波行时间
================================================================================
W0-b 实测：L=2000 在 ΔQ≈14.5 起爆（三档全 A3）；而 L=6000 到 ΔQ=50、L=12000 到 40 都不起爆。
解析式 C3（只含执行器能力 k_g·r）在 L=2000 上高估 2.9 倍 ⇒ C3 漏掉了真正的机制。

候选解释：L=2000 的池体波行时间 L/c ≈ 227 s（c≈8.8 m/s），而**控制周期是 300 s**
—— 控制器采样比池体能响应还慢 ⇒ 无论执行器多快，环路都守不住。这是执行器能力之外的
**第二个"因一"子机制：控制采样周期 vs 池体动力学**。

可证伪设计：固定 L=2000 与 ΔQ，只改控制周期 ctrl_period_s ∈ {60, 300}：
  - 若把周期从 300 s 缩到 60 s 后起爆点显著上移（或不复现）⇒ 子机制 (1b) 成立，
    "因一"必须写成两项（执行器能力 + 控制采样相对池体动力学）；
  - 若起爆点不动 ⇒ (1b) 不成立，L=2000 的失败另有原因。

设置：L=2000（ncell=80, dx=25 m），12 h，矩形，ΔQ ∈ {14, 16, 20, 30, 45}，三档取样。
输出：raw/bm01_pred_ctrlperiod.json
"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

L = 2000.0
NCELL = int(round(L / 25.0))
BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
DQS = [14.0, 16.0, 20.0, 30.0, 45.0]
PERIODS = [60.0, 300.0]
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
    out = {"purpose": "W0-c 控制周期 vs 池体动力学：因一是否含第二个子机制",
           "config": {"L_m": L, "ncell": NCELL, "duration_h": 12.0, "dqs": DQS,
                      "periods_s": PERIODS},
           "results": {}, "wall_total_s": None}
    for period in PERIODS:
        rows = []
        for dq in DQS:
            arms = {}
            for name, kw in ARMS:
                r = V2.sve_1d_v2(base_p(period), NCELL, dq, m_slope=0.0, **kw)
                arms[name] = round(r["peak_mean_mm"], 2)
            n_a3 = sum(1 for a in ARMS if arms[a[0]] > BAND_MM)
            rows.append({"dq": dq, "peaks_mm": arms, "n_A3_arms": n_a3, "ALL_A3": n_a3 == 3})
            print("[CTRL] period=%5.0fs ΔQ=%4.1f | point=%8.1f avg100=%8.1f mean=%8.1f | A3=%d %s"
                  % (period, dq, arms["point"], arms["avg100"], arms["mean"], n_a3,
                     "ALL_A3" if n_a3 == 3 else "     "), flush=True)
        a3 = [r["dq"] for r in rows if r["ALL_A3"]]
        ok = [r["dq"] for r in rows if not r["ALL_A3"]]
        out["results"]["period%d" % int(period)] = {
            "rows": rows, "all_A3_dqs": a3,
            "onset_bracket_m3s": [max(ok), min(a3)] if (ok and a3) else None,
            "mean_arm_max_mm": max(r["peaks_mm"]["mean"] for r in rows)}
        print("\n>>> ctrl_period=%.0f s：全 A3 的 ΔQ=%s；未全 A3 的 ΔQ=%s；"
              "全池平均档最大峰值 %.1f mm\n"
              % (period, a3, ok, out["results"]["period%d" % int(period)]["mean_arm_max_mm"]),
              flush=True)

    p60 = out["results"]["period60"]
    p300 = out["results"]["period300"]
    # 判据：周期缩短后，原本全 A3 的 ΔQ 是否被救回
    rescued = [d for d in p300["all_A3_dqs"] if d not in p60["all_A3_dqs"]]
    out["wall_total_s"] = round(time.time() - t0, 1)
    if rescued:
        out["verdict"] = (
            "**因一存在第二个子机制：控制采样周期 vs 池体波行时间。**"
            "L=2000、ctrl_period=300 s（池体波行时间 L/c≈227 s < 周期）时，ΔQ≥%.1f 三档全失败；"
            "把周期缩到 60 s 后，ΔQ=%s 被救回（不再全 A3），全池平均档最大峰值由 %.1f mm 降到 %.1f mm。"
            "⇒ 短池的失效**不是执行器能力（k_g·r）问题**（执行器参数一字未改），"
            "而是「控制器采样比池体动力学还慢」。"
            "因此因一应写成两项：① 执行器能力；② 控制采样周期相对池体波行时间。"
            % (min(p300["all_A3_dqs"]) if p300["all_A3_dqs"] else float('nan'), rescued,
               p300["mean_arm_max_mm"], p60["mean_arm_max_mm"]))
    else:
        out["verdict"] = ("缩短控制周期**未能**改变起爆点 ⇒ 子机制 (1b) 不成立，"
                          "L=2000 的失败另有原因，因一暂不增加第二项。")

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_pred_ctrlperiod.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
