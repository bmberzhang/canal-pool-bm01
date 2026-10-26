# -*- coding: utf-8 -*-
"""
(W9) 网格无关性检查 —— 必须先于 W3'' 执行
================================================================================
三模型第三轮 B 仲裁指出的依赖疏漏：W0-c 的"救回"结论（控制周期 300 s→60 s）
是在 **dx = 25 m** 单一网格上得到的；若峰值对网格敏感（>15%），
则 W0-c 的结论与 W3'' 的设计前提都要重做。

本脚本先做 W9：在 **dx ∈ {50, 25, 12.5} m** 三档网格上重跑关键算例，
检验峰值与判定码是否网格稳定。

算例：
  A) L=2000, ΔQ=45, ctrl_period=300 s（那片"三档全 A3"的代表格）
  B) L=2000, ΔQ=45, ctrl_period=60 s（W0-c 的"救回"格）
  C) L=6000, ΔQ=45, ctrl_period=300 s（翻转格，论文主案例）

判据（预先写死）：若同一算例在三档网格上的全池平均峰值**最大相对偏差 ≤15%**
且**判定码一致** → W9 通过，W0-c 与 W3'' 可继续；
否则 W9 不通过，必须先解决网格敏感性问题，且 W0-c 的"救回"结论降格。

输出：raw/bm01_w9_grid_check.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BAND_MM = 150.0
B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
DXS = (50.0, 25.0, 12.5)
CASES = ((2000.0, 45.0, 300.0, "A_短池_原周期"),
         (2000.0, 45.0, 60.0, "B_短池_快周期"),
         (6000.0, 45.0, 300.0, "C_主案例"))
ARMS = (("point", dict(ctrl_on="gate")),
        ("mean", dict(ctrl_on="mean")))


def base_p(L, period):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": period, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def main():
    t0 = time.time()
    out = {"purpose": "W9 网格无关性检查（W3'' 的前置）", "tolerance_pct": 15.0,
           "cases": {}, "wall_total_s": None}

    for L, dq, period, tag in CASES:
        rows = []
        for dx in DXS:
            ncell = int(round(L / dx))
            arms = {}
            for name, kw in ARMS:
                r = V2.sve_1d_v2(base_p(L, period), ncell, dq, m_slope=0.0, **kw)
                arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                              "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3",
                              "mass_err": float("%.2e" % r["mass_balance_rel_err"])}
            rows.append({"dx_m": dx, "ncell": ncell, "arms": arms})
            print("[W9] %-16s L=%5.0f ΔQ=%4.1f period=%5.0fs dx=%5.1f ncell=%4d | "
                  "point=%8.1f %s | mean=%8.1f %s"
                  % (tag, L, dq, period, dx, ncell, arms["point"]["peak_mean_mm"],
                     arms["point"]["judgement"], arms["mean"]["peak_mean_mm"],
                     arms["mean"]["judgement"]), flush=True)
        # 网格稳定性：以最细网格为基准
        ref = rows[-1]["arms"]["mean"]["peak_mean_mm"]
        devs = [{"dx_m": r["dx_m"],
                 "mean_peak_mm": r["arms"]["mean"]["peak_mean_mm"],
                 "rel_dev_pct": round(abs(r["arms"]["mean"]["peak_mean_mm"] - ref) / ref * 100, 2),
                 "code": r["arms"]["mean"]["judgement"]} for r in rows]
        codes = set(d["code"] for d in devs)
        maxdev = max(d["rel_dev_pct"] for d in devs)
        out["cases"][tag] = {"pool_len_m": L, "dq": dq, "ctrl_period_s": period,
                             "rows": rows, "deviations": devs,
                             "max_rel_dev_pct": maxdev,
                             "codes_consistent": len(codes) == 1,
                             "pass": bool(maxdev <= 15.0 and len(codes) == 1)}
        print("\n>>> %s：最大相对偏差 %.2f%%，判定码 %s → %s\n"
              % (tag, maxdev, sorted(codes),
                 "通过" if out["cases"][tag]["pass"] else "**不通过**"), flush=True)

    all_pass = all(v["pass"] for v in out["cases"].values())
    worst = max(v["max_rel_dev_pct"] for v in out["cases"].values())
    out["wall_total_s"] = round(time.time() - t0, 1)
    out["W9_pass"] = all_pass
    out["verdict"] = (
        "W9 **通过**：三个关键算例在三档网格（dx=50/25/12.5 m）上最大相对偏差 %.2f%%（阈值 15%%）、"
        "判定码三档一致 ⇒ W0-c 的『救回』结论与 W3'' 的设计前提成立，可继续。"
        % worst if all_pass else
        "W9 **不通过**：最大相对偏差 %.2f%%（阈值 15%%）或判定码不一致 ⇒ "
        "必须先解决网格敏感性；W0-c 的『救回』结论与 W3'' 设计须重做。" % worst)

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_w9_grid_check.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
