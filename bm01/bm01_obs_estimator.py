# -*- coding: utf-8 -*-
"""
(E-OBS) 可观测性实验：理想"全池平均"能不能被工程上真能装的量测替代？
================================================================================
这是三模型 Q&A 列出的第二硬伤（R2）：
  "全池蓄量等效水位是**理想量测**。真实渠池没有'全池平均水位计'。
   你的 A1/A3 翻转建立在一个不可实现的量测上 → 工程价值存疑。"

可证伪设计：判定码**始终以真实全池平均峰值**为准（判据不变），
只把**控制器实际喂进去的反馈量**换成四种工程可实现的估计量：

  midpoint : 池中点单点水位计（1 支）
  mp3      : 池内 3 点等距平均（0.25L / 0.5L / 0.75L，3 支）
  mp5      : 池内 5 点等距平均（5 支）
  volbal   : 容积平衡估计 ĥ = h0 + ∫(q_in − q_gate)dt /(B0·L)
             ——**只需上下游流量计，完全不需要水位计**，
               而且与 0D 积分蓄量模型同构（呼应"0D ≡ 1D 全池平均"主线）

判据：
  - 若有任一可实现估计量能把判定码还原成 A1（与理想全池平均一致）
    → R2 被闭合：翻转**不需要理想量测**，工程上可落地；
  - 若全部失败 → 必须诚实下调结论（翻转依赖理想量测，工程意义有限）。

同时报 est_track_err_mm：估计量偏离真实池均的最大值，给出"要多少支水位计才够"的定量答案。

设置：BM-01，L=6000 m，240 格（dx=25 m），12 h，矩形，+45 m³/s，kp=1.6，ki=0.0016
      （与 E-MC / E-TUNE 完全同构，可比）。
输出：raw/bm01_obs_estimator.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

L = 6000.0
NCELL = 240
DQ = 45.0
BAND_MM = 150.0


def base_p():
    return {"pool_len_m": L, "pool_width_m": 50.0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": 360.0, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


# 参照档（复现 E-MC 的 L=6000/ΔQ=45 行，做回归自检）
REF_ARMS = (("point", dict(ctrl_on="gate")),
            ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
            ("mean", dict(ctrl_on="mean")))
REF_EXPECT = {"point": 470.8, "avg100": 658.8, "mean": 75.0}   # E-MC 实测值

# 被检验的工程估计量
EST_ARMS = (("midpoint", dict(ctrl_on="est", estimator="midpoint")),
            ("mp3", dict(ctrl_on="est", estimator="mp3")),
            ("mp5", dict(ctrl_on="est", estimator="mp5")),
            ("volbal", dict(ctrl_on="est", estimator="volbal")))


def code(pk):
    return "A1" if pk <= BAND_MM else "A3"


def main():
    t0 = time.time()
    out = {"purpose": "E-OBS：理想全池平均能否被工程可实现的估计量替代",
           "config": {"pool_len_m": L, "ncell": NCELL, "dx_m": L / NCELL, "duration_h": 12.0,
                      "step_dq": DQ, "band_mm": BAND_MM, "kp": 1.6, "ki": 0.0016},
           "reference_arms": {}, "estimated_arms": {}, "regression": {}}

    print("--- 参照档（回归自检，应与 E-MC 逐位一致）---", flush=True)
    ref = {}
    for name, kw in REF_ARMS:
        r = V2.sve_1d_v2(base_p(), NCELL, DQ, m_slope=0.0, **kw)
        pk = r["peak_mean_mm"]
        ref[name] = {"peak_mean_mm": round(pk, 2), "final_mean_mm": round(r["final_mean_mm"], 2),
                     "judgement": code(pk)}
        exp = REF_EXPECT[name]
        ok = abs(pk - exp) <= 0.15
        out["regression"][name] = {"expect_mm": exp, "got_mm": round(pk, 2), "match": ok}
        print("[OBS-REF] %-7s peak=%8.2f mm %s  (期望 %.1f, %s)"
              % (name, pk, ref[name]["judgement"], exp, "OK" if ok else "★不一致"), flush=True)
    out["reference_arms"] = ref

    print("--- 工程可实现估计量（判定仍以真实池均峰值为准）---", flush=True)
    est = {}
    for name, kw in EST_ARMS:
        r = V2.sve_1d_v2(base_p(), NCELL, DQ, m_slope=0.0, **kw)
        pk = r["peak_mean_mm"]
        est[name] = {"peak_mean_mm": round(pk, 2), "final_mean_mm": round(r["final_mean_mm"], 2),
                     "judgement": code(pk), "est_track_err_mm": round(r["est_track_err_mm"], 2),
                     "mass_err": float("%.2e" % r["mass_balance_rel_err"])}
        print("[OBS-EST] %-9s peak=%8.2f mm %s  估计偏差<=%.1f mm  final=%7.2f"
              % (name, pk, est[name]["judgement"], r["est_track_err_mm"],
                 r["final_mean_mm"]), flush=True)
    out["estimated_arms"] = est

    ok = [k for k, v in est.items() if v["judgement"] == "A1"]
    ideal = ref["mean"]["peak_mean_mm"]
    if ok:
        verdict = ("R2 被闭合：理想全池平均峰值 %.1f mm（A1）；在 4 种工程可实现估计量中，"
                   "%s 同样判 A1（%s），最佳可实现估计偏差 %.1f mm。"
                   "⇒ 判定码翻转**不依赖理想量测**，可用%s实现。"
                   % (ideal, "、".join(ok),
                      "、".join("%.1f mm" % est[k]["peak_mean_mm"] for k in ok),
                      min(est[k]["est_track_err_mm"] for k in ok),
                      "上下游流量计（容积平衡）" if "volbal" in ok else "有限支数水位计"))
    else:
        verdict = ("R2 未能闭合：4 种工程可实现估计量**全部**判 A3（最小 %.1f mm）→ "
                   "翻转依赖理想全池平均量测，工程可落地性必须下调声明。"
                   % min(v["peak_mean_mm"] for v in est.values()))
    out["verdict"] = verdict
    out["wall_total_s"] = round(time.time() - t0, 1)

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_obs_estimator.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))
    print("回归自检:", "全部通过" if all(v["match"] for v in out["regression"].values())
          else "★有不一致，需查")
    print("A1 可实现估计量:", ok if ok else "无")
    print("verdict:", verdict)


if __name__ == "__main__":
    main()
