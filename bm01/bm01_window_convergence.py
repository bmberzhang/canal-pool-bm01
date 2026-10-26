# -*- coding: utf-8 -*-
"""
(W-win) 判定窗口收敛性：同一算例跑到 48 h，检验"窗口取多长才够"。
回应审稿意见 M12（12 h 是否截断慢发展振荡）与 M3（窗口规则的开环/闭环混淆）。

关键区分（这是修正旧表述的核心）：
  - **峰值捕获窗**：多长窗口才能把最坏 excursion 抓住？→ 用本脚本的峰值-窗口曲线回答。
  - **收敛确认窗**：开环（闸门固定）水位多久趋于终值？→ 由 τ_ol 回答（另见 bm01_tau_identify.py）。
  两者不是同一个量；旧稿用 τ_ol 去论证判定窗口，属混淆。

输出 raw/bm01_window_convergence.json
"""
import json, os, sys, time
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


V2 = load_mod("bm01_sve1d_v2.py", "bm01_sve1d_v2")

B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
L, DQ, PERIOD = 6000.0, 45.0, 300.0
NCELL = 240
DURATION_S = 172800.0   # 48 h
WINDOWS_H = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 36.0, 48.0]
BAND = 150.0


def base_p():
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": PERIOD, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": DURATION_S}


def main():
    t0 = time.time()
    out = {"case_id": "window_convergence",
           "purpose": "判定窗口收敛：C_主案例跑到48h，峰值-窗口曲线",
           "L_m": L, "dq": DQ, "ctrl_period_s": PERIOD, "dx_m": L / NCELL, "ncell": NCELL,
           "duration_h": DURATION_S / 3600.0, "band_mm": BAND, "arms": {}}

    for name, kw in [("mean", dict(ctrl_on="mean")), ("point", dict(ctrl_on="gate"))]:
        print("\n=== %s 臂（48 h）===" % name, flush=True)
        r = V2.sve_1d_v2(base_p(), NCELL, DQ, m_slope=0.0, **kw)
        s = r["series"]
        t_h = [x / 3600.0 for x in s["t"]]
        dev = [x * 1000.0 for x in s["dev_mean"]]

        curve = []
        for w in WINDOWS_H:
            vals = [(ti, abs(v)) for ti, v in zip(t_h, dev) if ti <= w]
            tp, pk = max(vals, key=lambda z: z[1])   # 注意：按 |dev| 取最大，返回 (t, |dev|)
            curve.append({"window_h": w, "peak_mm": round(pk, 2), "t_peak_h": round(tp, 2),
                          "judgement": "A1" if pk <= BAND else "A3"})
            print("  T=%5.1f h: 峰值 %8.2f mm (t=%.2f h) → %s"
                  % (w, pk, tp, curve[-1]["judgement"]), flush=True)

        full_peak = max(abs(v) for v in dev)
        # 逐小时包络（降采样，便于画"峰值随窗口收敛"图）
        env = []
        for h in range(1, int(DURATION_S / 3600.0) + 1):
            seg = [abs(v) for ti, v in zip(t_h, dev) if h - 1 < ti <= h]
            if seg:
                env.append({"h": h, "env_max_mm": round(max(seg), 2)})
        out["arms"][name] = {
            "peak_curve": curve,
            "hourly_envelope": env,
            "peak_48h_mm": round(full_peak, 2),
            "t_peak_48h_h": round(t_h[max(range(len(dev)), key=lambda i: abs(dev[i]))], 2),
            "peak_at_4h_over_48h": round(curve[1]["peak_mm"] / full_peak, 4),
            "peak_at_12h_over_48h": round(curve[5]["peak_mm"] / full_peak, 4),
            "peak_at_24h_over_48h": round(curve[7]["peak_mm"] / full_peak, 4),
            "converged_by_12h_1pct": abs(curve[5]["peak_mm"] - full_peak) / full_peak < 0.01,
        }
        print("  48 h 全局峰值 %.2f mm；12 h 窗达其 %.1f%%"
              % (full_peak, curve[5]["peak_mm"] / full_peak * 100.0), flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_window_convergence.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_s"]), flush=True)


if __name__ == "__main__":
    main()
