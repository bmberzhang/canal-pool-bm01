#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""图1 数据重跑：C_主案例（L=6000, ΔQ=45, T_ctrl=300s, dx=25m, 12h）三臂完整波形。

目的：修复图1 的图-数据-正文三方不一致。
现有图1 用 bm01_sve1d_compare.json（4h 窗口、另一算例），峰值 142.9/270.2mm 与
正文的 75.05/470.78mm（W9 C_主案例 12h）对不上，且 W9 未存波形。

本脚本用 sve_1d_v2 重跑 C_主案例 12h，输出三个臂的完整波形：
  - 0D 臂：无（0D 是另一套集总模型，此处用 mean 臂代替池均基准，另跑 0D 见注）
  - mean 臂：ctrl_on="mean"  → 池均偏差 dev_mean（判定量，峰值应≈75.05mm）
  - point 臂：ctrl_on="gate" → 池均偏差 dev_mean（判定量，峰值应≈470.78mm）
                   + 闸前水位 dev_gate（诊断量，展示"闸前单点"的振荡形态）

判定量一律用池均峰值（与 W9 一致），避免图1 再混入闸前诊断量。
产物：raw/bm01_fig1_rebuilt.json
"""
import json, os, sys, time
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# 加载 V2（bm01_sve1d_v2.py）与 V1（其依赖）
def load_mod(path, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod   # 注册进 sys.modules，dataclass 才能工作
    spec.loader.exec_module(mod)
    return mod

V2 = load_mod("bm01_sve1d_v2.py", "bm01_sve1d_v2")
PS = load_mod("bm01_pool_step.py", "bm01_pool_step")  # 0D 蓄量模型（需 numpy）

B0, KG, RATE = 50.0, 360.0, 1.0 / 3000.0
L, DQ, PERIOD = 6000.0, 45.0, 300.0
NCELL = 240  # dx = 25 m

def base_p():
    return {"pool_len_m": L, "pool_width_m": B0, "q0": 180.0, "h_target": 7.0, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": PERIOD, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}

def main():
    t0 = time.time()
    out = {"case_id": "fig1_rebuilt", "purpose":
           "图1 受控反例三臂波形（C_主案例 L=6000/ΔQ=45/T=300s/12h）",
           "L_m": L, "dq": DQ, "ctrl_period_s": PERIOD, "dx_m": 25.0, "ncell": NCELL,
           "duration_h": 12.0, "band_mm": 150.0, "arms": {}}

    # 0D 臂（集总蓄量模型，反馈=池水位，判定量=蓄量等效水位偏差）
    p0d = PS.Params(duration_s=43200, dt=5.0, t_step_s=900.0, step_dq=DQ)
    sc0d = PS.Scenario(key="main", label="0D", pools=1, step_dq=DQ,
                       feedforward=False, note="")
    r0d = PS.simulate(p0d, sc0d)
    s0d = r0d["series"]
    out["arms"]["0D"] = {
        "ctrl_on": "lumped",
        "t_h": [float(x) / 3600.0 for x in s0d["t"]],
        "dev_mean_mm": [float(x) * 1000.0 for x in s0d["h1_dev"]],
        "dev_gate_mm": None,
        "peak_mean_mm": round(float(max(abs(x) for x in s0d["h1_dev"])) * 1000.0, 2),
        "peak_gate_mm": None,
        "mass_balance_rel_err": None,
    }
    print("[fig1_rebuild] 0D 臂: peak=%.2fmm" % out["arms"]["0D"]["peak_mean_mm"], flush=True)

    for name, kw in [("mean", dict(ctrl_on="mean")), ("point", dict(ctrl_on="gate"))]:
        r = V2.sve_1d_v2(base_p(), NCELL, DQ, m_slope=0.0, **kw)
        series = r["series"]
        t_h = [x / 3600.0 for x in series["t"]]
        out["arms"][name] = {
            "ctrl_on": kw["ctrl_on"],
            "t_h": t_h,
            "dev_mean_mm": [x * 1000.0 for x in series["dev_mean"]],
            "dev_gate_mm": [x * 1000.0 for x in series["dev_gate"]],
            "peak_mean_mm": round(r["peak_mean_mm"], 2),
            "peak_gate_mm": round(r["peak_gate_mm"], 2),
            "mass_balance_rel_err": float("%.2e" % r["mass_balance_rel_err"]),
        }
        print("[fig1_rebuild] %s 臂: peak_mean=%.2fmm peak_gate=%.2fmm" %
              (name, r["peak_mean_mm"], r["peak_gate_mm"]), flush=True)

    out["wall_total_s"] = round(time.time() - t0, 1)
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_fig1_rebuilt.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("写出:", jp)
    print("0D %.2f | mean %.2f | point %.2f（正文应为 72.7 / 75.05 / 470.78）" %
          (out["arms"]["0D"]["peak_mean_mm"],
           out["arms"]["mean"]["peak_mean_mm"],
           out["arms"]["point"]["peak_mean_mm"]))

if __name__ == "__main__":
    main()
