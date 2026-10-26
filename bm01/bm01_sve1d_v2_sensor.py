# -*- coding: utf-8 -*-
"""
BM-01 一维 v2 附加诊断：**传感器通道的网格无关性**

背景（2026-10-01 发现）
----------------------
v2 此前用**单点传感器** h[-1] 送控制器。网格加密时出现**非单调、非收敛**的峰值序列：
    矩形档 120/240/480/960 格 → 142.9 / 156.6 / 158.7 / 212.2 mm
其中 960 格的 212 mm 不是物理峰：定位到它出现在 t≈216 min（物理峰在 t≈36 min），
且闸门开度末段呈**逐拍反号的增幅极限环**（0.530→0.588→0.488→0.408→0.508→… →0.808）。

机理
----
控制器**测点**和**动作点**都在最末单元。末单元蓄量 ∝ dx，网格越细蓄量越小：
    dh[-1]/dt ≈ (Q_in − Q_gate) / (B·dx·h)
⇒ 从开度 a 到测点 h[-1] 的回路增益 **∝ 1/dx**。控制周期固定 300 s，增益随 dx 升高，
采样 PI 回路越过稳定边界 → 增幅极限环 → 晚期伪峰。

修法
----
把测点改成**闸前固定物理长度内的平均水位**（`gate_avg`，默认 100 m）：
    meas = mean(h[后 ceil(100/dx) 个单元])
该平均的响应速率 (Q_in−Q_gate)/(B·L_s) **与 dx 无关** ⇒ 回路增益与网格无关。

本脚本做什么
------------
对若干断面形状，逐网格跑 `gate`（单点）与 `gate_avg`（固定长度平均）两路，输出：
  · 全程峰（peak_mean_mm）—— 是否收敛
  · 前 2 h 峰（early_peak_mm）—— 物理瞬变窗，作为交叉核对
  · 末段值（final_mean_mm）—— 看是否稳定回落
  · 闸门开度 min/max —— 判断是否出现极限环

用法
----
    python bm01_sve1d_v2_sensor.py
    python bm01_sve1d_v2_sensor.py --grids 240,480,960 --mslopes 0,1.5

输出：raw/bm01_sve1d_v2_sensor.json
证据等级：synthetic-1d-sve（诊断项，不改变 BM-01 判定码）
"""
import argparse
import importlib.util
import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


V2 = _load("V2", "bm01_sve1d_v2.py")          # V1 由 V2 内部导入


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 传感器通道网格无关性诊断")
    ap.add_argument("--case", default=".")
    ap.add_argument("--grids", default="240,480,960")
    ap.add_argument("--mslopes", default="0,1.5")
    ap.add_argument("--sensor-len", type=float, default=100.0)
    args = ap.parse_args()

    case = os.path.abspath(args.case)
    with open(os.path.join(case, "raw", "bm01_main.json"), encoding="utf-8") as fh:
        raw = json.load(fh)
    p = raw["params"]
    step_dq = float(raw["scenario"]["step_dq"])
    grids = [int(x) for x in args.grids.split(",")]
    mslopes = [float(x) for x in args.mslopes.split(",")]
    t0 = time.time()

    print("=" * 78)
    print(f"BM-01 传感器通道诊断　闸前固定长度 = {args.sensor_len:.0f} m　阶跃 +{step_dq:.0f} m³/s")
    print("=" * 78)
    print(f"{'m':>5}{'传感器':>10}{'格':>6}{'dt(s)':>8}{'全程峰':>10}{'前2h峰':>10}"
          f"{'末段':>10}{'开度min':>10}{'开度max':>10}")

    res = {}
    for m in mslopes:
        sec = V2.Section(p["pool_width_m"], m)
        s0 = V2.section_normal_slope(sec, p["h_target"], p["q0"], V2.MANNING_N)
        for chan in ("gate", "gate_avg"):
            rows = []
            for nc in grids:
                r = V2.sve_1d_v2(p, nc, step_dq, ctrl_on=chan, m_slope=m, s0=s0,
                                 sensor_len_m=args.sensor_len)
                a = r["series"]["gate"]
                dm = r["series"]["dev_mean"][-1] * 1000
                print(f"{m:>5g}{chan:>10}{nc:>6}{r['dt_s']:>8}{r['peak_mean_mm']:>10.1f}"
                      f"{r['early_peak_mm']:>10.1f}{dm:>10.1f}{min(a):>10.4f}{max(a):>10.4f}")
                rows.append({"grid": nc, "dt_s": r["dt_s"],
                             "peak_mean_mm": r["peak_mean_mm"],
                             "early_peak_mm": r["early_peak_mm"],
                             "final_mean_mm": r["final_mean_mm"],
                             "peak_gate_mm": r["peak_gate_mm"],
                             "over_band_s": r["over_band_s"],
                             "gate_min": min(a), "gate_max": max(a),
                             "mass_balance_rel_err": r["mass_balance_rel_err"]})
            res[f"m={m:g}|{chan}"] = rows
        print()

    # 收敛性判定：用「全程峰」在细网格上是否仍随网格漂移
    verdict = []
    for key, rows in res.items():
        pk = [r["peak_mean_mm"] for r in rows]
        ep = [r["early_peak_mm"] for r in rows]
        drift_pk = (max(pk) - min(pk)) / max(pk) * 100.0 if pk else 0.0
        drift_ep = (max(ep) - min(ep)) / max(ep) * 100.0 if ep else 0.0
        verdict.append({"channel": key, "peak_drift_pct": drift_pk,
                        "early_peak_drift_pct": drift_ep,
                        "converged_peak_ok": drift_pk < 5.0})
        print(f"{key:>16}：全程峰网格漂移 {drift_pk:6.2f}%　前2h峰漂移 {drift_ep:6.2f}%"
              f"　→ {'收敛 ✓' if drift_pk < 5.0 else '不收敛 ✗'}")

    payload = {
        "case_id": "BM-01",
        "purpose": "sensor_channel_grid_invariance",
        "evidence_grade": "synthetic-1d-sve",
        "note": "诊断项：只回答「传感器通道是否网格无关」，不改变 BM-01 判定码与证据等级。",
        "root_cause": "控制器测点/动作点都在最末单元，回路增益 ∝ 1/dx；"
                      "300 s 采样周期下，网格加密越过稳定边界 → 增幅极限环。",
        "fix": f"测点改为闸前固定物理长度 {args.sensor_len:.0f} m 内平均（gate_avg），"
               "响应速率分母为 B·L_s（dx 无关）⇒ 回路增益与网格无关。",
        "experiment": {"grids": grids, "mslopes": mslopes,
                       "sensor_len_m": args.sensor_len, "step_dq": step_dq},
        "results": res,
        "verdict": verdict,
        "runtime_s": round(time.time() - t0, 2),
    }
    out = os.path.join(case, "raw", "bm01_sve1d_v2_sensor.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(f"\n总用时 {time.time()-t0:.1f} s")
    print(f"输出：{out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
