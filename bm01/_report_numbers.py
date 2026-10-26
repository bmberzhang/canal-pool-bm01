# -*- coding: utf-8 -*-
"""汇报材料里要引用的全部数字，一次跑齐、自洽。

用法：PYTHONPATH=<venv>/pylibs-extra python _report_numbers.py
"""
from __future__ import annotations

import os
import sys
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bm01_pool_step import Params, build_scenarios          # noqa: E402
from bm01_single_pool_ff import id_tau, run, run_delayed, stats  # noqa: E402

base = Params()
sc = next(s for s in build_scenarios(base) if s.key == "main")
tau = id_tau(base, "dynamic")


def row(tag, r, p=base):
    s = stats(p, r["series"]["h1_dev"], r["series"]["a1"])
    print("%-34s peak=%6.1f mm  travel=%.4f  code=%s" %
          (tag, s["peak"], s["travel"], r["judgement"]["code"]))
    return s


print("== 集总 vs ID（单池，阶跃 +45 m3/s，tau=%.0f s）==" % tau)
ro = row("集总·纯反馈", run(base, sc, 0.0, 0.0))
row("集总·前馈", run(base, sc, 1.0, 0.0))
rio = row("ID·纯反馈", run_delayed(base, sc, tau, 0.0, 0.0))
row("ID·前馈 提前682s", run_delayed(base, sc, tau, 1.0, 0.0))
s382 = row("ID·前馈 提前382s", run_delayed(base, sc, tau, 1.0, 382.0))

print()
print("== 瓶颈分解（ID，基线取 ID 纯反馈 %.1f mm）==" % rio["peak"])
for tag, p, g, d in [
        ("只加快控制周期(300->5 s)", replace(base, ctrl_period_s=5.0), 0.0, 0.0),
        ("只放开闸门限速", replace(base, rate_limit_per_s=1.0), 0.0, 0.0),
        ("前馈382s（限速照旧）", base, 1.0, 382.0),
        ("前馈382s + 放开限速", replace(base, rate_limit_per_s=1.0), 1.0, 382.0)]:
    s = row(tag, run_delayed(p, sc, tau, g, d))
    print("%-34s 相对 ID 纯反馈 %+.0f%%" % ("", 100 * (1 - s["peak"] / rio["peak"]) - 100))

print()
print("== 尺度 ==")
print("自由上涨 %.1f mm/拍 = %.0f x %.0f / %.0f" %
      (sc.step_dq * base.ctrl_period_s / base.area_m2 * 1000, sc.step_dq,
       base.ctrl_period_s, base.area_m2))
print("理想前馈开度 = %.3f ; 开到位时间 = %.0f s ; tau/控制周期 = %.2f ; tau/开到位 = %.2f" %
      (sc.step_dq / base.kg, (sc.step_dq / base.kg) / base.rate_limit_per_s,
       tau / base.ctrl_period_s, tau / ((sc.step_dq / base.kg) / base.rate_limit_per_s)))

print()
print("== 真实来水（含噪）三种取法 ==")
z = np.load(os.path.join(HERE, "raw", "bm01_inflow_integral.npz"))
for k, lab in (("h_off", "纯反馈"), ("h_inst", "瞬时前馈"), ("h_win", "窗口积分前馈")):
    pk = float(np.abs(z[k]).max()) * 1000
    print("%-14s peak=%.2f mm  相对纯反馈 %+.0f%%" %
          (lab, pk, -(100 * (1 - pk / (float(np.abs(z["h_off"]).max()) * 1000)))))

print()
print("== 前馈收益 vs 传播时间（窗口/瞬时）==")
tm = np.load(os.path.join(HERE, "raw", "bm01_ff_tau_map.npz"))
for i, t in enumerate(tm["tau"]):
    print("tau=%5.0f s (%.2f 倍控制周期)  纯反馈 %.1f  瞬时 %+5.1f%%  窗口 %+5.1f%%" %
          (t, t / base.ctrl_period_s, tm["off"][i],
           -100 * (1 - tm["inst"][i] / tm["off"][i]),
           -100 * (1 - tm["win"][i] / tm["off"][i])))
print("真实来水下 纯反馈峰值区间：%.1f ~ %.1f mm" %
      (tm["off"].min(), tm["off"].max()))
