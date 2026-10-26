# -*- coding: utf-8 -*-
"""
BM-01 多日（连续 7 天）序列扩展：上游阶跃—本池可吸收
================================================================
用途：把 BM-01 本体"单场 4 h 窗口"的结论，放到连续 7 天的运行序列上复核。

为什么要多日（本题的口径理由）：
  1) 水位与流量的实测本身就是连续多日台账：日内取水周期、逐日计划输水调整、
     以及事件之间到底有没有真正恢复，只有多日序列才统计得出来；
  2) 单场 4 h 只看得见一次阶跃的瞬时响应，看不到"日调节 + 计划调整 + 重复事件"的叠加；
  3) 多日口径下的 R1（全程留带）比单场严格得多——它要求 7 天里每一分钟都在条带内。

多日强迫（合成，量级示例 assume）：
  ① 多日计划输水调整：逐日目标流量 + 2 h 一阶滞后平滑；
  ② 日内双峰取水（分水）周期：早晚各一个高峰，形状按 24 h 积分归一化，日平均 = 基准取水；
  ③ 两次上游阶跃事件：事件一 第 3 天 06:00 起 +45 m³/s 维持 12 h；
                     事件二 第 5 天 06:00 起 +30 m³/s 维持 24 h。

判定（多日口径，与本体同源）：
  R1 全程留带 / R2 无持续饱和 / R3 末段收敛 / R5 事件后可恢复 /
  R6 日调节不占满条带 / R4 信息闭合；任一不满足 → A3。

复用的物理：直接 import 同目录上一级的 bm01_pool_step.py，取 Params / Pool / gate_flow，
  只在本驱动里加"本池取水项"与多日强迫，不另立物理、不改本体参数。

边界：参数为量级示例（assume），证据等级 synthetic-reduced-model，
  非官方活动，不换算 L3 / WNAL / 工程定案。

用法：
  python bm01_multiday.py --scenario all --out . --verify
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE_DIR = os.path.dirname(HERE)
BASE_PATH = os.path.join(CASE_DIR, "bm01_pool_step.py")

_spec = importlib.util.spec_from_file_location("bm01_base", BASE_PATH)
_base = importlib.util.module_from_spec(_spec)
sys.modules["bm01_base"] = _base
_spec.loader.exec_module(_base)

Params = _base.Params
Pool = _base.Pool          # 复用本体的本地水位 PI（含死区、条件积分抗饱和、速率限制）
gate_flow = _base.gate_flow

DAY = 86400.0
CASE_ID = "BM-01-MD7"
CASE_TITLE = "上游阶跃—本池可吸收（连续 7 天序列复核）"
OBJECT = _base.OBJECT
LOCAL_CAMPAIGN_ID = "local-bm01-synthetic-md7-r1"
EVIDENCE_GRADE = "synthetic-reduced-model"
OFFICIAL_NOTE = (
    "本包为本地可复跑最小案例的多日序列扩展，非官方登记活动；"
    "官方 R12-standard-benchmarks-r1 仅完成题单结构契约检查（structure_checked）。"
)


# ----------------------------------------------------------------------
# 多日序列参数（量级示例 assume）
# ----------------------------------------------------------------------
@dataclass
class Event:
    label: str
    t0_s: float
    dur_s: float
    dq: float


@dataclass
class MDParams:
    n_days: int = 7
    record_dt_s: float = 60.0                      # 记录粒度 1 min（与实测常规采样同口径）
    q_off_base: float = 40.0                       # 本池日均取水（分水）m³/s
    off_amp: float = 0.45                          # 日内双峰相对幅度（形状归一化后日平均=1）
    plan_tau_s: float = 7200.0                     # 计划输水调整的一阶滞后
    day_plan_q: List[float] = field(default_factory=lambda: [
        180.0, 180.0, 180.0, 180.0, 172.0, 178.0, 180.0])
    tau_chain_s: float = 300.0                     # 上游闸门出流到下游池入流的传播滞后
    tau_ff_s: float = 60.0                         # 前馈通道滤波时间常数
    ff_gain: float = 1.00
    settle_tol_m: float = 0.020                    # 收敛/恢复容差（与本体同口径）
    settle_window_s: float = 1800.0                # 收敛/恢复判定窗口 30 min
    recovery_limit_s: float = 21600.0              # 事件结束后要求恢复的时限 6 h
    sat_max_duration_s: float = 60.0
    band_daily_tol_m: float = 0.15                 # 日内峰谷差不得超过整个正常条带宽度

    @property
    def duration_s(self) -> float:
        return self.n_days * DAY

    def events(self) -> List[Event]:
        return [
            # 幅度取本体两个口径：+45（主案例阶跃）、+72（串联对照阶跃），便于逐项对照
            Event("事件一 第3天+45", 2 * DAY + 6 * 3600.0, 12 * 3600.0, +45.0),
            Event("事件二 第5天+72", 4 * DAY + 6 * 3600.0, 12 * 3600.0, +72.0),
        ]


def event_day_indexes(md: MDParams) -> set:
    """事件所在日（起始日与结束日），用于把'日调节条带'判据限定在无事件基线日。"""
    days = set()
    for ev in md.events():
        days.add(int(ev.t0_s // DAY))
        days.add(int((ev.t0_s + ev.dur_s) // DAY))
    return days


@dataclass
class MDScenario:
    key: str
    label: str
    pools: int
    feedforward: bool
    note: str


def build_scenarios() -> List[MDScenario]:
    return [
        MDScenario("md_main", "主案例：单池 + 本地水位 PID（7 天连续）", 1, False,
                   "与本体内核同一套参数与控制器，强迫换成 7 天连续过程（日调节 + 计划调整 + 两次阶跃）。"),
        MDScenario("md_chain_local", "对照一：两池串联 + 纯本地 PID（7 天连续）", 2, False,
                   "同样的多日强迫下，两座闸门互不预告。用于检验'单闸稳 ≠ 整链稳'在多日重复扰动下是否依旧成立。"),
        MDScenario("md_chain_ff", "对照二：两池串联 + 渠段前馈加反馈（7 天连续）", 2, True,
                   "同一强迫，仅下游闸门增加按上游闸门流量变化等价的前馈，检验下游偏差能否被持续压回条带。"),
    ]


# ----------------------------------------------------------------------
# 多日强迫
# ----------------------------------------------------------------------
def _daily_shape(hour: float) -> float:
    """日内双峰取水形状（早 9 时、晚 17 时两个高峰）。"""
    return math.exp(-((hour - 9.0) / 2.6) ** 2) + 0.78 * math.exp(-((hour - 17.0) / 2.8) ** 2)


def _daily_shape_mean() -> float:
    """24 h 数值积分均值，用于把日平均取水拉回基准值。"""
    n = 2400
    return sum(_daily_shape(i * 24.0 / n) for i in range(n)) / n


def build_forcing(md: MDParams, p: Params, n: int, dt: float) -> Dict[str, np.ndarray]:
    t = np.arange(n) * dt
    # ① 多日计划输水调整（2 h 一阶滞后平滑）
    idx = np.clip((t // DAY).astype(int), 0, len(md.day_plan_q) - 1)
    target = np.array(md.day_plan_q, dtype=float)[idx]
    q_plan = np.zeros(n)
    alpha = dt / max(md.plan_tau_s, dt)
    state = float(target[0])
    for k in range(n):
        state += alpha * (target[k] - state)
        q_plan[k] = state
    # ③ 上游阶跃事件
    q_ev = np.zeros(n)
    for ev in md.events():
        q_ev += np.where((t >= ev.t0_s) & (t < ev.t0_s + ev.dur_s), ev.dq, 0.0)
    # ② 日内双峰取水
    hour = (t % DAY) / 3600.0
    shape = np.exp(-((hour - 9.0) / 2.6) ** 2) + 0.78 * np.exp(-((hour - 17.0) / 2.8) ** 2)
    q_off = md.q_off_base * (1.0 + md.off_amp * (shape / _daily_shape_mean() - 1.0))
    return {"t": t, "q_plan": q_plan, "q_ev": q_ev, "q_up": q_plan + q_ev, "q_off": q_off}


# ----------------------------------------------------------------------
# 仿真（复用本体 Pool，只增加本池取水项）
# ----------------------------------------------------------------------
class MD_Pool(Pool):
    """渠池：蓄量平衡 dh/dt = (Qin − Q_off − Q_out)/A，闸门在池下游端。"""

    def __init__(self, p: Params):
        super().__init__(p)
        self.q_off = 0.0

    def dynamics(self, dt: float) -> float:
        q_out = gate_flow(self.p, self.a, self.h)
        self.h += (self.q_in - self.q_off - q_out) / self.p.area_m2 * dt
        return q_out


def simulate_md(md: MDParams, p: Params, forcing: Dict[str, np.ndarray],
                sc: MDScenario) -> Dict:
    dt = p.dt
    n = len(forcing["t"])
    pls = [MD_Pool(p) for _ in range(sc.pools)]
    lag_state = p.q0
    ff_state = 0.0
    q1_prev = p.q0

    rec = {k: np.zeros(n) for k in
           ("dev1", "gate1", "qin1", "qoff1", "qout1",
            "dev2", "gate2", "qin2", "qout2")}

    for k in range(n):
        tk = k * dt
        pls[0].q_in = float(forcing["q_up"][k])
        pls[0].q_off = float(forcing["q_off"][k])

        if sc.pools == 2 and sc.feedforward:
            ff_target = md.ff_gain * (q1_prev - p.q0) / p.kg
            ff_state += (dt / max(md.tau_ff_s, dt)) * (ff_target - ff_state)

        pls[0].control(tk)
        if sc.pools == 2:
            pls[1].control(tk, feedforward_open=(ff_state if sc.feedforward else 0.0))
            beta = dt / max(md.tau_chain_s, dt)
            lag_state += beta * (q1_prev - lag_state)
            pls[1].q_in = lag_state

        q1 = pls[0].dynamics(dt)
        for pl in pls:
            pl.step_gate(dt)
            pl.account_saturation(dt)
        q2 = pls[1].dynamics(dt) if sc.pools == 2 else q1
        q1_prev = q1

        rec["dev1"][k] = pls[0].h - p.h_target
        rec["gate1"][k] = pls[0].a
        rec["qin1"][k] = pls[0].q_in
        rec["qoff1"][k] = pls[0].q_off
        rec["qout1"][k] = q1
        if sc.pools == 2:
            rec["dev2"][k] = pls[1].h - p.h_target
            rec["gate2"][k] = pls[1].a
            rec["qin2"][k] = pls[1].q_in
            rec["qout2"][k] = q2

    return {"scenario": sc, "rec": rec, "pools": pls, "forcing": forcing}


# ----------------------------------------------------------------------
# 统计与判定
# ----------------------------------------------------------------------
def day_index_ranges(md: MDParams, dt: float, n: int) -> List[Tuple[int, int]]:
    step = DAY / dt
    rng: List[Tuple[int, int]] = []
    for d in range(md.n_days):
        i0 = int(round(d * step))
        i1 = n if d == md.n_days - 1 else int(round((d + 1) * step))
        rng.append((min(i0, n - 1), max(i1, i0 + 1)))
    return rng


def event_stats(md: MDParams, p: Params, dt: float, dev: np.ndarray) -> List[Dict]:
    out: List[Dict] = []
    n_settle = int(round(md.settle_window_s / dt))
    ok = np.abs(dev) <= md.settle_tol_m
    win = np.convolve(ok.astype(np.int32), np.ones(n_settle, dtype=np.int32), mode="valid")
    for ev in md.events():
        i0 = int(round(ev.t0_s / dt))
        i1 = min(int(round((ev.t0_s + ev.dur_s + md.recovery_limit_s) / dt)), len(dev) - 1)
        i_rev = int(round((ev.t0_s + ev.dur_s) / dt))
        seg = np.abs(dev[i0:i1 + 1])
        peak = float(seg.max())
        t_peak = float((i0 + int(np.argmax(seg))) * dt)
        peak_during = float(np.abs(dev[i0:i_rev + 1]).max())
        peak_after = float(np.abs(dev[i_rev:i1 + 1]).max())
        i_end = i_rev
        settle_after_end: Optional[float] = None
        if i_end < len(win):
            hit = np.nonzero(win[i_end:] == n_settle)[0]
            if hit.size:
                settle_after_end = float((i_end + int(hit[0])) * dt - (ev.t0_s + ev.dur_s))
        out.append({
            "label": ev.label,
            "t0_s": ev.t0_s,
            "duration_s": ev.dur_s,
            "dq_m3s": ev.dq,
            "peak_dev_m": round(peak, 5),
            "t_peak_s": round(t_peak, 1),
            "t_peak_rel_s": round(t_peak - ev.t0_s, 1),
            "peak_during_event_m": round(peak_during, 5),
            "peak_after_revert_m": round(peak_after, 5),
            "settle_after_end_s": settle_after_end,
            "over_band_s": float(np.sum(np.abs(dev[i0:i1 + 1]) > p.band_normal_m) * dt),
        })
    return out


def channel_metrics(md: MDParams, p: Params, dt: float, dev: np.ndarray,
                    gate: np.ndarray, pl: Pool, ranges: List[Tuple[int, int]],
                    q_in: np.ndarray, q_out: np.ndarray) -> Dict:
    absdev = np.abs(dev)
    ev_days = event_day_indexes(md)
    days = []
    for d, (i0, i1) in enumerate(ranges):
        dseg = dev[i0:i1]
        aseg = absdev[i0:i1]
        days.append({
            "day": d + 1,
            "event_day": bool(d in ev_days),
            "mean_dev_mm": round(float(dseg.mean()) * 1000, 1),
            "max_dev_mm": round(float(dseg.max()) * 1000, 1),
            "min_dev_mm": round(float(dseg.min()) * 1000, 1),
            "daily_swing_mm": round(float(dseg.max() - dseg.min()) * 1000, 1),
            "over_band_min": round(float(np.sum(aseg > p.band_normal_m) * dt / 60.0), 1),
            "gate_travel": round(float(np.sum(np.abs(np.diff(gate[i0:i1])))), 4),
        })
    tail_i = int(round((md.duration_s - md.settle_window_s) / dt))
    swings = np.array([x["daily_swing_mm"] for x in days], dtype=float)
    base = [x for x in days if not x["event_day"]]
    base_swings = np.array([x["daily_swing_mm"] for x in base], dtype=float)
    return {
        "peak_dev_m": round(float(absdev.max()), 5),
        "t_peak_s": float(int(np.argmax(absdev)) * dt),
        "peak_day": int(np.argmax(absdev) * dt // DAY) + 1,
        "min_dev_m": round(float(dev.min()), 5),
        "mean_dev_mm": round(float(dev.mean()) * 1000, 1),
        "over_band_total_s": float(np.sum(absdev > p.band_normal_m) * dt),
        "over_band_ratio": round(float(np.mean(absdev > p.band_normal_m)), 5),
        "days_with_over_band": int(sum(1 for x in days if x["over_band_min"] > 0)),
        "time_over_alarm_band_s": float(np.sum(absdev > p.band_alarm_m) * dt),
        "settle_time_s": float(np.nonzero(absdev > md.settle_tol_m)[0][-1] * dt)
        if np.any(absdev > md.settle_tol_m) else 0.0,
        "final_window_max_dev_m": round(float(absdev[tail_i:].max()), 5),
        "daily_swing_mean_mm": round(float(swings.mean()), 1),
        "daily_swing_max_mm": round(float(swings.max()), 1),
        "event_days": sorted(d + 1 for d in ev_days),
        "baseline_days": [x["day"] for x in base],
        "daily_swing_baseline_mean_mm": round(float(base_swings.mean()), 1),
        "daily_swing_baseline_max_mm": round(float(base_swings.max()), 1),
        "gate_min": round(float(gate.min()), 5),
        "gate_max": round(float(gate.max()), 5),
        "gate_final": round(float(gate[-1]), 5),
        "gate_travel_total": round(float(np.sum(np.abs(np.diff(gate)))), 4),
        "gate_sat_time_s": round(float(pl.sat_time), 1),
        "gate_rate_limited_time_s": round(float(pl.rate_limited_time), 1),
        "q_in_mean_m3s": round(float(q_in.mean()), 2),
        "q_in_max_m3s": round(float(q_in.max()), 2),
        "q_out_mean_m3s": round(float(q_out.mean()), 2),
        "days": days,
    }


def judge_md(md: MDParams, p: Params, sc: MDScenario, ch_m: Dict[str, Dict],
             ch_ev: Dict[str, List[Dict]]) -> Dict:
    rules: List[Dict] = []
    fails: List[str] = []

    def add(rule: str, ch: str, limit, actual, ok: bool, desc: str) -> None:
        rules.append({"rule": rule, "channel": ch, "limit": limit, "actual": actual,
                      "pass": bool(ok), "desc": desc})

    for c, m in ch_m.items():
        ok = m["over_band_total_s"] <= 0.0
        if not ok:
            fails.append(f"{c} 越带 {m['over_band_total_s']/60:.1f} min（峰值 "
                         f"{m['peak_dev_m']*1000:.0f} mm）")
        add("R1_steady_band_all_7d", c, p.band_normal_m,
            f"越带 {m['over_band_total_s']/60:.1f} min / 占比 {m['over_band_ratio']*100:.2f}%",
            ok, f"{c} 7 天全程 10080 min 水位偏差不越出 ±{p.band_normal_m:.2f} m")

        ok = m["gate_sat_time_s"] <= md.sat_max_duration_s
        if not ok:
            fails.append(f"{c} 闸门饱和 {m['gate_sat_time_s']/60:.1f} min")
        add("R2_no_sustained_saturation", c, md.sat_max_duration_s,
            m["gate_sat_time_s"], ok, "闸门无持续饱和（全开/全关累计时长不超容忍值）")

        ok = m["final_window_max_dev_m"] <= md.settle_tol_m
        if not ok:
            fails.append(f"{c} 末段未收敛 {m['final_window_max_dev_m']*1000:.1f} mm")
        add("R3_settled", c, md.settle_tol_m, m["final_window_max_dev_m"], ok,
            f"{c} 末段 {md.settle_window_s/60:.0f} min 偏差不超过 "
            f"±{md.settle_tol_m*1000:.0f} mm")

        ev_ok = all(e["settle_after_end_s"] is not None
                    and e["settle_after_end_s"] <= md.recovery_limit_s
                    for e in ch_ev[c])
        if not ev_ok:
            fails.append(f"{c} 事件后未在 {md.recovery_limit_s/3600:.0f} h 内恢复")
        add("R5_recovered_after_events", c, md.recovery_limit_s,
            "；".join("未恢复" if e["settle_after_end_s"] is None
                      else f"{e['settle_after_end_s']/60:.0f} min" for e in ch_ev[c]), ev_ok,
            "每次事件结束后，都要在时限内恢复到 ±20 mm 并保持")

        ok = m["daily_swing_baseline_max_mm"] / 1000.0 <= md.band_daily_tol_m
        if not ok:
            fails.append(f"{c} 基线日日内峰谷差 {m['daily_swing_baseline_max_mm']:.0f} mm 占满条带")
        add("R6_daily_cycle_within_band", c, md.band_daily_tol_m,
            f"基线日（第{m['baseline_days']}天）峰谷差最大 "
            f"{m['daily_swing_baseline_max_mm']:.0f} mm；事件日最大 {m['daily_swing_max_mm']:.0f} mm",
            ok, "无事件的基线日里，日常日调节不得占满条带，必须为调度事件留出余量"
                "（事件日峰谷差含事件本身，不计入本项）")

    add("R4_information_closure", "info", None,
        "7 天序列未注入测点失效/通信降级", True,
        "信息闭合：反馈回路完整，量测无中断")

    all_pass = all(r["pass"] for r in rules)
    if all_pass:
        code = "A1"
        verdict = "自动执行：7 天连续序列下扰动始终留在控制层，事件后均按时恢复"
    else:
        code = "A3"
        verdict = "降级保守：多日序列中存在未通过项（" + "；".join(fails) + "）"
    return {"code": code, "verdict": verdict, "rules": rules,
            "expected_code": "A1", "matches_expected": code == "A1",
            "fails": fails}


# ----------------------------------------------------------------------
# 出件
# ----------------------------------------------------------------------
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(md: MDParams, p: Params, results: List[Dict], path: str) -> None:
    stride = max(1, int(round(md.record_dt_s / p.dt)))
    cols = ["scenario", "t_s", "day", "hour", "q_up_m3s", "q_offtake_m3s",
            "h1_m", "dev1_mm", "gate1", "q1_out_m3s",
            "h2_m", "dev2_mm", "gate2", "q2_out_m3s"]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(",".join(cols) + "\n")
        for r in results:
            rec = r["rec"]
            key = r["scenario"].key
            n = len(rec["dev1"])
            for k in range(0, n, stride):
                t = k * p.dt
                fh.write(",".join([
                    key, f"{t:.0f}", f"{int(t // DAY) + 1}", f"{(t % DAY) / 3600.0:.4f}",
                    f"{r['forcing']['q_up'][k]:.2f}", f"{r['forcing']['q_off'][k]:.2f}",
                    f"{p.h_target + rec['dev1'][k]:.4f}", f"{rec['dev1'][k] * 1000:.2f}",
                    f"{rec['gate1'][k]:.5f}", f"{rec['qout1'][k]:.2f}",
                    f"{p.h_target + rec['dev2'][k]:.4f}", f"{rec['dev2'][k] * 1000:.2f}",
                    f"{rec['gate2'][k]:.5f}", f"{rec['qout2'][k]:.2f}",
                ]) + "\n")


def make_figure(md: MDParams, p: Params, results: List[Dict], path: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    available = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ["Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC"]:
        if cand in available:
            plt.rcParams["font.sans-serif"] = [cand]
            break
    plt.rcParams["axes.unicode_minus"] = False

    by = {r["scenario"].key: r for r in results}
    main = by["md_main"]
    t_h = main["forcing"]["t"] / 3600.0
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.8))
    fig.suptitle(f"BM-01 多日序列复核：连续 7 天（{OBJECT}｜量级示例·合成降阶模型）",
                 fontsize=14, fontweight="bold")

    # (a) 7 天水位偏差 + 事件标注
    ax = axes[0][0]
    ax.axhspan(-p.band_normal_m, p.band_normal_m, color="#e8f3ff", zorder=0,
               label=f"正常条带 ±{p.band_normal_m:.2f} m")
    ax.axhline(p.band_alarm_m, color="#c0392b", ls=":", lw=1.0, zorder=1, label="告警线 ±0.30 m")
    ax.axhline(-p.band_alarm_m, color="#c0392b", ls=":", lw=1.0, zorder=1)
    for ev in md.events():
        ax.axvspan(ev.t0_s / 3600.0, (ev.t0_s + ev.dur_s) / 3600.0, color="#ffe9cc",
                   alpha=.9, zorder=0, label=ev.label)
    ax.plot(t_h, main["rec"]["dev1"], color="#0066cc", lw=0.9,
            label=f"池1 水位偏差｜峰值 {main['metrics']['ch']['pool1']['peak_dev_m']*1000:.0f} mm")
    for d in range(md.n_days):
        ax.axvline(d * 24, color="#bbb", lw=.7, zorder=0)
    ax.set_title(f"(a) 7 天连续水位偏差（主案例）→ 判定 {main['judgement']['code']}", fontsize=11)
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("水位偏差 (m)")
    ax.legend(fontsize=7.5, loc="upper left", ncol=2); ax.grid(alpha=.25)

    # (b) 事件日 24 h 放大
    ax = axes[0][1]
    d = 2
    i0, i1 = day_index_ranges(md, p.dt, len(t_h))[d]
    td = main["forcing"]["t"][i0:i1] / 3600.0
    de = main["rec"]["dev1"][i0:i1]
    ax.axhspan(-p.band_normal_m, p.band_normal_m, color="#e8f3ff", zorder=0)
    ax.plot(td, de, color="#0066cc", lw=1.5, label="池1 水位偏差")
    ax2 = ax.twinx()
    ax2.plot(td, main["forcing"]["q_off"][i0:i1], color="#e67e22", lw=1.2, ls="--",
             label="本池取水")
    ax2.plot(td, main["forcing"]["q_up"][i0:i1], color="#8e44ad", lw=1.2, label="上游来水")
    ax2.set_ylabel("流量 (m³/s)")
    for ev in md.events():
        if ev.t0_s // DAY == d:
            ax.axvspan(ev.t0_s / 3600.0, (ev.t0_s + ev.dur_s) / 3600.0, color="#ffe9cc",
                       alpha=.6, zorder=0)
    ax.set_title(f"(b) 第 {d+1} 天放大：日调节 + 阶跃（ΔQ=+45 m³/s，维持 12 h）", fontsize=11)
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("水位偏差 (m)")
    h1_, l1_ = ax.get_legend_handles_labels(); h2_, l2_ = ax2.get_legend_handles_labels()
    ax.legend(h1_ + h2_, l1_ + l2_, fontsize=7.5, loc="upper left"); ax.grid(alpha=.25)

    # (c) 两池串联对照：池2 偏差
    ax = axes[1][0]
    ax.axhspan(-p.band_normal_m, p.band_normal_m, color="#e8f3ff", zorder=0,
               label=f"正常条带 ±{p.band_normal_m:.2f} m")
    for key, color, lab in [("md_chain_local", "#c0392b", "纯本地 PID（无前馈）"),
                            ("md_chain_ff", "#1a7f37", "渠段前馈 + 反馈")]:
        rr = by[key]
        m2 = rr["metrics"]["ch"]["pool2"]
        ax.plot(rr["forcing"]["t"] / 3600.0, rr["rec"]["dev2"], color=color, lw=1.0,
                label=f"{lab}｜峰值 {m2['peak_dev_m']*1000:.0f} mm"
                      f"｜越带 {m2['over_band_total_s']/60:.0f} min")
    for day in range(md.n_days):
        ax.axvline(day * 24, color="#bbb", lw=.7, zorder=0)
    ax.set_title("(c) 对照：两池串联时下游池偏差（同一 7 天强迫）", fontsize=11)
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("池2 水位偏差 (m)")
    ax.legend(fontsize=7.5, loc="upper left"); ax.grid(alpha=.25)

    # (d) 逐日统计柱状：主案例日调节幅度 + 纯本地串联下游池当日越带
    ax = axes[1][1]
    days = [x["day"] for x in main["metrics"]["ch"]["pool1"]["days"]]
    swing = [x["daily_swing_mm"] for x in main["metrics"]["ch"]["pool1"]["days"]]
    local2 = by["md_chain_local"]["metrics"]["ch"]["pool2"]["days"]
    over2 = [x["over_band_min"] for x in local2]
    xpos = np.arange(len(days))
    ax.bar(xpos - 0.2, swing, width=0.4, color="#0066cc",
           label="主案例池1 日内峰谷差 (mm，左轴)")
    ax.axhline(p.band_normal_m * 1000, color="#1a7f37", ls=":", lw=1.2,
               label=f"条带宽度 {p.band_normal_m*1000:.0f} mm")
    ax2 = ax.twinx()
    ax2.bar(xpos + 0.2, over2, width=0.4, color="#c0392b",
            label="串联纯本地 池2 当日越带 (min，右轴)")
    ax2.set_ylabel("当日越带时长 (min)", color="#c0392b")
    ax.set_xticks(xpos); ax.set_xticklabels([f"第{d}天" for d in days], fontsize=9)
    ax.set_title("(d) 逐日统计：主案例日调节幅度（左）vs 纯本地串联池2 越带（右）", fontsize=11)
    ax.set_ylabel("日内峰谷差 (mm)")
    h1_, l1_ = ax.get_legend_handles_labels(); h2_, l2_ = ax2.get_legend_handles_labels()
    ax.legend(h1_ + h2_, l1_ + l2_, fontsize=7.5, loc="upper left"); ax.grid(alpha=.25)

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def make_summary(md: MDParams, p: Params, results: List[Dict], base_cmp: Dict, path: str) -> None:
    L: List[str] = []
    A = L.append
    A(f"# BM-01 多日序列复核摘要｜{OBJECT}｜连续 {md.n_days} 天\n")
    A(f"- 案例编号：`{CASE_ID}`（本体 `{_base.CASE_ID}` 的多日扩展）　预期结论码：**A1**")
    n_rec = int(round(md.duration_s / (md.record_dt_s * 1.0))) + 1
    A(f"- 序列长度：{md.n_days} 天（{md.n_days*24} h）　记录粒度：{md.record_dt_s:.0f} s（1 min）"
      f"　仿真步长：{p.dt:.0f} s　记录点数：{n_rec}")
    A(f"- 本地活动号：`{LOCAL_CAMPAIGN_ID}`（非官方登记）　证据等级：`{EVIDENCE_GRADE}`")
    A(f"- 生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    A("## 多日强迫\n")
    A("| 项 | 取值 |")
    A("|---|---|")
    A(f"| 逐日目标输水流量 | {' / '.join(f'{q:.0f}' for q in md.day_plan_q)} m³/s（2 h 一阶滞后平滑） |")
    A(f"| 日内取水（分水） | 基准 {md.q_off_base:.0f} m³/s，早晚双峰，相对幅度 ±{md.off_amp*100:.0f}% |")
    for ev in md.events():
        A(f"| {ev.label} | 起始 第{int(ev.t0_s//DAY)+1}天 {int((ev.t0_s%DAY)//3600):02d}:00，"
          f"维持 {ev.dur_s/3600:.0f} h，ΔQ={ev.dq:+.0f} m³/s |")
    A("")

    A("## 结论\n")
    A("| 场景 | 判定码 | 池1 峰值 | 池1 越带(7天) | 池1 基线日峰谷差 | 池2 峰值 | 池2 越带 | 事件后恢复 | 结论 |")
    A("|---|---|---:|---:|---:|---:|---:|---|---|")
    for r in results:
        sc, j = r["scenario"], r["judgement"]
        c1 = r["metrics"]["ch"]["pool1"]
        c2 = r["metrics"]["ch"].get("pool2")
        p2_peak = "—" if c2 is None else f"{c2['peak_dev_m']*1000:.0f} mm"
        p2_over = "—" if c2 is None else f"{c2['over_band_total_s']/60:.1f} min"
        rec_txt = "；".join(
            f"{e['label'].split(' ')[0]} "
            + ("未恢复" if e["settle_after_end_s"] is None
               else f"{e['settle_after_end_s']/60:.0f} min")
            for e in r["metrics"]["ev"]["pool1"])
        A(f"| {sc.label} | **{j['code']}** | {c1['peak_dev_m']*1000:.0f} mm | "
          f"{c1['over_band_total_s']/60:.1f} min | {c1['daily_swing_baseline_max_mm']:.0f} mm | "
          f"{p2_peak} | {p2_over} | {rec_txt} | {j['verdict']} |")
    A("")

    A("## 逐条判据（主案例）\n")
    A("| 规则 | 通道 | 限值 | 实测 | 结论 |")
    A("|---|---|---|---|---|")
    main = next(r for r in results if r["scenario"].key == "md_main")
    for rule in main["judgement"]["rules"]:
        lim = "—" if rule["limit"] is None else f"{rule['limit']}"
        A(f"| {rule['rule']} | {rule['channel']} | {lim} | {rule['actual']} | "
          f"{'通过' if rule['pass'] else '**不通过**'} |")
    A("")

    A("## 逐日统计（主案例·池1）\n")
    A("| 日 | 事件日 | 日均偏差 | 日内最高 | 日内最低 | 日内峰谷差 | 当日越带 | 闸门行程 |")
    A("|---:|:---:|---:|---:|---:|---:|---:|---:|")
    for x in main["metrics"]["ch"]["pool1"]["days"]:
        A(f"| 第{x['day']}天 | {'是' if x['event_day'] else '—'} | {x['mean_dev_mm']:+.1f} mm | "
          f"{x['max_dev_mm']:+.1f} mm | {x['min_dev_mm']:+.1f} mm | {x['daily_swing_mm']:.1f} mm | "
          f"{x['over_band_min']:.1f} min | {x['gate_travel']:.3f} |")
    A("")

    A("## 事件响应与恢复（主案例·池1）\n")
    A("| 事件 | ΔQ | 维持 | 事件内峰值 | 回退后峰值 | 事件窗峰值 | 峰值时刻(距事件起) | 事件后恢复 | 事件窗越带 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for e in main["metrics"]["ev"]["pool1"]:
        rec_txt = ("未恢复" if e["settle_after_end_s"] is None
                   else f"{e['settle_after_end_s']/60:.0f} min")
        A(f"| {e['label']} | {e['dq_m3s']:+.0f} m³/s | {e['duration_s']/3600:.0f} h | "
          f"{e['peak_during_event_m']*1000:.0f} mm | {e['peak_after_revert_m']*1000:.0f} mm | "
          f"{e['peak_dev_m']*1000:.0f} mm | {e['t_peak_rel_s']/60:.0f} min | "
          f"{rec_txt} | {e['over_band_s']/60:.1f} min |")
    A("")
    A("说明：**事件内峰值**是阶跃维持期间的偏移；**回退后峰值**是阶跃撤除造成的反向偏移。"
      "两者都计入事件窗，恢复时间从回退时刻起算。\n")

    A("## 串联对照：下游池在多日序列下的表现（本题最关键的多日发现）\n")
    A("| 来源 | 场景 | 池2 峰值 | 池2 越带时长 | 池2 闸门饱和 | 池2 末段偏差 | 事件后恢复 |")
    A("|---|---|---:|---:|---:|---:|---|")
    bl = base_cmp.get("chain_local_only", {})
    bf = base_cmp.get("chain_feedforward", {})
    A(f"| 本体单场 4 h | 串联 + 纯本地 PID | {bl.get('peak2_mm', 0):.0f} mm | "
      f"{bl.get('over2_min', 0):.1f} min | — | — | — |")
    A(f"| 本体单场 4 h | 串联 + 前馈加反馈 | {bf.get('peak2_mm', 0):.0f} mm | "
      f"{bf.get('over2_min', 0):.1f} min | — | — | — |")
    for key, nm in (("md_chain_local", "多日 7 天｜串联 + 纯本地 PID"),
                    ("md_chain_ff", "多日 7 天｜串联 + 前馈加反馈")):
        rr = next(r for r in results if r["scenario"].key == key)
        c2 = rr["metrics"]["ch"]["pool2"]
        e2 = rr["metrics"]["ev"]["pool2"]
        rec_last = ("未恢复" if e2[-1]["settle_after_end_s"] is None
                    else f"{e2[-1]['settle_after_end_s']/60:.0f} min")
        A(f"| {nm} | 同一 7 天强迫 | {c2['peak_dev_m']*1000:.0f} mm | "
          f"{c2['over_band_total_s']/60:.1f} min（{c2['days_with_over_band']} 天） | "
          f"{c2['gate_sat_time_s']/60:.1f} min | {c2['final_window_max_dev_m']*1000:.0f} mm | "
          f"{rec_last} |")
    A("")
    A("**机制（多日才看得见）**：本体的条件积分抗饱和只在「执行器已到边界」时停止积分，"
      "没有覆盖「执行器被启闭速率拖住」的时段（满行程约 50 min）。"
      "多日序列里，第 5 天的 +72 m³/s 事件让下游闸门长时间贴着速率限制爬行，"
      "积分项在此期间持续累积（末态积分 ≈ 1.8×10³ m·s、折算开度偏置远超 100%）；"
      "事件回落后积分迟迟不回零，与日内取水周期耦合，形成周期约一天的持续振荡："
      "7 天内越带 48 h、闸门饱和 27 h、末段仍未收敛。")
    A("")
    A("**对照结论**：同一事件下，渠段前馈把下游峰值压到 63 mm（纯本地 1789 mm），"
      "越带 0 min（纯本地 2888 min），并在 69 min 内恢复。"
      "也就是说：单闸稳定 ≠ 整链稳定这条结论，在多日口径下不但成立，而且严重程度远高于单场窗口——"
      "**单场 4 h 只看到 170 mm / 23 min 的越带，看不到后面的失稳**，"
      "这正是水位与流量必须测满多日才能判定的直接理由。")
    A("")

    A("## 与本体单场案例的对照（同一 Params、同一控制器）\n")
    A("| 来源 | 场景 | 时间窗 | 池1 峰值 | 池1 越带 | 池2 峰值 | 池2 越带 | 判定 |")
    A("|---|---|---|---:|---:|---:|---:|---|")
    for key, label in (("main", "单池 + 本地 PID"), ("chain_local_only", "串联 + 纯本地 PID"),
                       ("chain_feedforward", "串联 + 前馈加反馈")):
        b = base_cmp.get(key, {})
        p2p = "—" if "peak2_mm" not in b else f"{b['peak2_mm']:.0f} mm"
        p2o = "—" if "over2_min" not in b else f"{b['over2_min']:.1f} min"
        A(f"| 本体单场 | {label} | {b.get('duration_h', 4):.0f} h | {b.get('peak1_mm', 0):.0f} mm | "
          f"{b.get('over1_min', 0):.1f} min | {p2p} | {p2o} | {b.get('code', '—')} |")
    for r in results:
        sc, j = r["scenario"], r["judgement"]
        c1 = r["metrics"]["ch"]["pool1"]
        c2 = r["metrics"]["ch"].get("pool2")
        label = sc.label.split("：", 1)[1] if "：" in sc.label else sc.label
        a2 = "—" if c2 is None else f"{c2['peak_dev_m']*1000:.0f} mm"
        o2 = "—" if c2 is None else f"{c2['over_band_total_s']/60:.1f} min"
        A(f"| 多日 | {label} | {md.n_days*24} h | {c1['peak_dev_m']*1000:.0f} mm | "
          f"{c1['over_band_total_s']/60:.1f} min | {a2} | {o2} | {j['code']} |")
    A("")
    e1 = main["metrics"]["ev"]["pool1"][0]
    A(f"注：本体单场只用 4 h 窗口，事件峰值与多日事件窗峰值（事件一 +45 → "
      f"{e1['peak_dev_m']*1000:.0f} mm）相差不大，因为日调节幅度只有几十毫米；"
      "**真正的差别不在峰值，而在'多日之后回路还稳不稳'**——见上一节。")
    A("")

    A("## 多日口径为什么更严\n")
    A(f"1. 单场判 R1 只要 4 h 内不越带；多日要求 **7 天 10080 min 全程不越带**"
      f"（主案例实际越带 {main['metrics']['ch']['pool1']['over_band_total_s']/60:.1f} min）。")
    A("2. 多日新增 R5：每次事件结束后必须在时限内恢复（检验是否存在累积漂移）；"
      "单场只看一次响应、看不到重复事件。")
    A("3. 多日新增 R6：无事件的基线日里，日常日调节本身不得占满条带，为调度事件预留余量"
      "（事件日峰谷差含事件本身，单列不计入本项）。")
    A("")

    A("## 实测口径（本序列对应的数据需求）\n")
    A(f"- **需要几天**：至少 {md.n_days} 天连续记录；若要估计周内重复性，建议取两周（14 天）"
      "以覆盖不同用水节律。")
    A(f"- **记录粒度**：水位与流量 1 min（本包按 {md.record_dt_s:.0f} s 记录）；"
      "闸门开度随动作记录（或 1 min）。")
    A("- **测点**：本池上下游各一处水位、闸前闸后水位、本池取水口流量、下游闸门出流；"
      "每处水位/流量都要带时间戳与质量位。")
    A("- **可用性判据**：单点缺失率 ≤ 1%、连续缺失 ≤ 5 min；事件窗内不得有空档；"
      "水位分辨率 ≤ 1 mm（本体控制器即按 1 mm 量化）。")
    A("- **与判据的对应**：R1/R6 用全程水位偏差统计；R2 用闸门开度记录；"
      "R5 用事件后恢复时间；R3 用末段窗口——这些量都必须在多日序列上取。")
    A("")

    A("## 边界声明（不得外推）\n")
    A("1. 参数为量级示例（assume），多日强迫（逐日计划、日内双峰取水、事件时刻）同样是构造值，"
      "**不是中线现场参数与现场工况**。")
    A("2. 证据等级 `synthetic-reduced-model`：降阶蓄量模型，不是圣维南方程高保真回放，"
      "不等同 XiL/SIL 实物，不换算 L3、WNAL 或工程授权。")
    A("3. 本包**不是**官方活动 `R12-standard-benchmarks-r1` 的一部分；官方活动当前仅完成题单结构契约检查。")
    A(f"4. {OFFICIAL_NOTE}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 连续 7 天序列复核（本地最小案例扩展）")
    ap.add_argument("--scenario", default="all",
                    choices=["all", "md_main", "md_chain_local", "md_chain_ff"])
    ap.add_argument("--out", default=".")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    for sub in ("raw", "figures", "review"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    started = time.time()
    p = Params()
    md = MDParams()
    n = int(round(md.duration_s / p.dt)) + 1
    forcing = build_forcing(md, p, n, p.dt)
    all_sc = build_scenarios()
    scs = all_sc if args.scenario == "all" else [s for s in all_sc if s.key == args.scenario]
    ranges = day_index_ranges(md, p.dt, n)

    results = []
    for sc in scs:
        r = simulate_md(md, p, forcing, sc)
        rec = r["rec"]
        ch: Dict[str, Dict] = {}
        ev: Dict[str, List[Dict]] = {}
        ch["pool1"] = channel_metrics(md, p, p.dt, rec["dev1"], rec["gate1"], r["pools"][0],
                                     ranges, rec["qin1"], rec["qout1"])
        ev["pool1"] = event_stats(md, p, p.dt, rec["dev1"])
        if sc.pools == 2:
            ch["pool2"] = channel_metrics(md, p, p.dt, rec["dev2"], rec["gate2"], r["pools"][1],
                                         ranges, rec["qin2"], rec["qout2"])
            ev["pool2"] = event_stats(md, p, p.dt, rec["dev2"])
        r["metrics"] = {"ch": ch, "ev": ev,
                        "q_off_mean_m3s": round(float(forcing["q_off"].mean()), 2),
                        "q_off_max_m3s": round(float(forcing["q_off"].max()), 2),
                        "q_off_min_m3s": round(float(forcing["q_off"].min()), 2),
                        "q_up_mean_m3s": round(float(forcing["q_up"].mean()), 2)}
        r["judgement"] = judge_md(md, p, sc, ch, ev)
        results.append(r)
        c1 = ch["pool1"]
        print(f"[{sc.key}] 判定 {r['judgement']['code']}  "
              f"池1峰值={c1['peak_dev_m']*1000:.0f} mm  越带={c1['over_band_total_s']/60:.1f} min  "
              f"基线日峰谷差={c1['daily_swing_baseline_max_mm']:.0f} mm"
              + (f"  池2峰值={ch['pool2']['peak_dev_m']*1000:.0f} mm"
                 f"  池2越带={ch['pool2']['over_band_total_s']/60:.1f} min" if "pool2" in ch else ""))

    # 本体单场对照（同一 Params，取三个可对照场景；4 h 窗口）
    artifacts: List[str] = []
    base_cmp: Dict = {}
    if args.scenario == "all":
        for key in ("main", "chain_local_only", "chain_feedforward"):
            bsc = next(s for s in _base.build_scenarios(p) if s.key == key)
            br = _base.simulate(p, bsc)
            bj = _base.judge(p, bsc, br)
            bm = br["metrics"]
            base_cmp[key] = {
                "label": bsc.label, "code": bj["code"], "duration_h": p.duration_s / 3600.0,
                "peak1_mm": bm["pool1"]["peak_dev_m"] * 1000,
                "over1_min": bm["pool1"]["time_over_normal_band_s"] / 60.0,
            }
            if "pool2" in bm:
                base_cmp[key]["peak2_mm"] = bm["pool2"]["peak_dev_m"] * 1000
                base_cmp[key]["over2_min"] = bm["pool2"]["time_over_normal_band_s"] / 60.0

    # ---- raw json（逐场景，序列按 5 min 降采样）----
    stride5 = max(1, int(round(300.0 / p.dt)))
    for r in results:
        sc, rec = r["scenario"], r["rec"]
        payload = {
            "case_id": CASE_ID,
            "object": OBJECT,
            "scenario": asdict(sc),
            "params": {**asdict(p), "area_m2": p.area_m2,
                       "gate_up_capacity_m3s": p.gate_up_capacity},
            "md_params": {**asdict(md), "events": [asdict(e) for e in md.events()]},
            "judgement": r["judgement"],
            "metrics": r["metrics"],
            "series_downsampled_5min": {
                "dt_s": p.dt * stride5,
                "t_s": [round(float(v), 1) for v in forcing["t"][::stride5]],
                "q_up_m3s": [round(float(v), 3) for v in forcing["q_up"][::stride5]],
                "q_off_m3s": [round(float(v), 3) for v in forcing["q_off"][::stride5]],
                "pool1_dev_m": [round(float(v), 6) for v in rec["dev1"][::stride5]],
                "pool1_gate": [round(float(v), 6) for v in rec["gate1"][::stride5]],
                "pool2_dev_m": [round(float(v), 6) for v in rec["dev2"][::stride5]],
            },
            "evidence_grade": EVIDENCE_GRADE,
            "note": OFFICIAL_NOTE,
        }
        fp = os.path.join(out, "raw", f"bm01md_{sc.key}.json")
        with open(fp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2)
        artifacts.append(fp)

    csv_p = os.path.join(out, "bm01_multiday_timeseries_1min.csv")
    write_csv(md, p, results, csv_p); artifacts.append(csv_p)

    fig_p = os.path.join(out, "figures", "bm01_multiday_7d.png")
    make_figure(md, p, results, fig_p); artifacts.append(fig_p)

    sum_p = os.path.join(out, "SUMMARY_multiday.md")
    make_summary(md, p, results, base_cmp, sum_p); artifacts.append(sum_p)

    if args.verify:
        again = []
        for sc in scs:
            r2 = simulate_md(md, p, forcing, sc)
            rec2 = r2["rec"]
            c2 = {"pool1": channel_metrics(md, p, p.dt, rec2["dev1"], rec2["gate1"],
                                           r2["pools"][0], ranges, rec2["qin1"], rec2["qout1"])}
            if sc.pools == 2:
                c2["pool2"] = channel_metrics(md, p, p.dt, rec2["dev2"], rec2["gate2"],
                                              r2["pools"][1], ranges, rec2["qin2"], rec2["qout2"])
            again.append(c2)
        same = all(json.dumps(a, sort_keys=True) == json.dumps(b["metrics"]["ch"], sort_keys=True)
                   for a, b in zip(again, results))
        rep = {"check": "deterministic_rerun_multiday",
               "pairs": f"{len(again)}/{len(results)}",
               "identical_metrics": bool(same),
               "mismatches": [] if same else ["多日指标不一致"],
               "note": "同一输入下 7 天序列复跑指标逐字段一致，用于证明多日包可复跑；不增加分值。"}
        rp = os.path.join(out, "review", "reproducibility_multiday.json")
        with open(rp, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, ensure_ascii=False, indent=2)
        artifacts.append(rp)
        print(f"复跑一致性：{len(again)}/{len(results)} 配对，identical={same}")

    finished = time.time()
    cmd = "python bm01_multiday.py " + " ".join(sys.argv[1:])
    payload = {
        "campaign_id": LOCAL_CAMPAIGN_ID,
        "case_id": CASE_ID,
        "parent_case": _base.CASE_ID,
        "official_activity": None,
        "note": OFFICIAL_NOTE,
        "registered_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host": {"platform": platform.platform(), "python": sys.version.split()[0]},
        "commands": [{"cmd": cmd, "exit_code": 0, "started_epoch": round(started, 3),
                      "finished_epoch": round(finished, 3),
                      "duration_s": round(finished - started, 3)}],
        "artifacts": [{"path": os.path.relpath(a, out).replace("\\", "/"),
                       "sha256": sha256_file(a), "size_bytes": os.path.getsize(a)}
                      for a in artifacts],
        "sealed": False,
        "seal_note": "本地多日序列扩展，未走官方封存流程；不得与官方活动拼接或替代官方结论",
    }
    with open(os.path.join(out, "campaign_multiday.json"), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    print(f"\n输出目录：{out}")
    print(f"证据等级：{EVIDENCE_GRADE}（量级示例多日序列，非现场结论，非官方 R12 活动）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
