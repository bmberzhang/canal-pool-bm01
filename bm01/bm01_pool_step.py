# -*- coding: utf-8 -*-
"""
BM-01 上游阶跃—本池可吸收
================================================================
标准案例集 BM-01：对象 = 中线渠池；场景 = 上游阶跃—本池可吸收；预期结论码 = A1。

依据（口径来源）：
  自主水网标准案例集_BM-01至BM-10_文档脚本与物理原理_2026-09-09.md  第 5 节 BM-01

物理口径（与上述文档一致）：
  上游流量阶跃进入本池后，水位按调蓄关系上升：   Δh ≈ ΔQ · Δt / A
  闸门按本池水位偏差加大下泄，把多余水量送走。
  判定不是"完全无波动"，而是波动留在运行设计域内并收敛。

本脚本是"最小可运行案例"：合成降阶模型（synthetic-reduced-model），
参数为量级示例（assume），不是中线现场参数，不构成现场结论，不换算 L3 / 工程授权。

用法：
  python bm01_pool_step.py --scenario all --out .
  python bm01_pool_step.py --scenario all --out . --verify     # 复跑一致性核验
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from typing import Dict, List, Tuple

import numpy as np

CASE_ID = "BM-01"
CASE_TITLE = "上游阶跃—本池可吸收"
OBJECT = "中线渠池"
LOCAL_CAMPAIGN_ID = "local-bm01-synthetic-r1"
EVIDENCE_GRADE = "synthetic-reduced-model"
OFFICIAL_NOTE = (
    "本包为本地可复跑最小案例，非官方登记活动；"
    "官方 R12-standard-benchmarks-r1 仅完成题单结构契约检查（structure_checked），"
    "不等于物理回放已经完成。"
)


# ----------------------------------------------------------------------
# 参数：量级示例（assume）。全部为教材/规范量级的自洽取值，不是中线现场参数。
# ----------------------------------------------------------------------
@dataclass
class Params:
    # --- 渠池几何（中线典型渠池量级示例）---
    pool_len_m: float = 6000.0        # 池长
    pool_width_m: float = 50.0        # 平均水面宽
    # --- 水力工况 ---
    q0: float = 180.0                 # 稳态输水流量 m3/s
    h_target: float = 7.00            # 目标水位（相对渠底基准）m
    a0: float = 0.50                  # 稳态闸门开度（0=全关，1=全开）
    kg: float = 360.0                 # 闸门流量系数：Q = kg * a * sqrt(h/h_target)
    # --- 本池水位控制器（本地 PI）---
    kp: float = 1.60                  # 比例增益 [开度/m]
    ki: float = 1.60e-3               # 积分增益 [开度/(m·s)]
    ctrl_period_s: float = 300.0      # 控制周期（中线闸门水位控制为慢回路）
    deadband_m: float = 0.010         # 控制死区
    meas_resolution_m: float = 0.001  # 水位计量化分辨率
    rate_limit_per_s: float = 0.02 / 60.0   # 闸门启闭速率上限（满行程约 50 min）
    # --- 仿真 ---
    dt: float = 5.0
    duration_s: float = 14400.0       # 4 h
    t_step_s: float = 900.0           # 阶跃施加时刻
    step_dq: float = 30.0             # 上游流量阶跃幅度 m3/s
    # --- 判定条带（运行设计域）---
    band_normal_m: float = 0.15       # 正常调节条带 ±0.15 m
    band_alarm_m: float = 0.30        # 告警/吸收上限 ±0.30 m
    settle_tol_m: float = 0.020       # 收敛判据：末段偏差
    settle_window_s: float = 1800.0   # 末段窗口 30 min
    sat_tol: float = 1e-9             # 饱和判定容差
    sat_max_duration_s: float = 60.0  # 持续饱和容忍时长
    # --- 串联渠系 ---
    tau_chain_s: float = 300.0        # 上游闸门出流到下游池入流的传播滞后
    tau_ff_s: float = 60.0            # 前馈通道滤波时间常数
    ff_gain: float = 1.00             # 前馈增益（=1 表示按流量等价前馈）

    @property
    def area_m2(self) -> float:
        return self.pool_len_m * self.pool_width_m

    @property
    def gate_up_capacity(self) -> float:
        """开度从 a0 到全开可额外通过的流量（上游来水可吸收能力的执行侧上限）。"""
        return self.kg * (1.0 - self.a0)

    @property
    def gate_down_capacity(self) -> float:
        return self.kg * self.a0


# ----------------------------------------------------------------------
# 场景定义
# ----------------------------------------------------------------------
@dataclass
class Scenario:
    key: str
    label: str
    pools: int
    step_dq: float
    feedforward: bool
    note: str


def build_scenarios(base: Params) -> List[Scenario]:
    return [
        Scenario(
            key="main",
            label="主案例：单池 + 本地水位 PID",
            pools=1,
            step_dq=+45.0,
            feedforward=False,
            note=(
                "BM-01 本体。上游来水阶跃 +45 m³/s（约 +25%），"
                "闸门只用本池水位反馈，检验水位是否留在正常调节条带内并收敛。"
            ),
        ),
        Scenario(
            key="chain_local_only",
            label="对照一：两池串联 + 纯本地 PID（无前馈）",
            pools=2,
            step_dq=+72.0,
            feedforward=False,
            note=(
                "两池串联，两座闸门各自只看本池水位、互不预告。"
                "用于检验'单闸稳定不等于整链稳定'：上游池镇得住，下游池是否越幅。"
            ),
        ),
        Scenario(
            key="chain_feedforward",
            label="对照二：两池串联 + 渠段前馈加反馈",
            pools=2,
            step_dq=+72.0,
            feedforward=True,
            note=(
                "同一阶跃、同一对象，仅在下游闸门增加渠段前馈"
                "（按上游闸门流量变化等价前馈下游闸门开度）。用于验证越幅是否被压回条带内。"
            ),
        ),
        Scenario(
            key="overload",
            label="对照三：阶跃超出本池吸收能力",
            pools=1,
            step_dq=+200.0,
            feedforward=False,
            note=(
                "上游阶跃 +200 m³/s 已超过闸门余量（{cap:.0f} m³/s）。"
                "用于验证'超出吸收能力集就不再是控制问题'，应给出 A3 降级而非继续声称自动。"
            ),
        ),
    ]


# ----------------------------------------------------------------------
# 单池水力 + 本地控制器
# ----------------------------------------------------------------------
def gate_flow(p: Params, a: float, h: float) -> float:
    """闸门过流（平板闸孔自由出流的降阶形式，保留水位水头依赖）。"""
    head = max(h, 1e-3)
    return p.kg * a * math.sqrt(head / p.h_target)


def clamp(v: float, lo: float, hi: float) -> float:
    return lo if v < lo else (hi if v > hi else v)


class Pool:
    """渠池：蓄量平衡 dh/dt = (Qin - Qout)/A，闸门在池下游端。"""

    def __init__(self, p: Params):
        self.p = p
        self.h = p.h_target
        self.q_in = p.q0
        self.a = p.a0
        self.a_cmd = p.a0
        self.integral = 0.0
        self.last_ctrl_t = -1e18
        self.sat_time = 0.0
        self.rate_limited_time = 0.0
        self.sat_active = False

    # --- 本地水位 PI（含死区、条件积分抗饱和）---
    def control(self, t: float, feedforward_open: float = 0.0) -> None:
        p = self.p
        if t - self.last_ctrl_t < p.ctrl_period_s:
            return
        self.last_ctrl_t = t
        meas = round(self.h / p.meas_resolution_m) * p.meas_resolution_m
        e = meas - p.h_target
        if abs(e) >= p.deadband_m:
            # 条件积分：执行器已到边界且误差继续指向边界时不积分
            at_hi = self.a >= 1.0 - 1e-12 and e > 0
            at_lo = self.a <= 1e-12 and e < 0
            if not (at_hi or at_lo):
                self.integral += e * p.ctrl_period_s
        u = p.kp * e + p.ki * self.integral
        self.a_cmd = clamp(p.a0 + u + feedforward_open, 0.0, 1.0)

    def step_gate(self, dt: float) -> None:
        p = self.p
        step = p.rate_limit_per_s * dt
        delta = self.a_cmd - self.a
        if abs(delta) < step:
            self.a = self.a_cmd
        else:
            self.a += math.copysign(step, delta)
            self.rate_limited_time += dt
        self.a = clamp(self.a, 0.0, 1.0)

    def dynamics(self, dt: float) -> float:
        p = self.p
        q_out = gate_flow(p, self.a, self.h)
        self.h += (self.q_in - q_out) / p.area_m2 * dt
        return q_out

    def account_saturation(self, dt: float) -> None:
        if self.a <= 1e-9 or self.a >= 1.0 - 1e-9:
            self.sat_time += dt
            self.sat_active = True


# ----------------------------------------------------------------------
# 仿真
# ----------------------------------------------------------------------
def simulate(p: Params, sc: Scenario) -> Dict:
    n = int(round(p.duration_s / p.dt)) + 1
    t = np.arange(n) * p.dt

    pools = [Pool(p) for _ in range(sc.pools)]
    lag_state = p.q0          # 渠段传播一阶滞后状态（用于串联）
    ff_state = 0.0            # 前馈通道滤波状态

    series = {
        "t": t,
        "q_upstream": np.zeros(n),
        "h1_dev": np.zeros(n),
        "h2_dev": np.zeros(n),
        "a1": np.zeros(n),
        "a2": np.zeros(n),
        "q1_out": np.zeros(n),
        "q2_out": np.zeros(n),
    }

    for k in range(n):
        tk = float(t[k])

        # --- 上游阶跃注入 ---
        q_up = p.q0 + (sc.step_dq if tk >= p.t_step_s else 0.0)
        pools[0].q_in = q_up
        series["q_upstream"][k] = q_up

        # --- 前馈通道：按上游闸门流量变化等价折成下游闸门开度 ---
        if sc.pools == 2:
            q1_out_prev = series["q1_out"][k - 1] if k > 0 else p.q0
            ff_target = p.ff_gain * (q1_out_prev - p.q0) / p.kg if sc.feedforward else 0.0
            alpha = p.dt / max(p.tau_ff_s, p.dt)
            ff_state += alpha * (ff_target - ff_state)

        # --- 控制 ---
        if sc.pools == 1:
            pools[0].control(tk)
        else:
            pools[0].control(tk)
            pools[1].control(tk, feedforward_open=ff_state)

        # --- 串联入流：上游闸门出流经渠段传播滞后进入下游池 ---
        if sc.pools == 2:
            beta = p.dt / max(p.tau_chain_s, p.dt)
            lag_state += beta * (series["q1_out"][k - 1] - lag_state) if k > 0 else 0.0
            pools[1].q_in = lag_state

        # --- 水力推进 ---
        q1 = pools[0].dynamics(p.dt)
        series["q1_out"][k] = q1
        for pl in pools:
            pl.step_gate(p.dt)
            pl.account_saturation(p.dt)
        if sc.pools == 2:
            q2 = pools[1].dynamics(p.dt)
            series["q2_out"][k] = q2
        else:
            series["q2_out"][k] = q1

        series["h1_dev"][k] = pools[0].h - p.h_target
        series["a1"][k] = pools[0].a
        if sc.pools == 2:
            series["h2_dev"][k] = pools[1].h - p.h_target
            series["a2"][k] = pools[1].a

    return {
        "scenario": sc,
        "series": series,
        "pools": pools,
        "metrics": compute_metrics(p, sc, series),
    }


# ----------------------------------------------------------------------
# 判定
# ----------------------------------------------------------------------
def channel_metrics(p: Params, dev: np.ndarray, a: np.ndarray, dt: float) -> Dict:
    absdev = np.abs(dev)
    peak = float(absdev.max())
    t_peak = float(int(np.argmax(absdev)) * dt)
    tail = absdev[int(round((p.duration_s - p.settle_window_s) / dt)):]
    final_dev = float(tail.max())
    over = np.where(absdev > p.settle_tol_m)[0]
    settle_t = float(over[-1] * dt) if over.size else 0.0
    return {
        "peak_dev_m": round(peak, 5),
        "t_peak_s": t_peak,
        "final_window_max_dev_m": round(final_dev, 5),
        "settle_time_s": settle_t,
        "time_over_normal_band_s": float(np.sum(absdev > p.band_normal_m) * dt),
        "time_over_alarm_band_s": float(np.sum(absdev > p.band_alarm_m) * dt),
        "gate_min": round(float(a.min()), 5),
        "gate_max": round(float(a.max()), 5),
        "gate_final": round(float(a[-1]), 5),
    }


def judge(p: Params, sc: Scenario, r: Dict) -> Dict:
    """判据：R1 全程留条带内 / R2 闸门无持续饱和 / R3 末段收敛 / R4 信息闭合。"""
    m = r["metrics"]
    rules: List[Dict] = []
    chans = ["pool1"] + (["pool2"] if sc.pools == 2 else [])

    r1_ok = True
    for c in chans:
        cm = m[c]
        ok = cm["time_over_normal_band_s"] <= 0.0
        r1_ok &= ok
        rules.append({
            "rule": "R1_steady_band",
            "channel": c,
            "limit": p.band_normal_m,
            "actual": cm["peak_dev_m"],
            "pass": bool(ok),
            "desc": f"{c} 全程水位偏差不超过 ±{p.band_normal_m:.2f} m（正常调节条带）",
        })

    r2_ok = True
    for pl, c in zip(r["pools"], chans):
        ok = pl.sat_time <= p.sat_max_duration_s
        r2_ok &= ok
        rules.append({
            "rule": "R2_no_sustained_saturation",
            "channel": c,
            "limit": p.sat_max_duration_s,
            "actual": round(pl.sat_time, 1),
            "pass": bool(ok),
            "desc": "闸门无持续饱和（触及全开/全关的累计时长不超过容忍值）",
        })

    r3_ok = True
    for c in chans:
        cm = m[c]
        ok = cm["final_window_max_dev_m"] <= p.settle_tol_m
        r3_ok &= ok
        rules.append({
            "rule": "R3_settled",
            "channel": c,
            "limit": p.settle_tol_m,
            "actual": cm["final_window_max_dev_m"],
            "pass": bool(ok),
            "desc": f"{c} 末段 {p.settle_window_s/60:.0f} min 偏差不超过 ±{p.settle_tol_m*1000:.0f} mm（已收敛）",
        })

    rules.append({
        "rule": "R4_information_closure",
        "channel": "info",
        "limit": None,
        "actual": "无测点失效、无通信降级注入",
        "pass": True,
        "desc": "信息闭合：本案例不注入量测失效/通信退化，反馈回路完整",
    })

    all_pass = bool(r1_ok and r2_ok and r3_ok)
    exceed_alarm = any(m[c]["time_over_alarm_band_s"] > 0 for c in chans)

    if all_pass:
        code, verdict = "A1", "自动执行：扰动留在控制层，域内收敛，不必升级为调度事件"
    elif sc.step_dq > p.gate_up_capacity:
        code, verdict = "A3", "降级保守：阶跃超出本池可吸收能力，应收缩目标、邻段配合，不能继续声称自动"
    else:
        code, verdict = "A3", "降级保守：扰动越出运行设计域，应转降级而非硬顶"

    return {
        "code": code,
        "verdict": verdict,
        "rules": rules,
        "expected_code": "A1",
        "matches_expected": code == "A1",
        "exceeded_alarm_band": bool(exceed_alarm),
        "absorb_capacity_m3s": round(p.gate_up_capacity, 1),
        "step_over_capacity": bool(sc.step_dq > p.gate_up_capacity),
    }


def compute_metrics(p: Params, sc: Scenario, series: Dict) -> Dict:
    m = {"pool1": channel_metrics(p, series["h1_dev"], series["a1"], p.dt)}
    if sc.pools == 2:
        m["pool2"] = channel_metrics(p, series["h2_dev"], series["a2"], p.dt)
    return m


# ----------------------------------------------------------------------
# 输出与工件
# ----------------------------------------------------------------------
def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(p: Params, results: List[Dict], path: str) -> None:
    cols = ["scenario", "t", "q_upstream", "h1_dev", "a1", "q1_out", "h2_dev", "a2", "q2_out"]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(",".join(cols) + "\n")
        for r in results:
            s = r["series"]
            stride = max(1, int(round(60.0 / p.dt)))   # 每 60 s 一行，便于人读
            for k in range(0, len(s["t"]), stride):
                fh.write(",".join([
                    r["scenario"].key,
                    f"{s['t'][k]:.0f}",
                    f"{s['q_upstream'][k]:.2f}",
                    f"{s['h1_dev'][k]:.6f}",
                    f"{s['a1'][k]:.5f}",
                    f"{s['q1_out'][k]:.2f}",
                    f"{s['h2_dev'][k]:.6f}",
                    f"{s['a2'][k]:.5f}",
                    f"{s['q2_out'][k]:.2f}",
                ]) + "\n")


def make_figure(p: Params, results: List[Dict], path: str) -> None:
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

    by_key = {r["scenario"].key: r for r in results}
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.6))
    fig.suptitle(
        f"BM-01 上游阶跃—本池可吸收｜对象：{OBJECT}（量级示例·合成降阶模型）",
        fontsize=14, fontweight="bold",
    )

    def band(ax):
        ax.axhspan(-p.band_normal_m, p.band_normal_m, color="#e8f3ff", zorder=0,
                   label=f"正常条带 ±{p.band_normal_m:.2f} m")
        ax.axhline(p.band_alarm_m, color="#c0392b", ls=":", lw=1.1, zorder=1,
                   label=f"告警线 ±{p.band_alarm_m:.2f} m")
        ax.axhline(-p.band_alarm_m, color="#c0392b", ls=":", lw=1.1, zorder=1)
        ax.axvline(p.t_step_s / 60.0, color="#888", ls="--", lw=1.0, zorder=1,
                   label=f"阶跃时刻 {p.t_step_s/60:.0f} min")

    # (a) 主案例水位偏差
    ax = axes[0][0]
    r = by_key["main"]
    band(ax)
    ax.plot(r["series"]["t"] / 60, r["series"]["h1_dev"], color="#0066cc", lw=1.6,
            label=f"池1 水位偏差｜峰值 {r['metrics']['pool1']['peak_dev_m']*1000:.0f} mm"
                  f"｜末段 {r['metrics']['pool1']['final_window_max_dev_m']*1000:.1f} mm")
    ax.set_title(f"(a) 主案例：单池本地 PID，ΔQ=+{r['scenario'].step_dq:.0f} m³/s → 判定 "
                 f"{r['judgement']['code']}", fontsize=11)
    ax.set_xlabel("时间 (min)"); ax.set_ylabel("水位偏差 (m)")
    ax.legend(fontsize=8); ax.grid(alpha=.25)

    # (b) 主案例闸门开度
    ax = axes[0][1]
    ax.plot(r["series"]["t"] / 60, r["series"]["a1"], color="#1a7f37", lw=1.6, label="池1 闸门开度")
    ax.axhline(1.0, color="#c0392b", ls=":", lw=1.0)
    ax.axhline(p.a0, color="#888", ls="--", lw=1.0, label=f"初始开度 {p.a0:.2f}")
    ax.set_title(f"(b) 主案例闸门动作（余量 {p.gate_up_capacity:.0f} m³/s）", fontsize=11)
    ax.set_xlabel("时间 (min)"); ax.set_ylabel("开度 (0~1)")
    ax.legend(fontsize=8); ax.grid(alpha=.25)

    # (c) 两池串联对照：池2 水位偏差
    ax = axes[1][0]
    band(ax)
    dq_chain = by_key["chain_local_only"]["scenario"].step_dq
    for key, color, lab in [("chain_local_only", "#c0392b", "纯本地 PID（无前馈）"),
                            ("chain_feedforward", "#1a7f37", "渠段前馈 + 反馈")]:
        rr = by_key[key]
        ax.plot(rr["series"]["t"] / 60, rr["series"]["h2_dev"], color=color, lw=1.7,
                label=f"{lab}｜峰值 {rr['metrics']['pool2']['peak_dev_m']*1000:.0f} mm"
                      f"｜越带 {rr['metrics']['pool2']['time_over_normal_band_s']/60:.0f} min")
    ax.set_title(f"(c) 对照：两池串联时下游池水位偏差（同一 ΔQ=+{dq_chain:.0f} m³/s）", fontsize=11)
    ax.set_xlabel("时间 (min)"); ax.set_ylabel("池2 水位偏差 (m)")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=.25)

    # (d) 超出吸收能力对照
    ax = axes[1][1]
    rr = by_key["overload"]
    ax.plot(rr["series"]["t"] / 60, rr["series"]["h1_dev"], color="#8e44ad", lw=1.7,
            label=f"ΔQ=+{rr['scenario'].step_dq:.0f} m³/s（超能力）")
    ax2 = ax.twinx()
    ax2.plot(rr["series"]["t"] / 60, rr["series"]["a1"], color="#e67e22", lw=1.3, ls="-",
             label="闸门开度")
    ax2.set_ylabel("开度 (0~1)", color="#e67e22")
    ax2.set_ylim(-0.02, 1.05)
    ax.set_title(f"(d) 对照：阶跃超出可吸收能力 → 判定 {rr['judgement']['code']}", fontsize=11)
    ax.set_xlabel("时间 (min)"); ax.set_ylabel("池1 水位偏差 (m)")
    h1_, l1_ = ax.get_legend_handles_labels(); h2_, l2_ = ax2.get_legend_handles_labels()
    ax.legend(h1_ + h2_, l1_ + l2_, fontsize=8, loc="center right"); ax.grid(alpha=.25)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def make_summary(p: Params, results: List[Dict], path: str) -> None:
    L = []
    A = L.append
    A(f"# BM-01 判定摘要｜{OBJECT}｜{CASE_TITLE}\n")
    A(f"- 案例编号：`{CASE_ID}`　对象：{OBJECT}　场景：{CASE_TITLE}　预期结论码：**A1**")
    A(f"- 本地活动号：`{LOCAL_CAMPAIGN_ID}`（非官方登记）　证据等级：`{EVIDENCE_GRADE}`")
    A(f"- 生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    A("## 结论\n")
    A("| 场景 | 上游阶跃 | 判定码 | 与预期一致 | 池1 峰值偏差 | 池2 峰值偏差 | 收敛(末30min) | 结论 |")
    A("|---|---:|---|---|---:|---:|---:|---|")
    for r in results:
        sc, j, m = r["scenario"], r["judgement"], r["metrics"]
        p2 = f"{m['pool2']['peak_dev_m']:.3f} m" if "pool2" in m else "—"
        fin = m["pool1"]["final_window_max_dev_m"]
        A(f"| {sc.label} | +{sc.step_dq:.0f} m³/s | **{j['code']}** | "
          f"{'是' if j['matches_expected'] else '否（对照场景）'} | "
          f"{m['pool1']['peak_dev_m']:.3f} m | {p2} | {fin*1000:.1f} mm | {j['verdict']} |")
    A("")
    A("## 主案例判据逐条\n")
    main = next(r for r in results if r["scenario"].key == "main")
    A("| 规则 | 通道 | 限值 | 实测 | 结论 |")
    A("|---|---|---|---|---|")
    for rule in main["judgement"]["rules"]:
        lim = "—" if rule["limit"] is None else f"{rule['limit']}"
        A(f"| {rule['rule']} | {rule['channel']} | {lim} | {rule['actual']} | "
          f"{'通过' if rule['pass'] else '**不通过**'} |")
    A("")
    A("## 判据定义（本案例内生效）\n")
    A("| 规则 | 含义 | 限值 |")
    A("|---|---|---|")
    A(f"| R1 稳态条带 | 全程水位偏差不越出正常调节条带 | ±{p.band_normal_m:.2f} m |")
    A(f"| R2 无持续饱和 | 闸门触及全开/全关的累计时长不超过容忍值 | ≤{p.sat_max_duration_s:.0f} s |")
    A(f"| R3 已收敛 | 末段 {p.settle_window_s/60:.0f} min 内偏差不超过收敛容差 | ±{p.settle_tol_m*1000:.0f} mm |")
    A("| R4 信息闭合 | 无测点失效、无通信降级注入，反馈回路完整 | 布尔 |")
    A("")
    A("**判定合成**：R1∧R2∧R3∧R4 全部满足 → **A1 自动执行**；任一项不满足 → **A3 降级保守**"
      f"（若阶跃本身已超过闸门余量 {p.gate_up_capacity:.0f} m³/s，则属于能力不足，必须降级而非硬顶）。")
    A(f"告警线 ±{p.band_alarm_m:.2f} m 用于标记越限风险等级，不参与 A1/A3 合成。\n")
    A("## 对照解读（本题的核心论据）\n")
    chain_local = next(r for r in results if r["scenario"].key == "chain_local_only")
    chain_ff = next(r for r in results if r["scenario"].key == "chain_feedforward")
    ov = next(r for r in results if r["scenario"].key == "overload")
    A(f"1. **单池可吸收成立**：ΔQ=+{main['scenario'].step_dq:.0f} m³/s 时，池1 峰值偏差 "
      f"{main['metrics']['pool1']['peak_dev_m']*1000:.0f} mm，仅为条带的 "
      f"{main['metrics']['pool1']['peak_dev_m']/p.band_normal_m*100:.0f}%，末段收敛到 "
      f"{main['metrics']['pool1']['final_window_max_dev_m']*1000:.1f} mm → A1。")
    A(f"2. **单闸稳定 ≠ 整链稳定**：同样的两座闸、同样的本地 PID，ΔQ=+"
      f"{chain_local['scenario'].step_dq:.0f} m³/s 时上游池峰值 "
      f"{chain_local['metrics']['pool1']['peak_dev_m']*1000:.0f} mm 仍留在条带内，"
      f"但下游池峰值 {chain_local['metrics']['pool2']['peak_dev_m']*1000:.0f} mm、越带 "
      f"{chain_local['metrics']['pool2']['time_over_normal_band_s']/60:.0f} min → A3。"
      "每一池都在盯自己的水位，但没有任何一环预知上游多来了多少水。")
    A(f"3. **渠段前馈把下游越幅压下去**：同一阶跃下，下游闸门增加按上游闸门流量变化等价的"
      f"前馈后，下游池峰值由 {chain_local['metrics']['pool2']['peak_dev_m']*1000:.0f} mm 降到 "
      f"{chain_ff['metrics']['pool2']['peak_dev_m']*1000:.0f} mm、越带由 "
      f"{chain_local['metrics']['pool2']['time_over_normal_band_s']/60:.0f} min 降到 "
      f"{chain_ff['metrics']['pool2']['time_over_normal_band_s']/60:.0f} min → A1。")
    A(f"4. **能力不够时必须降级**：ΔQ=+{ov['scenario'].step_dq:.0f} m³/s 已超过闸门向上余量 "
      f"{p.gate_up_capacity:.0f} m³/s，闸门很快到全开并累计饱和 "
      f"{ov['pools'][0].sat_time/60:.0f} min，水位持续单向上抬至 "
      f"{ov['metrics']['pool1']['peak_dev_m']*1000:.0f} mm 且末段仍为 "
      f"{ov['metrics']['pool1']['final_window_max_dev_m']*1000:.0f} mm（不收敛）→ A3。"
      "此类情形不是控制参数问题，而是能力集问题。")
    A("")
    A("## 物理口径\n")
    A("上游流量阶跃进入本池后，水位按调蓄关系上升：\n")
    A("```\nΔh ≈ ΔQ · Δt / A\n\n")
    A(f"A = 池长 × 平均水面宽 = {p.pool_len_m:.0f} m × {p.pool_width_m:.0f} m = {p.area_m2:.3e} m²")
    dq_main = main["scenario"].step_dq
    A(f"ΔQ = {dq_main:.0f} m³/s → 无控制自由上升速率 ΔQ/A = {dq_main/p.area_m2:.3e} m/s "
      f"= {dq_main/p.area_m2*3600:.3f} m/h\n```\n")
    A(f"闸门按本池水位偏差加大下泄：`Q_gate = {p.kg:.0f} · a · √(h/h_target)`，"
      f"初始开度 {p.a0:.2f}，向上余量 {p.gate_up_capacity:.0f} m³/s、向下余量 {p.gate_down_capacity:.0f} m³/s。")
    A("判定不是“完全无波动”，而是波动留在运行设计域内并收敛。\n")
    A("## 边界声明（不得外推）\n")
    A("1. 参数为量级示例（assume）：池长、水面宽、流量、闸门系数、PID 增益、条带宽度均为自洽取值，"
      "**不是中线现场参数**，不构成现场结论。")
    A("2. 证据等级为 `synthetic-reduced-model`：降阶蓄量模型 + 线性化闸门，"
      "不是特征线法/圣维南方程高保真回放，不等同 XiL、SIL 实物或工程授权。")
    A("3. 本包是本地最小案例，**不是**官方活动 `R12-standard-benchmarks-r1` 的一部分；"
      "官方活动当前仅完成题单结构契约检查。")
    A("4. 判定码 A1/A3 只在本案例的条带与判据定义下成立，不得换算 L3、WNAL 或工程定案；"
      "不得与其他轨道（WNTS 八对象、DB、桌面、中线专项、大渡河）合并计分。")
    A(f"5. {OFFICIAL_NOTE}")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


def write_records(results: List[Dict], path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for r in results:
            sc, j, m = r["scenario"], r["judgement"], r["metrics"]
            rec = {
                "schema": "awn-bm-record/v1",
                "case_id": CASE_ID,
                "case_title": CASE_TITLE,
                "object": OBJECT,
                "scene": sc.key,
                "scene_label": sc.label,
                "injection": {"type": "upstream_flow_step", "delta_q_m3s": sc.step_dq,
                              "t_step_s": 900.0},
                "controller": "local_pool_level_PI" + ("_plus_chain_feedforward" if sc.feedforward else ""),
                "evidence_grade": EVIDENCE_GRADE,
                "verdict_code": j["code"],
                "expected_code": j["expected_code"],
                "matches_expected": j["matches_expected"],
                "metrics": m,
                "rules": j["rules"],
                "official_activity": None,
                "note": OFFICIAL_NOTE,
            }
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")


def write_campaign(out_dir: str, cmd: str, artifacts: List[str], started: float,
                   finished: float) -> str:
    path = os.path.join(out_dir, "campaign.json")
    payload = {
        "campaign_id": LOCAL_CAMPAIGN_ID,
        "case_id": CASE_ID,
        "official_activity": None,
        "note": OFFICIAL_NOTE,
        "registered_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host": {"platform": platform.platform(), "python": sys.version.split()[0]},
        "commands": [{
            "cmd": cmd,
            "exit_code": 0,
            "started_epoch": round(started, 3),
            "finished_epoch": round(finished, 3),
            "duration_s": round(finished - started, 3),
        }],
        "artifacts": [
            {"path": os.path.relpath(a, out_dir).replace("\\", "/"),
             "sha256": sha256_file(a),
             "size_bytes": os.path.getsize(a)}
            for a in artifacts
        ],
        "sealed": False,
        "seal_note": "本地最小案例，未走官方封存流程；不得与官方活动拼接或替代官方结论",
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    return path


# ----------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="BM-01 上游阶跃—本池可吸收（本地最小可运行案例）")
    ap.add_argument("--scenario", default="all",
                    choices=["all", "main", "chain_local_only", "chain_feedforward", "overload"])
    ap.add_argument("--out", default=".")
    ap.add_argument("--verify", action="store_true", help="复跑一致性核验并写 review/reproducibility.json")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    for sub in ("raw", "normalized", "figures", "review"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    started = time.time()
    p = Params()
    all_sc = build_scenarios(p)
    scs = all_sc if args.scenario == "all" else [s for s in all_sc if s.key == args.scenario]

    results = []
    for sc in scs:
        r = simulate(p, sc)
        r["judgement"] = judge(p, sc, r)
        results.append(r)

    # ---- raw json（逐场景）----
    artifacts: List[str] = []
    for r in results:
        sc, s = r["scenario"], r["series"]
        stride = max(1, int(round(30.0 / p.dt)))
        payload = {
            "case_id": CASE_ID,
            "object": OBJECT,
            "scenario": asdict(sc),
            "params": {**asdict(p), "area_m2": p.area_m2,
                       "gate_up_capacity_m3s": p.gate_up_capacity,
                       "gate_down_capacity_m3s": p.gate_down_capacity},
            "judgement": r["judgement"],
            "metrics": r["metrics"],
            "series_downsampled": {
                "dt_s": p.dt * stride,
                "t_s": [round(float(v), 1) for v in s["t"][::stride]],
                "q_upstream_m3s": [round(float(v), 3) for v in s["q_upstream"][::stride]],
                "pool1_dev_m": [round(float(v), 6) for v in s["h1_dev"][::stride]],
                "pool1_gate": [round(float(v), 6) for v in s["a1"][::stride]],
                "pool1_out_m3s": [round(float(v), 3) for v in s["q1_out"][::stride]],
                "pool2_dev_m": [round(float(v), 6) for v in s["h2_dev"][::stride]],
                "pool2_gate": [round(float(v), 6) for v in s["a2"][::stride]],
            },
            "evidence_grade": EVIDENCE_GRADE,
            "note": OFFICIAL_NOTE,
        }
        fp = os.path.join(out, "raw", f"bm01_{sc.key}.json")
        with open(fp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2)
        artifacts.append(fp)
        print(f"[{sc.key}] 判定 {r['judgement']['code']}  "
              f"池1峰值={r['metrics']['pool1']['peak_dev_m']:.4f} m  "
              f"池1末段={r['metrics']['pool1']['final_window_max_dev_m']*1000:.1f} mm"
              + (f"  池2峰值={r['metrics']['pool2']['peak_dev_m']:.4f} m"
                 f"  池2末段={r['metrics']['pool2']['final_window_max_dev_m']*1000:.1f} mm"
                 if "pool2" in r["metrics"] else ""))

    csv_p = os.path.join(out, "raw", "bm01_timeseries.csv")
    write_csv(p, results, csv_p); artifacts.append(csv_p)

    fig_p = os.path.join(out, "figures", "bm01_level_response.png")
    make_figure(p, results, fig_p); artifacts.append(fig_p)

    rec_p = os.path.join(out, "normalized", "records.jsonl")
    write_records(results, rec_p); artifacts.append(rec_p)

    sum_p = os.path.join(out, "summary.md")
    make_summary(p, results, sum_p); artifacts.append(sum_p)

    # ---- 复跑一致性核验 ----
    if args.verify:
        again = []
        for sc in scs:
            r2 = simulate(p, sc)
            again.append(r2["metrics"])
        same = all(
            json.dumps(a, sort_keys=True) == json.dumps(b["metrics"], sort_keys=True)
            for a, b in zip(again, results)
        )
        rep = {
            "check": "deterministic_rerun",
            "pairs": f"{len(again)}/{len(results)}",
            "identical_metrics": bool(same),
            "mismatches": [] if same else ["metrics 不一致"],
            "note": "同一输入下复跑指标逐字段一致，用于证明本案例可复跑；不增加任何分值。",
        }
        rp = os.path.join(out, "review", "reproducibility.json")
        with open(rp, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, ensure_ascii=False, indent=2)
        artifacts.append(rp)
        print(f"复跑一致性：{len(again)}/{len(results)} 配对，identical={same}")

    finished = time.time()
    cmd = "python bm01_pool_step.py " + " ".join(sys.argv[1:])
    write_campaign(out, cmd, artifacts, started, finished)

    print(f"\n输出目录：{out}")
    print(f"证据等级：{EVIDENCE_GRADE}（量级示例，非现场结论，非官方 R12 活动）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
