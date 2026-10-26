# -*- coding: utf-8 -*-
"""
BM-01 一维圣维南模型 —— 精简核心版（约 120 行）
==================================================
只保留"核心求解器 sve_1d + 最少依赖"，去掉打印报表、开关包装、命令行参数。
适合贴进论文附录或当教学示例。跑法和结果与完整版完全一致。

解：一维圣维南（浅水）方程，守恒形式 + Rusanov 通量 + 显式欧拉
    ∂A/∂t + ∂Q/∂x = 0
    ∂Q/∂t + ∂(Q²/A + gA²/2B)/∂x = gA(S0 − Sf)

只用 Python 标准库，无 numpy / matplotlib。
"""
import math

G = 9.81
MANNING_N = 0.014
CFL = 0.40
DT_CANDIDATES = [5, 4, 3.6, 3, 2.5, 2.25, 2, 1.8, 1.5, 1.25, 1.125, 1,
                 0.9, 0.75, 0.6, 0.5, 0.45, 0.4, 0.36, 0.3, 0.25, 0.225, 0.2]


class Section:
    def __init__(self, b, m=0.0):
        self.b, self.m = b, m
    def A(self, h): return h * (self.b + self.m * h)
    def B(self, h): return self.b + 2.0 * self.m * h
    def P(self, h): return self.b + 2.0 * h * math.sqrt(1.0 + self.m * self.m)


class Ctrl:
    def __init__(self, p):
        self.p = p
        self.a = p["a0"]; self.a_cmd = p["a0"]
        self.integral = 0.0; self.last_ctrl_t = -1e18
    def control(self, t, h):
        p = self.p
        if t - self.last_ctrl_t < p["ctrl_period_s"]: return
        self.last_ctrl_t = t
        meas = round(h / p["meas_resolution_m"]) * p["meas_resolution_m"]
        e = meas - p["h_target"]
        if abs(e) >= p["deadband_m"] and not (
            self.a >= 1.0 - 1e-12 and e > 0 or self.a <= 1e-12 and e < 0):
            self.integral += e * p["ctrl_period_s"]
        u = p["kp"] * e + p["ki"] * self.integral
        self.a_cmd = min(max(p["a0"] + u, 0.0), 1.0)
    def step_gate(self, dt):
        step = self.p["rate_limit_per_s"] * dt
        d = self.a_cmd - self.a
        if abs(d) < step: self.a = self.a_cmd
        else: self.a += step if d > 0 else -step
        self.a = min(max(self.a, 0.0), 1.0)


def normal_slope(p, n=MANNING_N):
    B, ht = p["pool_width_m"], p["h_target"]
    A = B * ht; R = A / (B + 2.0 * ht)
    return (n * p["q0"] / (A * R ** (2.0 / 3.0))) ** 2


def section_normal_slope(sec, h, Q, n):
    A = sec.A(h); P = sec.P(h); R = A / P if P > 0 else h
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


def pick_dt(dx, wmax):
    lo = CFL * dx / max(wmax, 1e-6)
    for c in DT_CANDIDATES:
        if c <= lo: return c
    raise ValueError("CFL too tight")


def sve_1d(p, ncell, step_dq, *, ctrl_on="gate_avg", m_slope=0.0,
           n=MANNING_N, s0=None, sensor_len_m=100.0, early_window_s=7200.0):
    B0, L, ht = p["pool_width_m"], p["pool_len_m"], p["h_target"]
    sec = Section(B0, m_slope)
    dx = L / ncell
    S0 = normal_slope(p, n) if s0 is None else s0
    h = [ht] * ncell; Q = [p["q0"]] * ncell
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))
    wmax = p["q0"] / (B0 * ht) + math.sqrt(G * ht)
    dt = pick_dt(dx, wmax)
    nt = int(round(p["duration_s"] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = Ctrl(p)
    sqrt = math.sqrt; hmin = 1e-3; c1 = dt / dx
    friction_correct = m_slope > 1e-12
    out = {"t": [], "dev_mean": [], "dev_gate": [], "gate": []}
    peak_mean = peak_gate = early_peak = tail_mean = 0.0
    t_tail = p["duration_s"] - p["settle_window_s"]
    over_band_dt = 0; vol_net = 0.0

    for k in range(nt):
        tk = k * dt
        q_up = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)
        q_gate = p["kg"] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)
        if ctrl_on != "off":
            if ctrl_on == "gate": meas = h[-1]
            elif ctrl_on == "gate_avg": meas = sum(h[ncell - nsens:]) / nsens
            else: meas = sum(h) / ncell
            ctrl.control(tk, meas); ctrl.step_gate(dt)

        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif; F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = max(hx[j], hmin), max(hx[j + 1], hmin)
            A1, A2 = max(sec.A(h1), 1e-6), max(sec.A(h2), 1e-6)
            B1, B2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / B1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / B2)
            sm = max(s1, s2)
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2
                           + 0.5 * G * A1 * A1 / B1 + 0.5 * G * A2 * A2 / B2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate

        h_new = [0.0] * ncell; Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]; A = max(sec.A(hi), 1e-6); Qi = Q[i]
            R = A / sec.P(hi) if friction_correct and sec.P(hi) > 0 else A / (B0 + 2.0 * hi)
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        vol_net += (q_up - q_gate) * dt

        dm = sum(h) / ncell - ht; dg = h[-1] - ht
        peak_mean = max(peak_mean, abs(dm)); peak_gate = max(peak_gate, abs(dg))
        if tk <= early_window_s: early_peak = max(early_peak, abs(dm))
        if tk >= t_tail: tail_mean = max(tail_mean, abs(dm))
        if abs(dm) > p["band_normal_m"]: over_band_dt += 1
        if k % rec_every == 0:
            out["t"].append(tk); out["dev_mean"].append(dm)
            out["dev_gate"].append(dg); out["gate"].append(ctrl.a)

    dv = sum((x - ht) * (sec.b + sec.m * (x + ht)) * dx for x in h)
    return {"peak_mean_mm": peak_mean * 1000, "peak_gate_mm": peak_gate * 1000,
            "early_peak_mm": early_peak * 1000, "final_mean_mm": tail_mean * 1000,
            "over_band_s": over_band_dt * dt, "dt_s": dt,
            "mass_balance_rel_err": abs(dv - vol_net) / max(abs(vol_net), 1.0),
            "series": out}


if __name__ == "__main__":
    P = {
        "pool_len_m": 6000.0, "pool_width_m": 50.0, "q0": 180.0,
        "h_target": 7.0, "a0": 0.5, "kg": 360.0, "kp": 1.6, "ki": 0.0016,
        "ctrl_period_s": 300.0, "deadband_m": 0.01, "meas_resolution_m": 0.001,
        "rate_limit_per_s": 1.0 / 3000.0, "duration_s": 14400.0,
        "t_step_s": 900.0, "band_normal_m": 0.15, "band_alarm_m": 0.3,
        "settle_tol_m": 0.02, "settle_window_s": 1800.0,
    }
    DQ = 45.0
    print("矩形基准（闸前 100 m 平均传感器）")
    for nc in (240, 480):
        r = sve_1d(P, nc, DQ, ctrl_on="gate_avg")
        print(f"  {nc} 格：前2h峰 {r['early_peak_mm']:.1f} mm｜"
              f"池均峰 {r['peak_mean_mm']:.1f} mm｜守恒 {r['mass_balance_rel_err']:.1e}")
    print("梯形开关（边坡 m）")
    for m in (1.0, 1.5, 2.5):
        sec = Section(P["pool_width_m"], m)
        s0 = section_normal_slope(sec, P["h_target"], P["q0"], MANNING_N)
        r = sve_1d(P, 480, DQ, ctrl_on="gate_avg", m_slope=m, s0=s0)
        print(f"  m={m}：前2h峰 {r['early_peak_mm']:.1f} mm")
