# -*- coding: utf-8 -*-
"""
(E-FD) 频域机制分离 + 环路增益**直接测量**
================================================================================
论文最大的逻辑漏洞：因果链「取样位置局部化 → 环路增益升高 → 失稳 → 极限环」
里，中间环节「环路增益升高」**从未被直接测量**，只有端点证据。
本实验用频率响应把它测出来。

Part A —— 极限环周期辨识（判定"振荡由谁设定"）
  对闭环保留三档取样（闸前单点 / 闸前100m平均 / 全池平均）跑 12 h，取尾段 dev_mean，
  用自相关 + 过零双估计器求周期，与三个候选时间尺度比较：
    池体固有：T_half = 2L/c、T_quarter = 4L/c（c = √(g·h) + v）
    控制器：ctrl_period_s = 300 s；闸门全行程 = 0.5 / rate_limit = 1500 s
  判据：若 T_LC 既不等于 2L/c 也不等于 4L/c，而与控制器特征时间同一量级 → 振荡由
        控制器非线性（限幅/速率）设定，而非池体水力学固有模态 → "非水力学不可吸收"。

Part B —— 开环频率响应（直接测"取样位置 → 增益"）
  断开控制器（闸门不再反馈），在闸门开度上叠加小正弦 a = a0 + dA·sin(2πt/T)，
  用 Goertzel 在稳态段测三个传感通道（池均 / 闸前100m平均 / 闸前单点）在同频处的振幅。
  → 直接给出「取样窗 → 测量通道增益」的频响曲线。
  预期：低频（准稳态）三通道增益一致（空间均匀）；高频时闸前局部通道增益显著高于池均
  → 这就是"取样位置局部化抬高环路增益"的**直接证据**。

输出：raw/bm01_fd_mech.json
"""
import json, os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2
from bm01_sve1d_core import Section, G, MANNING_N, pick_dt, normal_slope

L = 6000.0
B = 50.0
HT = 7.0
Q0 = 180.0
KG = 360.0
A0 = 0.5
BAND_MM = 150.0
SENSOR_LEN_M = 100.0


def closed_loop_periods():
    """Part A：三档取样，12 h 闭环，辨识极限环周期。"""
    P = {"pool_len_m": L, "pool_width_m": B, "q0": Q0, "h_target": HT, "a0": A0,
         "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": 300.0, "deadband_m": 0.01,
         "meas_resolution_m": 0.001, "rate_limit_per_s": 1.0 / 3000.0, "t_step_s": 900.0,
         "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
         "settle_window_s": 1800.0, "duration_s": 43200.0}
    arms = (("point", dict(ctrl_on="gate")),
            ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=SENSOR_LEN_M)),
            ("mean", dict(ctrl_on="mean")))
    out = {}
    nc = 480
    for name, kw in arms:
        t1 = time.time()
        r = V2.sve_1d_v2(P, nc, 45.0, m_slope=0.0, **kw)
        t = r["series"]["t"]
        y = [v * 1000.0 for v in r["series"]["dev_mean"]]  # mm
        per = estimate_period(t, y, t_start=21600.0)
        out[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                     "final_mean_mm": round(r["final_mean_mm"], 2),
                     "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3",
                     "period_s": per, "dt_s": r["dt_s"],
                     "tail_amp_mm": round(0.5 * (max(y[-len(y) // 8:] or [0]) - min(y[-len(y) // 8:] or [0])), 2),
                     "wall_s": round(time.time() - t1, 1)}
        print("[FD-A] %-7s peak=%8.1f mm %s  T=%s s" %
              (name, out[name]["peak_mean_mm"], out[name]["judgement"], per), flush=True)
    return out


def _acf(y, maxlag):
    n = len(y)
    m = sum(y) / n
    yc = [v - m for v in y]
    den = sum(v * v for v in yc)
    if den <= 0:
        return [0.0] * maxlag
    ac = [0.0] * maxlag
    for k in range(1, maxlag):
        s = 0.0
        for i in range(n - k):
            s += yc[i] * yc[i + k]
        ac[k] = s / den
    return ac


def estimate_period(t, y, t_start):
    """自相关第一显著峰 + 上穿零间隔，双估计。"""
    tt = [t[i] for i in range(len(t)) if t[i] >= t_start]
    yy = [y[i] for i in range(len(t)) if t[i] >= t_start]
    n = len(yy)
    if n < 60:
        return None
    dt = tt[1] - tt[0]
    # --- 过零（上穿）---
    m = sum(yy) / n
    zc = []
    for i in range(1, n):
        if yy[i - 1] - m <= 0.0 < yy[i] - m:
            zc.append(tt[i])
    T_zc = None
    if len(zc) >= 3:
        diffs = [zc[i + 1] - zc[i] for i in range(len(zc) - 1)]
        if diffs:
            T_zc = round(sum(diffs) / len(diffs), 1)
    # --- 自相关第一峰 ---
    maxlag = min(int(n / 2), 3000)
    ac = _acf(yy, maxlag)
    T_ac = None
    for k in range(2, maxlag - 1):
        if ac[k] > ac[k - 1] and ac[k] >= ac[k + 1] and ac[k] > 0.15:
            T_ac = round(k * dt, 1)
            break
    swing = max(yy) - min(yy)
    return {"T_zerocross_s": T_zc, "T_acf_s": T_ac, "tail_swing_mm": round(swing, 2),
            "n_zerocross": len(zc)}


def probe_open_loop(ncell, period_s, dA=0.010, n_periods=6.0, analyze_periods=2.0,
                    m_slope=0.0, sensor_len_m=SENSOR_LEN_M):
    """Part B：开环闸门正弦注入 → 三通道频响（Goertzel）。"""
    sec = Section(B, m_slope)
    dx = L / ncell
    S0 = normal_slope({"pool_width_m": B, "h_target": HT, "q0": Q0}, MANNING_N)
    h = [HT] * ncell
    Q = [Q0] * ncell
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))
    wmax = Q0 / (B * HT) + math.sqrt(G * HT)
    dt = pick_dt(dx, wmax)
    nt = int(round(n_periods * period_s / dt)) + 1
    t_an = (n_periods - analyze_periods) * period_s
    w = 2.0 * math.pi / period_s
    fric = m_slope > 1e-12
    hmin = 1e-3
    c1 = dt / dx
    ch = {"mean": [0.0, 0.0, 0], "avg100": [0.0, 0.0, 0], "point": [0.0, 0.0, 0]}
    sqrt = math.sqrt
    cos, sin = math.cos, math.sin

    for k in range(nt):
        tk = k * dt
        q_up = Q0
        a_act = A0 + dA * sin(w * tk)
        q_gate = KG * a_act * sqrt(max(h[-1], hmin)) / sqrt(HT)

        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif
        F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = max(hx[j], hmin), max(hx[j + 1], hmin)
            A1 = max(sec.A(h1), 1e-6)
            A2 = max(sec.A(h2), 1e-6)
            B1, B2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / B1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / B2)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2
                           + 0.5 * G * A1 * A1 / B1 + 0.5 * G * A2 * A2 / B2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate

        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            Aa = max(sec.A(hi), 1e-6)
            Qi = Q[i]
            R = Aa / sec.P(hi) if fric and sec.P(hi) > 0 else Aa / (B + 2.0 * hi)
            Sf = MANNING_N * MANNING_N * Qi * abs(Qi) / (Aa * Aa * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * Aa * (S0 - Sf)
        h, Q = h_new, Q_new

        if tk >= t_an:
            ph = w * tk
            cc, ss = cos(ph), sin(ph)
            dm = sum(h) / ncell - HT
            dg100 = sum(h[ncell - nsens:]) / nsens - HT
            dp = h[-1] - HT
            for key, val in (("mean", dm), ("avg100", dg100), ("point", dp)):
                ch[key][0] += val * cc
                ch[key][1] += val * ss
                ch[key][2] += 1

    res = {}
    for key, (cs, sn, cnt) in ch.items():
        if cnt == 0:
            res[key] = {"amp_mm": None, "phase_deg": None}
            continue
        amp = 2.0 / cnt * math.hypot(cs, sn) * 1000.0   # mm
        ph = math.degrees(math.atan2(-sn, cs))
        res[key] = {"amp_mm": round(amp, 4), "phase_deg": round(ph, 1)}
    return {"ncell": ncell, "dx_m": dx, "dt_s": dt, "steps": nt,
            "period_s": period_s, "dA": dA, "n_periods": n_periods, "channels": res}


def main():
    t0 = time.time()
    print("=== Part A：极限环周期辨识（闭环 12h, 480 格） ===", flush=True)
    A = closed_loop_periods()

    c = math.sqrt(G * HT) + Q0 / (B * HT)
    scales = {
        "celerity_m_s": round(c, 3),
        "T_half_pool_2L_c_s": round(2.0 * L / c, 1),
        "T_quarter_pool_4L_c_s": round(4.0 * L / c, 1),
        "ctrl_period_s": 300.0,
        "gate_full_stroke_s": round(0.5 / (1.0 / 3000.0), 1),
    }

    print("\n=== Part B：开环闸门正弦注入 → 三通道频响（240 格, dA=1%" + "） ===", flush=True)
    periods = [1200.0, 2400.0, 4800.0, 9600.0, 14400.0]
    B_res = {}
    for T in periods:
        t1 = time.time()
        r = probe_open_loop(240, T)
        B_res["T%.0f" % T] = r
        cm = r["channels"]
        ratio_p = cm["point"]["amp_mm"] / max(cm["mean"]["amp_mm"], 1e-9)
        ratio_a = cm["avg100"]["amp_mm"] / max(cm["mean"]["amp_mm"], 1e-9)
        print("[FD-B] T=%6.0fs  amp: mean=%8.3f  avg100=%8.3f  point=%8.3f mm | "
              "ratio(100m/mean)=%6.2f  ratio(point/mean)=%7.2f (%.0fs)"
              % (T, cm["mean"]["amp_mm"], cm["avg100"]["amp_mm"], cm["point"]["amp_mm"],
                 ratio_a, ratio_p, time.time() - t1), flush=True)

    out = {"purpose": "E-FD 频域机制分离 + 环路增益直接测量",
           "config": {"pool_len_m": L, "pool_width_m": B, "h_target": HT, "q0": Q0,
                      "kg": KG, "a0": A0, "sensor_len_m": SENSOR_LEN_M},
           "partA_limit_cycle": A, "reference_time_scales": scales,
           "partB_frequency_response": B_res,
           "wall_total_s": round(time.time() - t0, 1)}

    # 结论摘要
    ratios = {}
    for T in periods:
        cm = B_res["T%.0f" % T]["channels"]
        ratios["T%.0f" % T] = {
            "avg100_over_mean": round(cm["avg100"]["amp_mm"] / max(cm["mean"]["amp_mm"], 1e-9), 3),
            "point_over_mean": round(cm["point"]["amp_mm"] / max(cm["mean"]["amp_mm"], 1e-9), 3)}
    out["gain_ratio_summary"] = ratios

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_fd_mech.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))
    print("参考时间尺度:", json.dumps(scales, ensure_ascii=False))
    print("增益比:", json.dumps(ratios, ensure_ascii=False))


if __name__ == "__main__":
    main()
