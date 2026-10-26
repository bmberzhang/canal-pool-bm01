# -*- coding: utf-8 -*-
"""
任务 A（双模型共识最高优先）：反饱和四臂 12h → 24h 全臂延伸 + 五关判决
================================================================================
目的：E4 已证"闸门限幅+积分饱和是维持极限环的必要条件"（12h 窗口内的观测）。
本实验把同样的四臂延伸到 24h，回答：这条"必要条件"是**时间尺度稳健**的结论，
还是 12h 窗口内的暂态假象？

四臂（m=0 矩形 / 480 格 / 采样窗长 1500 m / 阶跃 +45 m³/s）：
  base      a∈[0.00, 1.00]  现有 PI（条件积分抗饱和）
  tightclip a∈[0.05, 0.95]  收紧物理限幅
  backcalc  a∈[0.00, 1.00]  back-calculation 抗饱和（未饱和命令才积分）
  noclip    a∈[-0.20, 1.20] 诊断臂：拿掉限幅非线性

五关判决（先写死，防事后拟合）：
  ① 叠合预检：24h 运行的前 12h 与既有 12h 结果必须一致（重合 ≥50%），否则本实验无效；
  ② 包络中线零趋势：对末 12h 各周期包络中线 ((max+min)/2) 做线性回归，
     斜率 95% CI 必须跨 0，且 CI 宽度 ≤0.10 mm/min 才继续（否则按方向判"真失败/衰减"）；
  ③ 周期数 N_cycle ≥ 3（末 12h）；
  ④ 峰峰变异系数 CV_p2p ≤ 0.4（>0.4 → 弱通过最低档 + base 延 36h，仍 >0.4 判真失败）；
  ⑤ μ_p2p 与 k_p2p 的 95% CI → 判"强通过 / 弱通过 / 真失败"。

输出：raw/bm01_a_24h_antiwindup.json
"""
import os, sys, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2
import bm01_sve1d_compare as V1

G = V1.G; MANNING_N = V1.MANNING_N; CFL = V1.CFL
LIMITS = {'base': (0.0, 1.0), 'tightclip': (0.05, 0.95), 'backcalc': (0.0, 1.0), 'noclip': (-0.2, 1.2)}
WIN_M = 1500.0
DUR_24H = 86400.0
DUR_12H = 43200.0
STEP_DQ = 45.0

# 既有 E4 的 12h 结果（用于叠合预检）
E4_12H = {'base': 1374.0, 'tightclip': 417.5, 'backcalc': 253.1, 'noclip': 2309.0}

T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
       8: 2.306, 9: 2.262, 10: 2.228, 12: 2.179, 15: 2.131, 20: 2.086, 30: 2.042}


def tq(df):
    if df <= 0: return float('inf')
    for d in sorted(T95):
        if df <= d: return T95[d]
    return 1.96


class CtrlX:
    def __init__(self, p, mode='base'):
        self.p = p; self.mode = mode; self.a = p['a0']; self.a_cmd = p['a0']
        self.integral = 0.0; self.last_ctrl_t = -1e18

    def control(self, t, h):
        p = self.p
        if t - self.last_ctrl_t < p['ctrl_period_s']:
            return
        self.last_ctrl_t = t
        meas = round(h / p['meas_resolution_m']) * p['meas_resolution_m']
        e = meas - p['h_target']
        if self.mode == 'backcalc':
            a_uns = p['a0'] + p['kp'] * e + p['ki'] * self.integral
            if 0.0 < a_uns < 1.0 and abs(e) >= p['deadband_m']:
                self.integral += e * p['ctrl_period_s']
        else:
            if abs(e) >= p['deadband_m']:
                at_hi = self.a >= 1.0 - 1e-12 and e > 0
                at_lo = self.a <= 1e-12 and e < 0
                if not (at_hi or at_lo):
                    self.integral += e * p['ctrl_period_s']
        u = p['kp'] * e + p['ki'] * self.integral
        lo, hi = LIMITS[self.mode]
        self.a_cmd = min(max(p['a0'] + u, lo), hi)

    def step_gate(self, dt):
        lo, hi = LIMITS[self.mode]
        step = self.p['rate_limit_per_s'] * dt
        delta = self.a_cmd - self.a
        if abs(delta) < step:
            self.a = self.a_cmd
        else:
            self.a += step if delta > 0 else -step
        self.a = min(max(self.a, lo), hi)


def run(p, ncell, mode):
    B0, L, ht = p['pool_width_m'], p['pool_len_m'], p['h_target']
    sec = V2.Section(B0, 0.0)
    dx = L / ncell
    S0 = V1.normal_slope(p)
    h = [ht] * ncell; Q = [p['q0']] * ncell
    nsens = max(1, min(ncell, int(round(WIN_M / dx))))
    wmax = p['q0'] / (B0 * ht) + math.sqrt(G * ht)
    dt = V1.pick_dt(dx, wmax)
    nt = int(round(p['duration_s'] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = CtrlX(p, mode)
    c1 = dt / dx; sqrt = math.sqrt; hmin = 1e-3
    ts, dev, gate = [], [], []
    for k in range(nt):
        tk = k * dt
        q_up = p['q0'] + (STEP_DQ if tk >= p['t_step_s'] else 0.0)
        q_gate = p['kg'] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)
        ctrl.control(tk, sum(h[ncell - nsens:]) / nsens)
        ctrl.step_gate(dt)
        hx = [h[0]] + h + [h[-1]]; Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1; F1 = [0.0] * nif; F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = hx[j], hx[j + 1]
            if h1 < hmin: h1 = hmin
            if h2 < hmin: h2 = hmin
            A1, A2 = sec.A(h1), sec.A(h2)
            if A1 < 1e-6: A1 = 1e-6
            if A2 < 1e-6: A2 = 1e-6
            Bc1, Bc2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / Bc1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / Bc2)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            P1 = 0.5 * G * A1 * A1 / Bc1; P2 = 0.5 * G * A2 * A2 / Bc2
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2 + P1 + P2) - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate
        h_new = [0.0] * ncell; Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]; A = sec.A(hi)
            if A < 1e-6: A = 1e-6
            Qi = Q[i]
            R = A / (B0 + 2.0 * hi)
            Sf = MANNING_N * MANNING_N * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        if k % rec_every == 0:
            ts.append(tk); dev.append(sum(h) / ncell - ht); gate.append(ctrl.a)
    return ts, dev, gate


def stats(xs):
    n = len(xs)
    if n == 0: return {'n': 0}
    m = sum(xs) / n
    if n == 1: return {'n': 1, 'mean': m, 'sd': 0.0, 'ci95': [m, m]}
    var = sum((x - m) ** 2 for x in xs) / (n - 1)
    sd = math.sqrt(var)
    h = tq(n - 1) * sd / math.sqrt(n)
    return {'n': n, 'mean': round(m, 4), 'sd': round(sd, 4),
            'ci95': [round(m - h, 4), round(m + h, 4)], 'ci_halfwidth': round(h, 4)}


def find_cycles(ts, dev, t_from, t_to):
    """在 [t_from,t_to] 内找 |dev| 的局部极大（周期峰），返回 (t_peak, amp) 列表。"""
    idx = [i for i, t in enumerate(ts) if t_from <= t <= t_to]
    peaks = []
    for j in range(1, len(idx) - 1):
        i = idx[j]
        a = abs(dev[i])
        if a >= abs(dev[idx[j - 1]]) and a > abs(dev[idx[j + 1]]) and a > 0.05:  # 5 cm 门槛
            peaks.append((ts[i], a))
    # 合并 180s 内相邻峰（同一周期重复检出）
    merged = []
    for t, a in peaks:
        if merged and t - merged[-1][0] < 180.0:
            if a > merged[-1][1]: merged[-1] = (t, a)
        else:
            merged.append((t, a))
    return merged


def window_metrics(ts, dev, t_from, t_to, label):
    pk = find_cycles(ts, dev, t_from, t_to)
    amps = [a * 1000.0 for _, a in pk]
    n_cycle = len(amps)
    sub = [(t, d) for t, d in zip(ts, dev) if t_from <= t <= t_to]
    peak = max(abs(d) for _, d in sub) * 1000.0
    tail_pts = [d for t, d in sub if t >= t_to - 1800.0]
    tail = max(abs(d) for d in tail_pts) * 1000.0 if tail_pts else 0.0
    # 末段斜率
    xs = [t - t_to for t, d in sub if t >= t_to - 1800.0]
    ys = [abs(d) * 1000.0 for t, d in sub if t >= t_to - 1800.0]
    slope = 0.0
    if len(xs) > 2:
        n = len(xs); sx = sum(xs); sy = sum(ys); sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
        den = n * sxx - sx * sx
        if abs(den) > 1e-9: slope = (n * sxy - sx * sy) / den * 60.0  # mm/min
    res = {'window': label, 'peak_mm': round(peak, 1), 'tail_mm': round(tail, 1),
           'tail_slope_mm_per_min': round(slope, 3), 'n_cycle': n_cycle}
    if n_cycle >= 2:
        p2p = [max(amps[i], amps[i - 1]) - min(amps[i], amps[i - 1]) for i in range(1, n_cycle)]
        mid = [0.5 * (amps[i] + amps[i - 1]) for i in range(1, n_cycle)]
        tmid = [0.5 * (pk[i][0] + pk[i - 1][0]) for i in range(1, n_cycle)]
        mu = stats(p2p); kp = stats(mid)
        # 中线趋势
        n = len(tmid); slope_ci = [0.0, 0.0]; width = 0.0
        if n >= 3:
            mx = sum(tmid) / n; my = sum(mid) / n
            sxx = sum((x - mx) ** 2 for x in tmid); sxy = sum((x - mx) * (y - my) for x, y in zip(tmid, mid))
            b = sxy / sxx if sxx > 1e-9 else 0.0
            resid = [my + b * (x - mx) - y for x, y in zip(tmid, mid)]
            s2 = sum(r * r for r in resid) / max(n - 2, 1)
            se = math.sqrt(s2 / sxx) if sxx > 1e-9 else 0.0
            hw = tq(n - 2) * se
            slope_ci = [round((b - hw) * 60.0, 4), round((b + hw) * 60.0, 4)]
            width = round(2 * hw * 60.0, 4)
            res['midline_slope_mm_per_min'] = round(b * 60.0, 4)
        res['midline_slope_ci95'] = slope_ci
        res['midline_slope_ci_width'] = width
        res['p2p_amps_mm'] = [round(x, 1) for x in p2p]
        res['mu_p2p_mm'] = mu
        res['k_p2p_mm'] = kp
        res['cv_p2p'] = round(mu['sd'] / mu['mean'], 3) if mu['mean'] > 1e-9 else None
        res['amp_min_mm'] = round(min(amps), 1)
        res['amp_max_mm'] = round(max(amps), 1)
    return res


def judge_last12h(m):
    """五关判决（末 12h）。"""
    reasons = []
    # ② 零趋势
    ci = m.get('midline_slope_ci95', [0, 0]); w = m.get('midline_slope_ci_width', 99)
    crosses = ci[0] <= 0.0 <= ci[1]
    trend_ok = crosses and w <= 0.10
    if not trend_ok:
        reasons.append(f'②零趋势不达标(CI={ci},宽度={w})')
    # ③ 周期数
    nc = m.get('n_cycle', 0)
    if nc < 3: reasons.append(f'③周期数不足(N={nc})')
    # ④
    cv = m.get('cv_p2p', None)
    if cv is not None and cv > 0.4: reasons.append(f'④CV_p2p>{cv}>0.4')
    # ⑤
    kp = m.get('k_p2p_mm', {})
    kp_lo = kp.get('ci95', [0, 0])[0] if kp else 0
    strong = trend_ok and nc >= 3 and (cv is not None and cv <= 0.4) and kp_lo > 0
    if not reasons:
        return ('强通过' if strong else '弱通过'), reasons
    if len(reasons) == 1 and ('④' in reasons[0] or '②' in reasons[0]):
        return '弱通过', reasons
    return '真失败', reasons


def main():
    t0 = time.time()
    out = {'task': 'A: anti-windup 四臂 12h→24h 延伸 + 五关判决',
           'setup': {'ncell': 480, 'win_m': WIN_M, 'step_dq': STEP_DQ,
                     'duration_h': 24.0, 'm_slope': 0.0,
                     'e4_12h_peaks_mm': E4_12H},
           'arms': {}}
    P = {'pool_len_m': 6000.0, 'pool_width_m': 50.0, 'q0': 180.0, 'h_target': 7.0, 'a0': 0.5,
         'kg': 360.0, 'kp': 1.6, 'ki': 0.0016, 'ctrl_period_s': 300.0, 'deadband_m': 0.01,
         'meas_resolution_m': 0.001, 'rate_limit_per_s': 1.0 / 3000.0, 't_step_s': 900.0,
         'band_normal_m': 0.15, 'band_alarm_m': 0.3, 'settle_tol_m': 0.02,
         'settle_window_s': 1800.0, 'duration_s': DUR_24H}

    for mode in ('base', 'tightclip', 'backcalc', 'noclip'):
        print(f"[A] 跑 {mode}（24h, 480格, 窗长{WIN_M:.0f}m）…", flush=True)
        t1 = time.time()
        ts, dev, gate = run(P, 480, mode)
        wall = time.time() - t1
        m0_12 = window_metrics(ts, dev, 0.0, DUR_12H, '0-12h')
        m12_24 = window_metrics(ts, dev, DUR_12H, DUR_24H, '12-24h')
        mfull = window_metrics(ts, dev, 0.0, DUR_24H, '0-24h')
        overlap = round(min(m0_12['peak_mm'], E4_12H[mode]) / max(m0_12['peak_mm'], E4_12H[mode]) * 100.0, 2)
        v, reasons = judge_last12h(m12_24)
        sat = sum(1 for a in gate if a <= 1e-9 or a >= 1 - 1e-9) / len(gate)
        out['arms'][mode] = {
            'precheck_overlap_pct_0_12h_vs_E4': overlap,
            'precheck_pass': overlap >= 50.0,
            'win_0_12h': m0_12, 'win_12_24h': m12_24, 'win_0_24h': mfull,
            'sat_frac_24h': round(sat, 4),
            'a_min': round(min(gate), 3), 'a_max': round(max(gate), 3),
            'verdict_24h': v, 'verdict_reasons': reasons,
            'wall_s': round(wall, 1),
        }
        print(f"    {mode}: 0-12h peak={m0_12['peak_mm']} (叠合{overlap}%) | "
              f"12-24h peak={m12_24['peak_mm']} tail={m12_24['tail_mm']} N={m12_24['n_cycle']} "
              f"CV={m12_24.get('cv_p2p')} → {v} {reasons}", flush=True)

    # 全局判决：四臂方向性是否在 24h 下保持
    b = out['arms']['base']['win_12_24h']
    tc = out['arms']['tightclip']['win_12_24h']
    bk = out['arms']['backcalc']['win_12_24h']
    nc = out['arms']['noclip']['win_12_24h']
    out['global_verdict'] = {
        'base_still_oscillating': b['n_cycle'] >= 3 and b['peak_mm'] > 100,
        'tightclip_bounded': abs(tc['peak_mm']) < abs(b['peak_mm']),
        'backcalc_bounded': abs(bk['peak_mm']) < abs(b['peak_mm']),
        'noclip_worse': abs(nc['peak_mm']) > abs(b['peak_mm']),
        'peaks_12_24h_mm': {'base': b['peak_mm'], 'tightclip': tc['peak_mm'],
                            'backcalc': bk['peak_mm'], 'noclip': nc['peak_mm']},
        'note': '若四项仍成立 → E4「限幅=因(必要条件)」在 24h 尺度稳健 → §3.x.4 可写「12h/24h 一致」',
    }
    out['wall_total_s'] = round(time.time() - t0, 1)

    os.makedirs('raw', exist_ok=True)
    jp = os.path.join('raw', 'bm01_a_24h_antiwindup.json')
    json.dump(out, open(jp, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
    print(f"\n已写入 {jp}（总耗时 {out['wall_total_s']}s）")
    gv = out['global_verdict']
    print(f"全局：base 仍振荡={gv['base_still_oscillating']} | tightclip有界={gv['tightclip_bounded']} | "
          f"backcalc有界={gv['backcalc_bounded']} | noclip更差={gv['noclip_worse']}")


if __name__ == "__main__":
    main()
