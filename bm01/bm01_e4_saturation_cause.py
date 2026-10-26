# -*- coding: utf-8 -*-
"""
E4：饱和是"因"还是"果"？—— 小信号线性检验 + 反饱和三臂
================================================================
决策备忘要求区分：
  H_因：闸门限幅（积分饱和/windup）导致失稳；H_果：线性不稳定先发生，闸门触限只是伴随。
判别法（两条互补）：
  (1) 小信号（阶跃幅值扫描）：把阶跃从 +45 缩到 +20/+8/+3 m³/s（偏差随之缩小、闸门远离边界），
      看振荡是否仍在增长：
        · 幅值越小仍增长（增长率与幅值无关）→ 线性不稳定 → 饱和是"果"；
        · 幅值小则衰减、只有大幅值才增长     → 非线性（饱和）驱动 → 饱和是"因"。
  (2) 三臂反饱和：
        base     : 现有 PI，a∈[0,1]
        L1 tight : a∈[0.05,0.95]（收紧物理限幅）
        L3 aw    : back-calculation 抗饱和（未饱和命令落在 [0,1] 内才积分）
        noclip   : 诊断臂，a∈[-0.2,1.2]（把限幅非线性拿掉）
      若 noclip 稳定而 base 失稳 → 限幅非线性是失稳的必要条件 → 饱和是"因"。
固定：m=0 矩形 / 480 格 / 12h / 采样窗长取最恶劣档 1500m。
"""
import os, sys, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2
import bm01_sve1d_compare as V1

G = V1.G; MANNING_N = V1.MANNING_N; CFL = V1.CFL
P = {
    'pool_len_m': 6000.0, 'pool_width_m': 50.0, 'q0': 180.0, 'h_target': 7.0, 'a0': 0.5,
    'kg': 360.0, 'kp': 1.6, 'ki': 0.0016, 'ctrl_period_s': 300.0, 'deadband_m': 0.01,
    'meas_resolution_m': 0.001, 'rate_limit_per_s': 1.0 / 3000.0, 't_step_s': 900.0,
    'band_normal_m': 0.15, 'band_alarm_m': 0.3, 'settle_tol_m': 0.02,
    'settle_window_s': 1800.0, 'duration_s': 43200.0,
}
LIMITS = {'base': (0.0, 1.0), 'tightclip': (0.05, 0.95), 'backcalc': (0.0, 1.0), 'noclip': (-0.2, 1.2)}


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
        u_unsat = p['kp'] * e + p['ki'] * self.integral
        a_uns = p['a0'] + u_unsat
        if self.mode == 'backcalc':
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


def sve_ctrl(p, ncell, step_dq, *, mode='base', win=1500.0, m_slope=0.0, s0=None):
    B0, L, ht = p['pool_width_m'], p['pool_len_m'], p['h_target']
    sec = V2.Section(B0, m_slope)
    dx = L / ncell
    S0 = V1.normal_slope(p) if s0 is None else s0
    h = [ht] * ncell; Q = [p['q0']] * ncell
    nsens = max(1, min(ncell, int(round(win / dx))))
    wmax = p['q0'] / (B0 * ht) + math.sqrt(G * ht)
    dt = V1.pick_dt(dx, wmax)
    nt = int(round(p['duration_s'] / dt)) + 1
    rec_every = max(1, int(round(30.0 / dt)))
    ctrl = CtrlX(p, mode)
    c1 = dt / dx; sqrt = math.sqrt; hmin = 1e-3
    friction_correct = not sec.is_rect()
    ts, dev, gate = [], [], []
    for k in range(nt):
        tk = k * dt
        q_up = p['q0'] + (step_dq if tk >= p['t_step_s'] else 0.0)
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
            R = (A / sec.P(hi) if sec.P(hi) > 0 else hi) if friction_correct else A / (B0 + 2.0 * hi)
            Sf = MANNING_N * MANNING_N * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        if k % rec_every == 0:
            ts.append(tk); dev.append(sum(h) / ncell - ht); gate.append(ctrl.a)
    peak = max(abs(d) for d in dev) * 1000
    tail = max(abs(d) for d in dev[-int(1800 / 30):]) * 1000
    sat = sum(1 for a in gate if a <= 1e-9 or a >= 1 - 1e-9) / len(gate)
    # 末段斜率（mm/min）：对最后 1800s 的 |dev| 线性拟合
    n_last = int(1800 / 30); xs = ts[-n_last:]; ys = [abs(d) for d in dev[-n_last:]]
    npt = len(xs); sx = sum(xs); sy = sum(ys); sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    den = npt * sxx - sx * sx
    slope = (npt * sxy - sx * sy) / den * 60 * 1000 if abs(den) > 1e-9 else 0.0
    return {'peak_mm': round(peak, 1), 'tail_mm': round(tail, 1), 'tail_slope_mm_per_min': round(slope, 3),
            'sat_frac': round(sat, 4), 'a_min': round(min(gate), 3), 'a_max': round(max(gate), 3)}


def verdict(x):
    grow = x['tail_slope_mm_per_min'] > 0.01 and x['tail_mm'] > 100
    return '增长/失稳' if grow else '衰减/有界'


rows = []
print('E4 饱和因/果判别（m=0 / 480格 / 12h / 采样窗长1500m）')
print('--- (1) 小信号阶跃幅值扫描（base PI）---')
print(f"{'阶跃':>6} | {'峰值mm':>9} {'末段mm':>8} {'末段斜率mm/min':>14} {'饱和%':>6} | 结论")
for dq in (45.0, 20.0, 8.0, 3.0):
    x = sve_ctrl(P, 480, dq, mode='base')
    x.update({'arm': 'amp', 'step_dq': dq})
    x['verdict'] = verdict(x); rows.append(x)
    print(f"{dq:>6.0f} | {x['peak_mm']:>9} {x['tail_mm']:>8} {x['tail_slope_mm_per_min']:>14} {x['sat_frac']*100:>5.1f}% | {x['verdict']}")

print('\n--- (2) 反饱和三臂（阶跃+45）---')
print(f"{'臂':>10} | {'峰值mm':>9} {'末段mm':>8} {'末段斜率mm/min':>14} {'饱和%':>6} {'a范围':>14} | 结论")
for mode in ('base', 'tightclip', 'backcalc', 'noclip'):
    x = sve_ctrl(P, 480, 45.0, mode=mode)
    x.update({'arm': 'antiwindup', 'mode': mode})
    x['verdict'] = verdict(x); rows.append(x)
    print(f"{mode:>10} | {x['peak_mm']:>9} {x['tail_mm']:>8} {x['tail_slope_mm_per_min']:>14} {x['sat_frac']*100:>5.1f}% "
          f"{str((x['a_min'],x['a_max'])):>14} | {x['verdict']}")

jp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw', 'bm01_e4_saturation_cause.json')
json.dump({'experiment': 'E4 saturation cause-vs-effect', 'ncell': 480, 'win_m': 1500,
           'rows': rows,
           'note': '小信号仍增长→线性不稳定（饱和=果）；小信号衰减仅大幅值增长→饱和=因。'
                   'noclip稳定而base失稳→限幅是必要条件（饱和=因）。'},
          open(jp, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(f'\n输出：{jp}')
