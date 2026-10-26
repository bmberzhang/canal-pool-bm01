# -*- coding: utf-8 -*-
"""
E2：空间模态分解 —— 失稳到底是"哪个空间模态"在长？
================================================================
动机：
  P1-1' 显示峰值随采样窗长**非单调**（中间窗长 500~1500m 最不稳），这推翻了
  "空间平均=低通=单调稳定化"的线性频域预期。要把它从"现象"升级为"机制"，
  必须回答：失稳时到底是哪个**空间波数**在增长？它的波长是否与采样窗长挂钩？

方法（不改物理，只加观测）：
  1) 用与 bm01_sve1d_v2 逐行一致的一维求解器（自检保证），每 300 s 抓一次全场 h(x)；
  2) 对偏差场 e(x,t)=h(x,t)-mean_x h(x,t) 做空间 DFT，取 j=1..JMAX（λ=L/j）；
  3) 追踪各模态振幅 A_j(t)：
     - 晚期主导模态 j*（λ*=L/j*）；
     - 该模态的指数增长率 σ（对 lnA 线性拟合）；
  4) 检验两个假设：
     H-phase：λ* 随采样窗长 L_s 单调/成比例变化（相位延迟机制 → 采样窗长调制回路相位）；
     H-sat  ：λ* 与 L_s 无关（饱和极限环 → 模态只是饱和后自然涌现的振荡，与窗长无关）。

判据：若 λ*(L_s) 与 L_s 显著同向变化（相关系数高）→ 支持 H-phase；
      若 λ* 基本恒定 → 支持 H-sat。
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


def sve_field(p, ncell, step_dq, *, ctrl_on='gate', m_slope=0.0, n=MANNING_N, cfl=CFL,
              s0=None, sensor_len_m=100.0, snap_every_s=300.0):
    """与 bm01_sve1d_v2.sve_1d_v2 逐行一致的求解器，额外记录全场快照。"""
    B0, L, ht = p['pool_width_m'], p['pool_len_m'], p['h_target']
    sec = V2.Section(B0, m_slope)
    dx = L / ncell
    S0 = V1.normal_slope(p, n) if s0 is None else s0
    h = [ht] * ncell
    Q = [p['q0']] * ncell
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))
    wmax = p['q0'] / (B0 * ht) + math.sqrt(G * ht)
    lo = cfl * dx / wmax
    dt = V1.pick_dt(dx, wmax) if abs(cfl - CFL) < 1e-15 else None
    if dt is None:
        dt = next(c for c in V1.DT_CANDIDATES if c <= lo)
    nt = int(round(p['duration_s'] / dt)) + 1
    snap_every = max(1, int(round(snap_every_s / dt)))
    ctrl = V1.Ctrl(p)
    c1 = dt / dx
    sqrt = math.sqrt
    hmin = 1e-3
    friction_correct = not sec.is_rect()

    snap_t, snap_h, snap_gate, snap_dev = [], [], [], []
    for k in range(nt):
        tk = k * dt
        q_up = p['q0'] + (step_dq if tk >= p['t_step_s'] else 0.0)
        q_gate = p['kg'] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)
        if ctrl_on != 'off':
            if ctrl_on == 'gate':
                meas = h[-1]
            elif ctrl_on == 'gate_avg':
                meas = sum(h[ncell - nsens:]) / nsens
            else:
                meas = sum(h) / ncell
            ctrl.control(tk, meas)
            ctrl.step_gate(dt)

        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif; F2 = [0.0] * nif
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
            P1 = 0.5 * G * A1 * A1 / Bc1
            P2 = 0.5 * G * A2 * A2 / Bc2
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2 + P1 + P2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate

        h_new = [0.0] * ncell; Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]; A = sec.A(hi)
            if A < 1e-6: A = 1e-6
            Qi = Q[i]
            R = (A / sec.P(hi) if sec.P(hi) > 0 else hi) if friction_correct else A / (B0 + 2.0 * hi)
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        if k % snap_every == 0:
            snap_t.append(tk); snap_h.append(list(h))
            snap_gate.append(ctrl.a); snap_dev.append(sum(h) / ncell - ht)
    return {'ncell': ncell, 'dx_m': dx, 'dt_s': dt, 'steps': nt,
            'snap_t': snap_t, 'snap_h': snap_h, 'snap_gate': snap_gate, 'snap_dev': snap_dev}


def spectrum(e, L, jmax):
    """对偏差场 e(x) 做空间 DFT，返回 j=1..jmax 的振幅 A_j（物理量：m）。"""
    n = len(e)
    dx = L / n
    amps = [0.0] * (jmax + 1)
    for j in range(1, jmax + 1):
        kk = 2.0 * math.pi * j / L
        re = im = 0.0
        for i in range(n):
            ang = kk * i * dx
            re += e[i] * math.cos(ang)
            im += e[i] * math.sin(ang)
        amps[j] = 2.0 * math.sqrt(re * re + im * im) / n
    return amps


def analyse(r, jmax=40):
    L = P['pool_len_m']; n = r['ncell']
    times = r['snap_t']; devs = r['snap_dev']
    A = []  # A[k][j]
    for hf in r['snap_h']:
        m = sum(hf) / n
        e = [v - m for v in hf]
        A.append(spectrum(e, L, jmax))
    # 晚期窗口：末 30% 快照；早期窗口：第 20%~40%
    K = len(A)
    late = range(int(K * 0.75), K)
    early = range(int(K * 0.2), int(K * 0.4))
    def dom(idx):
        s = [0.0] * (jmax + 1)
        for k in idx:
            for j in range(1, jmax + 1):
                s[j] += A[k][j] * A[k][j]
        jstar = max(range(1, jmax + 1), key=lambda j: s[j])
        return jstar, math.sqrt(s[jstar] / max(1, len(list(idx))))
    j_late, a_late = dom(late)
    j_early, a_early = dom(early)
    # 主导模态增长率 σ：对 ln(A[k][j_late]) 在晚期窗口线性拟合
    ts = [times[k] for k in late if A[k][j_late] > 1e-9]
    ys = [math.log(A[k][j_late]) for k in late if A[k][j_late] > 1e-9]
    sigma = None
    if len(ts) >= 3:
        npts = len(ts); sx = sum(ts); sy = sum(ys)
        sxx = sum(x * x for x in ts); sxy = sum(x * y for x, y in zip(ts, ys))
        den = npts * sxx - sx * sx
        if abs(den) > 1e-9:
            sigma = (npts * sxy - sx * sy) / den
    peak = max(abs(d) for d in devs) * 1000
    tail = max(abs(d) for d in devs[-int(1800 / 300):]) * 1000
    return {'j_star': j_late, 'lambda_star_m': round(L / j_late, 1),
            'amp_late_mm': round(a_late * 1000, 1), 'amp_early_mm': round(a_early * 1000, 1),
            'sigma_per_s': (None if sigma is None else round(sigma, 6)),
            'sigma_per_h': (None if sigma is None else round(sigma * 3600, 4)),
            'peak_mm': round(peak, 1), 'tail_mm': round(tail, 1)}


CASES = [(0.0, 0), (0.0, 100), (0.0, 500), (0.0, 1500), (0.0, 3000), (0.0, 6000),
         (1.5, 0), (1.5, 100), (1.5, 500), (1.5, 1500), (1.5, 3000), (1.5, 6000)]
rows = []
print('E2 空间模态分解（480格 / 12h / 快照每300s / j=1..40）')
print(f"{'断面':>5} {'窗长':>6} | {'λ*(m)':>8} {'j*':>4} {'晚期振幅mm':>10} {'早期振幅mm':>10} {'σ(/h)':>9} | {'峰值mm':>8} {'末段mm':>7}")
print('-' * 110)
for m, win in CASES:
    kw = dict(ctrl_on=('gate' if win == 0 else 'gate_avg'), m_slope=m)
    if win > 0:
        kw['sensor_len_m'] = float(win)
    if m > 0:
        kw['s0'] = V2.section_normal_slope(V2.Section(P['pool_width_m'], m), P['h_target'], P['q0'], MANNING_N)
    t0 = time.time()
    r = sve_field(P, 480, 45.0, **kw)
    x = analyse(r)
    x['m'] = m; x['win'] = win; x['wall_s'] = round(time.time() - t0, 1)
    rows.append(x)
    s = '—'.rjust(9) if x['sigma_per_h'] is None else f"{x['sigma_per_h']:.4f}"
    print(f"{('m='+str(m)):>5} {win:>6} | {x['lambda_star_m']:>8} {x['j_star']:>4} "
          f"{x['amp_late_mm']:>10} {x['amp_early_mm']:>10} {s:>9} | {x['peak_mm']:>8} {x['tail_mm']:>7}")

out = {'experiment': 'E2 spatial modal decomposition', 'ncell': 480, 'duration_h': 12,
       'jmax': 40, 'pool_len_m': P['pool_len_m'],
       'hypothesis': 'H-phase: λ*随窗长同向变化 → 相位延迟机制；H-sat: λ*与窗长无关 → 饱和极限环',
       'rows': rows}
jp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw', 'bm01_e2_modal_decomp.json')
json.dump(out, open(jp, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(f'\n输出：{jp}')
