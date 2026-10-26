# -*- coding: utf-8 -*-
"""
E1：闸门饱和是"因"还是"果"？——时间分辨的饱和起始点判别
================================================================
背景（P0-1 的遗留问题）：
  P0-1 发现"所有 A3 越带场景闸门都饱和、所有 A1 场景都不饱和"，
  但没判定因果方向：是"局部采样→闸门饱和→失去线性调节→极限环"（饱和=因），
  还是"闭环本身线性不稳定→振荡变大→把闸门顶到饱和"（饱和=果）。

判别逻辑（时间序）：
  对每个失稳场景，从 30 s 采样的时间序列里取：
    t_sat  = 闸门首次触界 (a<=eps 或 a>=1-eps) 的时刻
    t_band = 池均偏差首次越过正常条带 (|dev|>150mm) 的时刻
    t_grow = 池均偏差进入"持续增长"的时刻（用 3 个连续 5min 均值单调递增判）
  若 t_grow 明显早于 t_sat（或 t_sat 缺失但仍在增长）→ 饱和是"果"（先失稳后饱和）
  若 t_sat 早于 t_grow → 饱和是"因"（先饱和后失稳）

另外给出"首次触界时的偏差水平"：若首次触界时 |dev| 仍很小（<0.3m），
说明闸门在偏差还小的时候就被顶到边界，控制器窗口窄 → 支持"因"。
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

P = {
    'pool_len_m': 6000.0, 'pool_width_m': 50.0, 'q0': 180.0, 'h_target': 7.0, 'a0': 0.5,
    'kg': 360.0, 'kp': 1.6, 'ki': 0.0016, 'ctrl_period_s': 300.0, 'deadband_m': 0.01,
    'meas_resolution_m': 0.001, 'rate_limit_per_s': 1.0 / 3000.0, 't_step_s': 900.0,
    'band_normal_m': 0.15, 'band_alarm_m': 0.3, 'settle_tol_m': 0.02,
    'settle_window_s': 1800.0, 'duration_s': 43200.0,
}
EPS = 1e-9

def run(m, win, dur=43200.0):
    kw = dict(ctrl_on=('gate' if win == 0 else 'gate_avg'), m_slope=m)
    if win > 0:
        kw['sensor_len_m'] = float(win)
    if m > 0:
        kw['s0'] = V2.section_normal_slope(V2.Section(P['pool_width_m'], m), P['h_target'], P['q0'], V2.MANNING_N)
    p = dict(P); p['duration_s'] = dur
    t0 = time.time()
    r = V2.sve_1d_v2(p, 480, 45.0, **kw)
    r['_wall_s'] = round(time.time() - t0, 1)
    return r

def analyse(t, dev, gate, thr=0.15, floor=0.005, t_start=900.0):
    n = len(t)
    t_sat = None; dev_at_sat = None
    for i in range(n):
        if gate[i] <= EPS or gate[i] >= 1 - EPS:
            t_sat = t[i]; dev_at_sat = abs(dev[i]); break
    t_band = None
    for i in range(n):
        if abs(dev[i]) > thr:
            t_band = t[i]; break
    # t_grow：仅在阶跃发生之后、且已有可辨信号(>floor)时才判"持续增长"，避免噪声误触发
    t_grow = None
    if n >= 20:
        w = 10
        means = [sum(abs(x) for x in dev[i:i+w]) / w for i in range(0, n - w + 1)]
        for i in range(len(means) - 3):
            if t[i] < t_start or means[i] < floor:
                continue
            if means[i+1] > means[i] and means[i+2] > means[i+1] and means[i+3] > means[i+2] \
               and means[i+3] > 1.5 * means[i]:
                t_grow = t[i]; break
    peak = max(abs(x) for x in dev)
    tail = max(abs(x) for x in dev[-int(1800/30):])
    a_min = min(gate); a_max = max(gate)
    sat_frac = sum(1 for a in gate if a <= EPS or a >= 1 - EPS) / len(gate)
    return {
        't_sat_min': None if t_sat is None else round(t_sat/60, 1),
        'dev_at_first_sat_mm': None if dev_at_sat is None else round(dev_at_sat*1000, 1),
        't_band_min': None if t_band is None else round(t_band/60, 1),
        't_grow_min': None if t_grow is None else round(t_grow/60, 1),
        'peak_mm': round(peak*1000, 1), 'tail_mm': round(tail*1000, 1),
        'a_min': round(a_min, 3), 'a_max': round(a_max, 3), 'sat_frac': round(sat_frac, 4),
    }

def verdict(x):
    tg, tsat = x['t_grow_min'], x['t_sat_min']
    if x['peak_mm'] < 200:      # 未失稳（A1 档）
        return '稳定（未失稳）'
    if tsat is None:
        return '失稳但未饱和 → 饱和非必要（线性机制）'
    if tg is not None and tg < tsat - 1.0:
        return f'先增长后饱和（Δ={tsat-tg:.0f}min）→ 饱和是"果"，线性不稳定'
    if tg is not None and tsat < tg - 1.0:
        return f'先饱和后增长（Δ={tg-tsat:.0f}min）→ 饱和是"因"'
    return '同步（Δ<1min）→ 因/果耦合，需降步长试验'

CASES = [(0.0, 0), (0.0, 100), (0.0, 500), (0.0, 1500), (0.0, 3000), (0.0, 6000),
         (1.5, 0), (1.5, 100), (1.5, 500), (1.5, 1500), (1.5, 3000), (1.5, 6000)]

rows = []
print('E1 闸门饱和：因 or 果？（480格 / 12h）')
print(f"{'断面':>5} {'窗长':>6} | {'t_grow(min)':>11} {'t_sat(min)':>10} {'首触界偏差':>10} | {'峰值':>8} {'末段':>7} {'饱和%':>6} | 判决")
print('-' * 132)
for m, win in CASES:
    r = run(m, win)
    s = r['series']
    x = analyse(s['t'], s['dev_mean'], s['gate'])
    x['m'] = m; x['win'] = win; x['wall_s'] = r['_wall_s']
    x['verdict'] = verdict(x)
    rows.append(x)
    f = lambda v, w, d=1: ('—'.rjust(w) if v is None else f'{v:.{d}f}'.rjust(w))
    print(f"{('m='+str(m)):>5} {win:>6} | {f(x['t_grow_min'],11)} {f(x['t_sat_min'],10)} "
          f"{f(x['dev_at_first_sat_mm'],10)} | {f(x['peak_mm'],8)} {f(x['tail_mm'],7)} "
          f"{x['sat_frac']*100:5.1f}% | {x['verdict']}")

out = {'experiment': 'E1 gate saturation cause-vs-effect', 'ncell': 480, 'duration_h': 12,
       'band_mm': 150, 'rows': rows,
       'note': 't_grow=池均偏差持续增长起点；t_sat=闸门首次触界；dev_at_first_sat=首次触界时的偏差水平'}
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw'), exist_ok=True)
jp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw', 'bm01_e1_saturation_causality.json')
json.dump(out, open(jp, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(f'\n输出：{jp}')
