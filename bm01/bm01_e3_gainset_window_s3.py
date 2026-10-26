# -*- coding: utf-8 -*-
"""
E3：S3 生死判 —— "采样窗长非单调"是物理现象，还是基准 PI 整定的病态？
================================================================
决策备忘把这一条列为"唯一能在一个月内改变路线的动作"。
判据：换掉基准增益 (kp=1.6, ki=0.0016) 后，非单调性（中间窗长峰值最高）是否仍然出现？
  - 仍然出现（不止一组增益）→ S3 排除，非单调是稳健现象，一区增量才有意义；
  - 只有基准增益出现      → S3 成立，属于整定病态，立即转二区主线。
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BASE = {
    'pool_len_m': 6000.0, 'pool_width_m': 50.0, 'q0': 180.0, 'h_target': 7.0, 'a0': 0.5,
    'kg': 360.0, 'kp': 1.6, 'ki': 0.0016, 'ctrl_period_s': 300.0, 'deadband_m': 0.01,
    'meas_resolution_m': 0.001, 'rate_limit_per_s': 1.0 / 3000.0, 't_step_s': 900.0,
    'band_normal_m': 0.15, 'band_alarm_m': 0.3, 'settle_tol_m': 0.02,
    'settle_window_s': 1800.0, 'duration_s': 43200.0,
}
GAINS = [('V1_base', 1.6, 0.0016), ('V2_cons', 0.8, 0.0008), ('V3_aggr', 3.2, 0.0032)]
WINS = [0, 100, 500, 1500, 3000, 6000]


def run(m, win, kp, ki):
    kw = dict(ctrl_on=('gate' if win == 0 else 'gate_avg'), m_slope=m)
    if win > 0:
        kw['sensor_len_m'] = float(win)
    if m > 0:
        kw['s0'] = V2.section_normal_slope(V2.Section(BASE['pool_width_m'], m), BASE['h_target'], BASE['q0'], V2.MANNING_N)
    p = dict(BASE); p['kp'] = kp; p['ki'] = ki
    r = V2.sve_1d_v2(p, 480, 45.0, **kw)
    peak = r['peak_mean_mm']; over = r['over_band_s']
    code = 'A1' if peak < 150.0 else 'A3'
    return {'peak_mm': round(peak, 1), 'final_mm': round(r['final_mean_mm'], 1),
            'over_band_min': round(over / 60.0, 1), 'code': code}


def shape(pk):
    """返回窗长-峰值曲线的形状描述与是否非单调。"""
    imax = max(range(len(pk)), key=lambda i: pk[i])
    monotone_dec = all(pk[i] >= pk[i + 1] - 1e-9 for i in range(len(pk) - 1))
    mid_max = 0 < imax < len(pk) - 1           # 峰值出现在中间（两端都不是最大）
    return {'argmax_win': WINS[imax], 'monotone_decreasing': monotone_dec, 'mid_max': mid_max}


rows = []
print('E3 S3 生死判：增益集 × 采样窗长扫描（m=0 矩形 / 480格 / 12h / 阶跃+45）')
print(f"{'增益集':>8} {'kp':>6} {'ki':>8} | " + ' '.join(f'{w:>8}m' for w in WINS) + ' |  形状判定')
print('-' * 108)
for name, kp, ki in GAINS:
    pk = []; codes = []
    for win in WINS:
        r = run(0.0, win, kp, ki)
        pk.append(r['peak_mm']); codes.append(r['code'])
        rows.append({'gain_set': name, 'kp': kp, 'ki': ki, 'win_m': win, **r})
    sh = shape(pk)
    tag = ('非单调(中间窗长最高)' if sh['mid_max'] else
           '单调递减' if sh['monotone_decreasing'] else '其他(见数据)')
    print(f"{name:>8} {kp:>6} {ki:>8} | " + ' '.join(f'{v:>8.1f}' for v in pk) + f' |  {tag} (argmax={sh["argmax_win"]}m)')

# 汇总判定
nonmono = [g[0] for g in GAINS
           if shape([r2['peak_mm'] for r2 in rows if r2['gain_set'] == g[0]])['mid_max']]
verdict = ('S3 排除：多组增益均出现中间窗长最不稳 → 非单调是稳健现象（一区增量有意义）'
           if len(nonmono) >= 2 else
           'S3 成立：仅基准增益出现非单调 → 属于整定病态（立即转二区主线）')
print(f'\n出现非单调的增益集：{nonmono}')
print(f'判定：{verdict}')

jp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'raw', 'bm01_e3_gainset_window_s3.json')
json.dump({'experiment': 'E3 S3 gain-set x window scan', 'wins': WINS,
           'gains': [{'name': n, 'kp': k, 'ki': i} for n, k, i in GAINS],
           'rows': rows, 'nonmonotone_gain_sets': nonmono, 'verdict': verdict,
           'note': 'code: A1=峰值<150mm, A3=峰值>=150mm；mid_max=峰值出现在中间窗长'},
          open(jp, 'w', encoding='utf-8'), ensure_ascii=False, indent=2)
print(f'输出：{jp}')
