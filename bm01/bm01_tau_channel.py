# -*- coding: utf-8 -*-
"""
(W-tau-b) 通道对照：同一个开环阶跃，分别辨识
   通道1 = 闸前（下游端点）水深 h_end
   通道2 = 全池平均水深 h_mean
检验两者是否给出同一个开环整定时间常数 τ。

动机：正文 §3.5 声称 τ 是"上游阶跃→全池平均水位"通道，
     而原辨识脚本 bm01_ver_t5long.py 实际用的是 h_end。必须把两者都算出来对齐口径。

输出：raw/bm01_tau_channel.json
"""
import json, os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_ver_solver import (Section, G, MANNING_N, pick_dt, normal_slope,
                             P_of, backwater_profile, HT, B, Q0, KG, A0, L)

NCELL = 240
DQ_LIST = [15.0, 45.0]
CHECKPOINTS_H = [4.0, 8.0, 12.0, 24.0, 48.0, 96.0, 144.0]
FIT_LO, FIT_HI = 4.0, 48.0
TOL_MM = 3.0


def run_open_loop(ncell, cp_hours, step_dq):
    p = P_of(1.0, step_dq)
    dx = L / ncell
    S0 = normal_slope(p, MANNING_N)
    h = [HT] * ncell
    Q = [Q0] * ncell
    Qs_guess = Q0 + step_dq
    h_ss_guess = HT * (Qs_guess / (KG * A0)) ** 2
    wmax = max(Q0 / (B * HT) + math.sqrt(G * HT),
               Qs_guess / (B * h_ss_guess) + math.sqrt(G * h_ss_guess))
    dt = pick_dt(dx, wmax)
    total_s = max(cp_hours) * 3600.0
    nt = int(round(total_s / dt)) + 1
    cps = sorted(set(int(round(c * 3600.0 / dt)) for c in cp_hours))
    cpset = set(cps)
    cp_map = {c: round(c * dt / 3600.0, 2) for c in cps}
    sqrt = math.sqrt
    hmin = 1e-3
    c1 = dt / dx
    shots = {}
    for k in range(nt):
        tk = k * dt
        q_up = Q0 + (step_dq if tk >= 900.0 else 0.0)
        q_gate = KG * A0 * sqrt(max(h[-1], hmin)) / sqrt(HT)
        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif
        F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = max(hx[j], hmin), max(hx[j + 1], hmin)
            A1, A2 = B * h1, B * h2
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / B)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / B)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2
                           + 0.5 * G * A1 * A1 / B + 0.5 * G * A2 * A2 / B) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate
        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            A = B * hi
            Qi = Q[i]
            R = A / (B + 2.0 * hi)
            Sf = MANNING_N * MANNING_N * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            hn = hi - c1 * (F1[i + 1] - F1[i]) / B
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        if k in cpset:
            shots[cp_map[k]] = {"h_end_m": h[-1], "h_mean_m": sum(h) / ncell}
    return {"dx": dx, "dt": dt, "steps": nt, "shots": shots}


def linreg(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx else 0.0
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    return a, b, (1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"))


def identify(rows, key, h_ss):
    """对某通道做逐段 τ + log 回归。"""
    taus = []
    for i in range(len(rows) - 1):
        t1, t2 = rows[i]["t_h"], rows[i + 1]["t_h"]
        if t1 < FIT_LO - 1e-9 or t2 > FIT_HI + 1e-9:
            continue
        g1 = abs(h_ss - rows[i][key]) * 1000.0
        g2 = abs(h_ss - rows[i + 1][key]) * 1000.0
        if g1 > 0.05 and g2 > 0.05 and g1 > g2:
            taus.append({"seg_h": [t1, t2], "tau_h": round((t2 - t1) / math.log(g1 / g2), 3)})
    tv = [x["tau_h"] for x in taus]
    srt = sorted(tv)
    med = srt[len(srt) // 2] if tv else None
    xs = [r["t_h"] for r in rows if FIT_LO - 1e-9 <= r["t_h"] <= FIT_HI + 1e-9
          and abs(h_ss - r[key]) * 1000.0 > 0.05]
    ys = [math.log(abs(h_ss - r[key]) * 1000.0) for r in rows
          if FIT_LO - 1e-9 <= r["t_h"] <= FIT_HI + 1e-9
          and abs(h_ss - r[key]) * 1000.0 > 0.05]
    a, b, r2 = linreg(xs, ys)
    return {"tau_segments": taus, "tau_median_h": med,
            "tau_logfit_h": round(-1.0 / b, 3) if b else None,
            "logfit_R2": round(r2, 6)}


def main():
    t0 = time.time()
    out = {"purpose": "W-tau-b 通道对照：h_end（闸前） vs h_mean（全池平均）",
           "channel_defs": {"h_end": "闸前/下游端点水深", "h_mean": "全池空间平均水深"},
           "control_state": "开环（控制器关闭，闸门开度固定 a0=0.5 m2）",
           "fit_interval_h": [FIT_LO, FIT_HI], "settle_tolerance_mm": TOL_MM,
           "ncell": NCELL, "dx_m": L / NCELL, "cases": {}}

    for dq in DQ_LIST:
        print("\n=== ΔQ = %.0f m3/s ===" % dq, flush=True)
        res = run_open_loop(NCELL, CHECKPOINTS_H, dq)
        Qs = Q0 + dq
        h_ss = HT * (Qs / (KG * A0)) ** 2
        rows = []
        for ch in CHECKPOINTS_H:
            s = res["shots"][ch]
            rows.append({"t_h": ch, "h_end_m": round(s["h_end_m"], 6),
                         "h_mean_m": round(s["h_mean_m"], 6),
                         "gap_end_mm": round((h_ss - s["h_end_m"]) * 1000.0, 4),
                         # 池均的稳态值不等于 h_ss（有水面坡降），这里先记原始值，
                         # 用 144h 的 h_mean 作为其自稳态参考（下面重算）
                         })
        h_mean_ss = rows[-1]["h_mean_m"]
        for r in rows:
            r["gap_mean_mm"] = round((h_mean_ss - r["h_mean_m"]) * 1000.0, 4)
        e = identify(rows, "h_end_m", h_ss)
        m = identify(rows, "h_mean_m", h_mean_ss)
        print("  h_end : τ=%s h (R²=%.5f, logfit %.3f)" % (e["tau_median_h"], e["logfit_R2"], e["tau_logfit_h"]))
        print("  h_mean: τ=%s h (R²=%.5f, logfit %.3f)" % (m["tau_median_h"], m["logfit_R2"], m["tau_logfit_h"]))
        gap_end4 = next(r["gap_end_mm"] for r in rows if r["t_h"] == 4.0)
        gap_mean4 = next(r["gap_mean_mm"] for r in rows if r["t_h"] == 4.0)
        out["cases"]["dq_%02.0f" % dq] = {
            "dq_m3s": dq, "Qs_m3s": Qs,
            "h_ss_downstream_m": round(h_ss, 6), "h_ss_poolmean_m": round(h_mean_ss, 6),
            "rows": rows, "channel_h_end": e, "channel_h_mean": m,
            "gap_at_4h_end_mm": gap_end4, "gap_at_4h_mean_mm": gap_mean4,
            "channel_tau_diff_pct": round((m["tau_median_h"] - e["tau_median_h"]) / e["tau_median_h"] * 100.0, 2)
            if e["tau_median_h"] else None,
        }
        print("  两通道 τ 差 %.2f%% | 4h残余: 闸前 %.1f mm, 池均 %.1f mm"
              % (out["cases"]["dq_%02.0f" % dq]["channel_tau_diff_pct"], gap_end4, gap_mean4), flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_tau_channel.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_s"]), flush=True)


if __name__ == "__main__":
    main()
