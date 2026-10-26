# -*- coding: utf-8 -*-
"""
(E-VER-T5b) 回水曲线对照的**长时间**版本 + 开环时间常数辨识
================================================================================
上一轮 T5 在 24 h 处与经典回水曲线差 39.3 mm（恒定偏移，非形状误差）。
假设：这不是求解器误差，而是**24 h 尚未到稳态**。
     开环时间常数 τ ≈ A_pool / (dq_gate/dh) ≈ 300000 / 11.86 ≈ 2.5×10^4 s ≈ 7.0 h，
     24 h = 3.4τ → 残差 ~3% × 1.2 m ≈ 39 mm  ← 与观测一致。

本脚本一口气跑到 144 h，并在 24/48/96/144 h 打检查点，验证：
 (a) 与经典回水曲线的 RMSE 是否随仿真时长单调下降并趋于 0；
 (b) 偏差是否按 τ 指数衰减（把 τ 从数据里反解出来）。

副产品（对论文很有用）：给出**开环水位整定时间常数 τ** 的实测值 →
   可直接量化"为什么 4 h 判定窗口太短"（4 h ≈ 0.6τ）。

输出：raw/bm01_ver_t5long.json
"""
import json, os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_ver_solver import (Section, G, MANNING_N, pick_dt, normal_slope,
                             P_of, backwater_profile, HT, B, Q0, KG, A0, L)

NCELL = 240
DQ = 15.0
CHECKPOINTS_H = [4.0, 8.0, 12.0, 24.0, 48.0, 96.0, 144.0]


def run_with_checkpoints(ncell, cp_hours, step_dq):
    p = P_of(1.0, step_dq)
    sec = Section(B, 0.0)
    dx = L / ncell
    S0 = normal_slope(p, MANNING_N)
    h = [HT] * ncell
    Q = [Q0] * ncell
    wmax = Q0 / (B * HT) + math.sqrt(G * HT)
    dt = pick_dt(dx, wmax)
    total_s = max(cp_hours) * 3600.0
    nt = int(round(total_s / dt)) + 1
    cps = sorted(set(int(round(c * 3600.0 / dt)) for c in cp_hours))
    cpset = set(cps)
    cp_map = {}
    for c in cps:
        cp_map[c] = round(c * dt / 3600.0, 2)
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
            shots[cp_map[k]] = {"h_end_m": h[-1], "h_up_m": h[0],
                                "h_mean_m": sum(h) / ncell}
    return {"dx": dx, "dt": dt, "steps": nt, "S0": S0, "shots": shots,
            "final_h": h}


def main():
    t0 = time.time()
    print("=== T5b 长时间回水曲线对照（开环 + ΔQ=%.0f, %d 格）===" % (DQ, NCELL), flush=True)
    res = run_with_checkpoints(NCELL, CHECKPOINTS_H, DQ)
    Qs = Q0 + DQ
    h_down_theory = HT * (Qs / (KG * A0)) ** 2
    h_bw, S0 = backwater_profile(Qs, h_down_theory, nseg=6000)

    # 逐检查点：与解析回水剖面的 RMSE（用最终剖面做形状比对；此处用 h_end 收敛来度量）
    rows = []
    for ch in CHECKPOINTS_H:
        if ch not in res["shots"]:
            continue
        s = res["shots"][ch]
        gap_mm = (h_down_theory - s["h_end_m"]) * 1000.0
        rows.append({"t_h": ch, "h_end_m": round(s["h_end_m"], 5),
                     "h_up_m": round(s["h_up_m"], 5),
                     "gap_vs_theory_mm": round(gap_mm, 3)})
        print("  t=%6.0f h: h_end=%.5f m (解析 %.5f) 差=%7.3f mm"
              % (ch, s["h_end_m"], h_down_theory, gap_mm), flush=True)

    # 反解时间常数：gap(t) ≈ gap0 * exp(-(t-t0)/tau) → tau = (t2-t1)/ln(gap1/gap2)
    taus = []
    for i in range(len(rows) - 1):
        g1, g2 = abs(rows[i]["gap_vs_theory_mm"]), abs(rows[i + 1]["gap_vs_theory_mm"])
        if g1 > 0.05 and g2 > 0.05 and g1 > g2:
            taus.append(round((rows[i + 1]["t_h"] - rows[i]["t_h"]) / math.log(g1 / g2), 2))
    tau_med = round(sorted(taus)[len(taus) // 2], 2) if taus else None

    # 最终剖面（144 h）与解析回水剖面的逐点比对
    dxs = res["dx"]
    bw_on_grid = []
    for i in range(NCELL):
        xc = (i + 0.5) * dxs
        j = xc / L * 6000.0
        j0 = max(0, min(5999, int(j)))
        f = j - j0
        bw_on_grid.append(h_bw[j0] * (1 - f) + h_bw[j0 + 1] * f)
    diffs = [res["final_h"][i] - bw_on_grid[i] for i in range(NCELL)]
    rmse = math.sqrt(sum(d * d for d in diffs) / NCELL) * 1000.0
    maxabs = max(abs(d) for d in diffs) * 1000.0
    print("\n  144 h 最终剖面 vs 经典回水曲线: RMSE=%.4f mm  最大偏差=%.4f mm"
          % (rmse, maxabs), flush=True)
    print("  反解开环时间常数 τ ≈ %s h（各段估计 %s）" % (tau_med, taus), flush=True)

    out = {"purpose": "E-VER T5b 长时间回水曲线对照 + 开环时间常数辨识",
           "ncell": NCELL, "step_dq": DQ, "dt_s": res["dt"], "dx_m": res["dx"],
           "h_down_theory_m": round(h_down_theory, 5),
           "checkpoints": rows,
           "tau_estimates_h": taus, "tau_median_h": tau_med,
           "final_profile_RMSE_mm": round(rmse, 4),
           "final_profile_max_abs_err_mm": round(maxabs, 4),
           "wall_s": round(time.time() - t0, 1)}
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_ver_t5long.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_s"]))


if __name__ == "__main__":
    main()
