# -*- coding: utf-8 -*-
"""
(W-tau) 开环整定时间常数的系统辨识 + 不确定度 + 幅值依赖
================================================================================
目的：回应审稿意见"τ 的识别方法、拟合区间、拟合优度、不确定度；
      τ 是哪个输入—输出通道、在哪个控制状态下识别的"。

辨识设置（必须写进论文）：
  通道：输入 = 上游流量阶跃 ΔQ；输出 = 闸前（下游端点）水深 h_end
  控制状态：开环（控制器关闭，闸门开度固定 a0=0.5 m²）
  τ 定义：h_end(t) 相对其稳态值的残余 gap(t) 的指数衰减常数
          gap(t) = h_ss - h_end(t),  τ = (t2-t1)/ln(gap1/gap2)
  拟合区间：4–48 h（gap 已进入单指数区；96 h 后 gap→0，比值噪声放大，不参与）
  判稳指标：|gap| < tol，tol 取 3 mm（=150 mm 条带的 2%）

同时给出：
  (a) 逐段 τ 估计（4 段）+ 中位/均值/标准差（离散度即不确定度）
  (b) log(gap)-t 线性回归的 R²（指数律拟合优度）
  (c) 解析交叉验证 τ_lin = A_pool·2·h_ss/Qs（蓄量—泄流线性化）
  (d) 幅值依赖：ΔQ ∈ {10,15,30,45}，验证 τ ∝ (Q0+ΔQ)

输出：raw/bm01_tau_identify.json
"""
import json, os, sys, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_ver_solver import (Section, G, MANNING_N, pick_dt, normal_slope,
                             P_of, backwater_profile, HT, B, Q0, KG, A0, L)

NCELL = 240
DQ_LIST = [10.0, 15.0, 30.0, 45.0]
CHECKPOINTS_H = [4.0, 8.0, 12.0, 24.0, 48.0, 96.0, 144.0]
FIT_LO, FIT_HI = 4.0, 48.0     # 拟合区间
TOL_MM = 3.0                   # 判稳容差（150mm 条带的 2%）


def run_open_loop(ncell, cp_hours, step_dq):
    """开环：闸门开度固定 a0，控制器关闭；上游在 900 s 处阶跃 step_dq。"""
    p = P_of(1.0, step_dq)
    dx = L / ncell
    S0 = normal_slope(p, MANNING_N)
    h = [HT] * ncell
    Q = [Q0] * ncell
    # CFL：按该幅值的稳态波速取 dt（阶跃后 wmax 升高，用较小 dt 保证稳定）
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
            shots[cp_map[k]] = {"h_end_m": h[-1], "h_up_m": h[0]}
    return {"dx": dx, "dt": dt, "steps": nt, "shots": shots}


def linreg(xs, ys):
    """最小二乘 y = a + b x，返回 a, b, R²。"""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx else 0.0
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return a, b, r2


def main():
    t0 = time.time()
    out = {
        "purpose": "W-tau 开环整定时间常数系统辨识（通道/状态/拟合区间/优度/不确定度/幅值依赖）",
        "channel": {"input": "上游流量阶跃 ΔQ [m3/s]",
                    "output": "闸前（下游端点）水深 h_end [m]"},
        "control_state": "开环（控制器关闭，闸门开度固定 a0=0.5 m2）",
        "fit_interval_h": [FIT_LO, FIT_HI],
        "settle_tolerance_mm": TOL_MM,
        "pool": {"B_m": B, "L_m": L, "Q0_m3s": Q0, "HT_m": HT, "KG": KG, "A0_m2": A0,
                 "A_pool_m2": B * L},
        "ncell": NCELL, "dx_m": L / NCELL,
        "cases": {},
    }

    for dq in DQ_LIST:
        print("\n=== ΔQ = %.0f m3/s ===" % dq, flush=True)
        res = run_open_loop(NCELL, CHECKPOINTS_H, dq)
        Qs = Q0 + dq
        # 稳态解析解：q_gate = KG*a0*sqrt(h)/sqrt(HT) = Qs  =>  h_ss = HT*(Qs/(KG*a0))^2
        h_ss = HT * (Qs / (KG * A0)) ** 2
        rows = []
        for ch in CHECKPOINTS_H:
            s = res["shots"][ch]
            gap_mm = (h_ss - s["h_end_m"]) * 1000.0
            rows.append({"t_h": ch, "h_end_m": round(s["h_end_m"], 6),
                         "gap_mm": round(gap_mm, 4)})
            print("  t=%5.0f h  h_end=%.6f  gap=%9.4f mm" % (ch, s["h_end_m"], gap_mm),
                  flush=True)

        # (a) 逐段 τ（只用拟合区间内、且 gap 足够大的相邻对）
        taus = []
        for i in range(len(rows) - 1):
            t1, t2 = rows[i]["t_h"], rows[i + 1]["t_h"]
            if t1 < FIT_LO - 1e-9 or t2 > FIT_HI + 1e-9:
                continue
            g1, g2 = abs(rows[i]["gap_mm"]), abs(rows[i + 1]["gap_mm"])
            if g1 > 0.05 and g2 > 0.05 and g1 > g2:
                taus.append({"seg_h": [t1, t2],
                             "tau_h": round((t2 - t1) / math.log(g1 / g2), 3)})
        tv = [x["tau_h"] for x in taus]
        tv_sorted = sorted(tv)
        med = tv_sorted[len(tv_sorted) // 2] if tv else None
        mean = sum(tv) / len(tv) if tv else None
        std = math.sqrt(sum((x - mean) ** 2 for x in tv) / (len(tv) - 1)) if len(tv) > 1 else 0.0

        # (b) log(gap)-t 线性回归（拟合区间内）
        xs = [r["t_h"] for r in rows if FIT_LO - 1e-9 <= r["t_h"] <= FIT_HI + 1e-9
              and abs(r["gap_mm"]) > 0.05]
        ys = [math.log(abs(r["gap_mm"])) for r in rows if FIT_LO - 1e-9 <= r["t_h"] <= FIT_HI + 1e-9
              and abs(r["gap_mm"]) > 0.05]
        a, b, r2 = linreg(xs, ys)
        tau_fit = -1.0 / b if b else float("nan")

        # (c) 解析交叉验证
        dqdh = Qs / (2.0 * h_ss)
        tau_lin_h = (B * L) / dqdh / 3600.0

        # (d) 4h 处残余与"4h/τ"
        gap4 = next(r["gap_mm"] for r in rows if abs(r["t_h"] - 4.0) < 1e-9)
        out["cases"]["dq_%02.0f" % dq] = {
            "dq_m3s": dq, "Qs_m3s": Qs, "h_ss_m": round(h_ss, 6),
            "checkpoints": rows,
            "tau_segments": taus,
            "tau_median_h": med, "tau_mean_h": round(mean, 3) if mean else None,
            "tau_std_h": round(std, 3),
            "tau_range_h": [min(tv), max(tv)] if tv else None,
            "tau_logfit_h": round(tau_fit, 3), "logfit_R2": round(r2, 6),
            "tau_analytic_lin_h": round(tau_lin_h, 3),
            "analytic_vs_median_pct": round((tau_lin_h - med) / med * 100.0, 2) if med else None,
            "gap_at_4h_mm": round(gap4, 3),
            "new_4h_over_tau": round(4.0 / med, 4) if med else None,
            "gap_at_4h_over_band": round(gap4 / 150.0, 3),
            "settled_at_h": next((r["t_h"] for r in rows if abs(r["gap_mm"]) < TOL_MM), None),
        }
        print("  tau: 逐段 %s | 中位 %.3f | 均值 %.3f±%.3f | logfit %.3f (R2=%.5f)"
              % (tv, med, mean, std, tau_fit, r2), flush=True)
        print("  解析 tau_lin=%.3f h (与中位差 %.2f%%) | 4h/τ=%.3f | 4h残差 %.1f mm"
              % (tau_lin_h, out["cases"]["dq_%02.0f" % dq]["analytic_vs_median_pct"],
                 4.0 / med, gap4), flush=True)

    # 幅值依赖检验：τ ∝ (Q0+ΔQ) ?
    print("\n=== 幅值依赖检验（τ ∝ Q0+ΔQ）===", flush=True)
    chk = []
    for k, c in sorted(out["cases"].items()):
        pred = out["cases"]["dq_15"]["tau_median_h"] * c["Qs_m3s"] / out["cases"]["dq_15"]["Qs_m3s"]
        chk.append({"dq": c["dq_m3s"], "tau_obs_h": c["tau_median_h"],
                    "tau_pred_from_15_h": round(pred, 3),
                    "dev_pct": round((c["tau_median_h"] - pred) / pred * 100.0, 2)})
        print("  ΔQ=%2.0f: 观测 %.3f h, 按 τ∝Qs 预测 %.3f h, 偏差 %+.2f%%"
              % (c["dq_m3s"], c["tau_median_h"], pred, chk[-1]["dev_pct"]), flush=True)
    out["amplitude_scaling_check"] = chk

    out["wall_s"] = round(time.time() - t0, 1)
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_tau_identify.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_s"]), flush=True)


if __name__ == "__main__":
    main()
