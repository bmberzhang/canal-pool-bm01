# -*- coding: utf-8 -*-
"""
(W4'') 两型失效假设的预注册检验：跃升倍数是否随 T_ctrl/(L/c) 单调
================================================================================
三模型第三轮 C 判定：这是"唯一还抬机制证据上限的生死活"。要求**判据先预注册**再扫描。

=================== 预注册（先于任何计算，不得事后修改）===================
对象：BM-01 系列（矩形 B=50 m、h_t=7 m、q0=180、k_g=360、速率限 1/3000、kp=1.6、ki=1.6e-3）
波速取 c_eff = sqrt(g·h_t) + v0 = sqrt(9.81*7) + 180/(50*7) = 8.287 + 0.514 = 8.801 m/s
比值        rho = T_ctrl / (L / c_eff)

度量：跃升倍数  R = (起爆档的全池平均峰值) / (起爆前一档的全池平均峰值)，
      起爆档 = 使"三档全 A3"成立的最小 ΔQ；前一档 = 该 ΔQ 减去一个搜索步长 δ。
      因实测显示"三档全 A3" ⟺ 全池平均档越带（点/100m 档早已越带），
      故用**全池平均档**对 ΔQ 做二分以定位起爆点（步长内 ±2 m³/s），
      定位后再用**三档**复核一次。

预注册判据（双分支，二选一）：
  分支 A（两型假设存活）：在 rho ∈ [0.4, 1.8] 内取 ≥4 个水平，
       R 随 rho 单调不降（Spearman 秩相关 ρ_s ≥ +0.7，且两端 R 相差 ≥3 倍）。
  分支 B（砍掉两型）：不满足分支 A ⇒ 撤销"失稳型/边际型"两型叙事，
       退回单一机制（"控制环路设计限制"），不得在正文区分两型。

输出：raw/bm01_w4pp_two_type.json
=========================================================================
"""
import json, os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_sve1d_v2 as V2

BAND_MM = 150.0
B0, HT, Q0, KG, RATE = 50.0, 7.0, 180.0, 360.0, 1.0 / 3000.0
C_EFF = math.sqrt(9.81 * HT) + Q0 / (B0 * HT)
RHO_LO, RHO_HI = 0.4, 1.8
DELTA = 2.0
# (L, T_ctrl) 组合：按 rho = T/(L/c) 覆盖 [0.4, 1.8]
COMBOS = [(2000.0, 100.0), (2000.0, 150.0), (2000.0, 200.0), (2000.0, 250.0), (2000.0, 300.0),
          (6000.0, 300.0), (6000.0, 450.0), (6000.0, 600.0), (6000.0, 800.0), (6000.0, 1000.0),
          (6000.0, 1200.0)]
SEARCH_LO, SEARCH_HI = 6.0, 120.0
ARMS3 = (("point", dict(ctrl_on="gate")),
         ("avg100", dict(ctrl_on="gate_avg", sensor_len_m=100.0)),
         ("mean", dict(ctrl_on="mean")))


def base_p(L, period):
    return {"pool_len_m": L, "pool_width_m": B0, "q0": Q0, "h_target": HT, "a0": 0.5,
            "kg": KG, "kp": 1.6, "ki": 0.0016, "ctrl_period_s": period, "deadband_m": 0.01,
            "meas_resolution_m": 0.001, "rate_limit_per_s": RATE, "t_step_s": 900.0,
            "band_normal_m": 0.15, "band_alarm_m": 0.3, "settle_tol_m": 0.02,
            "settle_window_s": 1800.0, "duration_s": 43200.0}


def mean_peak(L, ncell, period, dq):
    r = V2.sve_1d_v2(base_p(L, period), ncell, dq, m_slope=0.0, ctrl_on="mean")
    return r["peak_mean_mm"]


def all_three(L, ncell, period, dq):
    arms = {}
    for name, kw in ARMS3:
        r = V2.sve_1d_v2(base_p(L, period), ncell, dq, m_slope=0.0, **kw)
        arms[name] = {"peak_mean_mm": round(r["peak_mean_mm"], 2),
                      "judgement": "A1" if r["peak_mean_mm"] <= BAND_MM else "A3"}
    n = sum(1 for a in ARMS3 if arms[a[0]]["judgement"] == "A3")
    return arms, n == 3


def find_onset(L, ncell, period):
    """二分找最小使全池平均越带的 ΔQ（±1）"""
    lo, hi = SEARCH_LO, SEARCH_HI
    if mean_peak(L, ncell, period, hi) <= BAND_MM:
        return None
    if mean_peak(L, ncell, period, lo) > BAND_MM:
        return lo
    while hi - lo > 1.0:
        mid = 0.5 * (lo + hi)
        if mean_peak(L, ncell, period, mid) > BAND_MM:
            hi = mid
        else:
            lo = mid
    return hi


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1.0
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx * dy > 0 else 0.0


def main():
    t0 = time.time()
    out = {"purpose": "W4'' 预注册检验：跃升倍数 vs rho = T_ctrl/(L/c_eff)",
           "preregistered": {"rho_range": [RHO_LO, RHO_HI], "delta_m3s": DELTA,
                             "branch_A": "rho 取 ≥4 水平、R 随 rho 单调不降（Spearman ρs ≥ +0.7 且两端 R 相差 ≥3 倍）→ 两型存活",
                             "branch_B": "不满足 A ⇒ 砍掉两型，退回单一机制"},
           "c_eff_m_s": round(C_EFF, 3), "rows": [], "wall_total_s": None}

    for L, period in COMBOS:
        ncell = int(round(L / 25.0))
        rho = period / (L / C_EFF)
        onset = find_onset(L, ncell, period)
        if onset is None:
            print("[W4''] L=%5.0f T=%5.0fs rho=%.2f | ΔQ≤%.0f 内无起爆" % (L, period, rho, SEARCH_HI),
                  flush=True)
            out["rows"].append({"L_m": L, "ctrl_period_s": period, "rho": round(rho, 3),
                                "onset_dq": None, "note": "区间内无起爆"})
            continue
        prev = max(SEARCH_LO, onset - DELTA)
        mp_on = mean_peak(L, ncell, period, onset)
        mp_prev = mean_peak(L, ncell, period, prev)
        arms_on, all3_on = all_three(L, ncell, period, onset)
        R = mp_on / mp_prev if mp_prev > 0 else float('inf')
        out["rows"].append({"L_m": L, "ctrl_period_s": period, "rho": round(rho, 3),
                            "onset_dq": round(onset, 2), "prev_dq": round(prev, 2),
                            "mean_peak_onset_mm": round(mp_on, 2),
                            "mean_peak_prev_mm": round(mp_prev, 2),
                            "jump_ratio_R": round(R, 2),
                            "arms_at_onset": arms_on, "all_three_A3_at_onset": all3_on,
                            "L_over_c_s": round(L / C_EFF, 1)})
        print("[W4''] L=%5.0f T=%5.0fs rho=%.2f | 起爆 ΔQ=%.1f（前档 %.1f：%.1f → %.1f mm）"
              " | 跃升 R=%.2f | 三档全A3=%s"
              % (L, period, rho, onset, prev, mp_prev, mp_on, R, all3_on), flush=True)

    ok = [r for r in out["rows"] if r.get("jump_ratio_R")]
    ok.sort(key=lambda r: r["rho"])
    out["n_levels"] = len(ok)
    if len(ok) >= 4:
        rhos = [r["rho"] for r in ok]
        Rs = [r["jump_ratio_R"] for r in ok]
        rho_s = spearman(rhos, Rs)
        ratio_ends = max(Rs) / min(Rs) if min(Rs) > 0 else float('inf')
        branch_A = (rho_s >= 0.7 and ratio_ends >= 3.0)
        out["spearman_rho"] = round(rho_s, 3)
        out["R_max_over_R_min"] = round(ratio_ends, 2)
        out["branch"] = "A（两型假设存活）" if branch_A else "B（砍掉两型）"
        out["rhos"] = rhos
        out["Rs"] = Rs
        out["verdict"] = (
            ("**分支 A：两型假设存活。** 在 rho ∈ [%.2f, %.2f] 的 %d 个水平上，"
             "跃升倍数 R 随 rho 单调不降（Spearman ρs = %+.3f ≥ +0.7），两端相差 %.1f 倍（≥3 倍）"
             "⇒ 可保留【失稳型（R 大）/ 边际型（R 小）】两型，并报告 rho 为该二分的控制参数"
             "（n=%d 水平，非 n=2）。" % (min(rhos), max(rhos), len(ok), rho_s, ratio_ends, len(ok)))
            if branch_A else
            ("**分支 B：砍掉两型。** 跃升倍数 R 随 rho 的 Spearman ρs = %+.3f（需 ≥ +0.7）、"
             "两端相差 %.1f 倍（需 ≥3 倍）⇒ 不满足预注册判据，**撤销失稳型/边际型之分**，"
             "退回单一机制【控制环路设计限制】，正文不得区分两型。"
             % (rho_s, ratio_ends)))
    else:
        out["branch"] = "样本不足"
        out["verdict"] = ("有效水平 %d 个（<4），**样本不足** ⇒ 按预注册的分支 B 处理："
                          "撤销两型叙事，退回单一机制。" % len(ok))
    out["wall_total_s"] = round(time.time() - t0, 1)

    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_w4pp_two_type.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n" + out["verdict"])
    print("已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
