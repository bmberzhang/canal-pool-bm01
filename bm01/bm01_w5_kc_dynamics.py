# -*- coding: utf-8 -*-
"""
(W5) 判据2 动力学核验：初始空间正弦扰动 → 全池平均通道的响应
================================================================================
预注册判据 v2（W5_理论澄清与诚实边界.md）判据2：
  给定空间波长 λ（经色散 λ=c/f 来自执行器时间尺度），全池平均对 sin(2πx/λ)
  的响应 = sin²(πx)/(πx)，x=L/λ。本脚本做动力学侧的干净核验。

方法（隔离边界注入与色散混杂）：
  1. 初始水深 h(x) = ht + δh0·sin(2πx/λ)，δh0 小振幅（线性区）；
  2. 闸门固定（ctrl_on="off"，或直接固定 a=a0），上游 q_up=q0 固定；
  3. 在 t=0 立即测量"全池平均响应" mean(δh)/δh0；
  4. 与解析 sin²(πx)/(πx) 对比 —— 这是 t=0 的精确检验（无演化污染）；
  5. 再短演化（t < 0.2·L/c，避免边界反射）看响应如何随色散演化。

关键诚实点：t=0 的响应就是恒等式（与静态核验同），真正的物理检验是
"色散如何改变这个响应"——即把 E-FD 的"闸门注入→传播→池均"串起来。

本脚本聚焦 t=0 与短时演化的响应，输出与 sin²(πx)/(πx) 的偏差随时间的轨迹。
"""
import json, math, time
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_sve1d_core import Section, G, MANNING_N, pick_dt, normal_slope

L = 6000.0
B = 50.0
HT = 7.0
Q0 = 180.0
KG = 360.0
A0 = 0.5
C = math.sqrt(G * HT) + Q0 / (B * HT)  # ≈ 8.80 m/s

def resp_exact(x):
    if abs(x) < 1e-9:
        return 0.0
    return math.sin(math.pi * x) ** 2 / (math.pi * x)

def run_initial_sine(L, lam, delta_h0, ncell=400, n_steps=200, m_slope=0.0):
    """初始正弦扰动，闸门固定，上游固定，短时演化。返回 mean(δh)/δh0 的时序。"""
    sec = Section(B, m_slope)
    dx = L / ncell
    S0 = normal_slope({"pool_width_m": B, "h_target": HT, "q0": Q0}, MANNING_N)
    # 初始水深正弦
    h = [HT + delta_h0 * math.sin(2 * math.pi * (i + 0.5) * dx / lam) for i in range(ncell)]
    Q = [Q0] * ncell
    wmax = Q0 / (B * HT) + math.sqrt(G * HT)
    dt = pick_dt(dx, wmax)
    sqrt = math.sqrt; hmin = 1e-3; c1 = dt / dx
    fric = m_slope > 1e-12
    q_gate = KG * A0 * sqrt(HT) / sqrt(HT)  # = KG*A0 = 180
    trace = []
    for k in range(n_steps):
        tk = k * dt
        q_up = Q0
        hx = [h[0]] + h + [h[-1]]
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1
        F1 = [0.0] * nif; F2 = [0.0] * nif
        for j in range(nif):
            h1, h2 = max(hx[j], hmin), max(hx[j + 1], hmin)
            A1, A2 = max(sec.A(h1), 1e-6), max(sec.A(h2), 1e-6)
            B1, B2 = sec.B(h1), sec.B(h2)
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / B1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / B2)
            sm = s1 if s1 > s2 else s2
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2
                           + 0.5 * G * A1 * A1 / B1 + 0.5 * G * A2 * A2 / B2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate
        h_new = [0.0] * ncell; Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]; A = max(sec.A(hi), 1e-6); Qi = Q[i]
            R = A / sec.P(hi) if fric and sec.P(hi) > 0 else A / (B + 2.0 * hi)
            Sf = MANNING_N * MANNING_N * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)
        h, Q = h_new, Q_new
        dm = sum(h) / ncell - HT
        trace.append({"t_s": round(tk, 3), "mean_dh_m": dm, "resp": dm / delta_h0})
    return trace

def main():
    t0 = time.time()
    # 极限环波纹对应的空间波长（来自 E-FD Part A 实测 T_LC 经色散）
    T_lc_list = [5355.0, 5910.0, 7710.0, 8820.0]
    lam_list = [C * T for T in T_lc_list]
    xs = [L / lam for lam in lam_list]
    delta_h0 = 0.01  # 1 cm 小扰动，线性区

    print("=" * 78)
    print("W5 判据2：初始空间正弦 → 全池平均通道响应（动力学核验）")
    print("=" * 78)
    print(f"极限环周期 → 波长 → x=L/λ → 解析 sin²(πx)/(πx) → 实测(t=0)响应")
    results = []
    for T, lam, x in zip(T_lc_list, lam_list, xs):
        trace = run_initial_sine(L, lam, delta_h0, ncell=400, n_steps=50)
        resp_t0 = trace[0]["resp"]
        exact = resp_exact(x)
        # 短时演化（t < 0.2 L/c ≈ 136 s）后响应变化
        t_ref = 0.2 * L / C
        resp_later = next((r["resp"] for r in trace if r["t_s"] >= t_ref), trace[-1]["resp"])
        results.append({"T_LC_s": T, "lambda_m": round(lam,1), "x": round(x,4),
                        "exact": round(exact,4), "resp_t0": round(resp_t0,4),
                        "resp_later": round(resp_later,4)})
        print(f"  T={T:>5.0f}s λ={lam:>6.0f}m x={x:.4f} "
              f"exact={exact:.4f} t0={resp_t0:.4f} t≈{t_ref:.0f}s={resp_later:.4f}")

    # 判断 t=0 响应是否 = 恒等式（应精确，因为初始条件就是正弦）
    max_err_t0 = max(abs(r["exact"] - r["resp_t0"]) for r in results)
    print(f"\n判据2 t=0 最大偏差 = {max_err_t0:.2e}  (应≈0，因初始条件=正弦)")
    # 短时演化后的偏差（色散效应）
    max_err_later = max(abs(r["exact"] - r["resp_later"]) for r in results)
    print(f"判据2 短时演化(t≈136s)后最大偏差 = {max_err_later:.4f}  (色散如何偏离静态响应)")

    out = {"purpose": "W5 判据2：初始正弦→池均响应动力学核验",
           "delta_h0_m": delta_h0, "c_m_s": round(C,3),
           "results": results,
           "max_err_t0": max_err_t0, "max_err_later": max_err_later,
           "wall_s": round(time.time()-t0,1)}
    json.dump(out, open("raw/bm01_w5_kc_dynamics.json","w",encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n已写入 raw/bm01_w5_kc_dynamics.json")

if __name__ == "__main__":
    main()
