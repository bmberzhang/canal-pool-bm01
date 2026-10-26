# -*- coding: utf-8 -*-
"""
(W5) 空间平均算子的直接数值核验 —— 构造空间正弦初始条件，测全池平均响应
================================================================================
预注册判据 v2（W5_理论澄清与诚实边界.md）判据1：
  全池平均算子对空间正弦的响应 = sin²(πx)/(πx)，x=L/λ，残差应 < 1e-3（恒等式）。

本脚本**不演化水力学**（避免色散/边界污染），直接：
  1. 构造空间正弦 δh(x) = δh0·sin(2πx/λ)（x 在 [0,L]，用 ncell 个格心采样）；
  2. 计算离散均值 mean = (1/ncell)Σ δh(x_i)；
  3. 与解析 sin²(πx)/(πx) 对比，验证"离散均值 → 连续均值"的收敛；
  4. 跨 L ∈ {2000,6000,12000,18000} 验证响应只依赖 x=L/λ（不依赖 L 本身）；
  5. 同时给出"闸前单点"读到的 δh(L)=sin(2πL/λ)，展示 point 通道为何不是干净参考。

关键：这是**静态几何核验**，隔离掉 E-FD 里"边界注入 + 色散传播"的混杂，
验证"全池平均算子"这一环的数学本质。色散链（判据2）已用 E-FD 数据核验过。
"""
import json, math, time

def resp_exact(x):
    if abs(x) < 1e-9:
        return 0.0
    return math.sin(math.pi * x) ** 2 / (math.pi * x)

def resp_discrete(L, lam, ncell):
    dx = L / ncell
    s = 0.0
    for i in range(ncell):
        x = (i + 0.5) * dx  # 格心采样
        s += math.sin(2 * math.pi * x / lam)
    return s / ncell  # 离散均值（相对振幅1）

def main():
    t0 = time.time()
    Ls = [2000.0, 6000.0, 12000.0, 18000.0]
    xs = [0.02, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5]
    rows = []
    print("=" * 78)
    print("W5 判据1：全池平均算子对空间正弦的响应 = sin²(πx)/(πx)（静态核验）")
    print("=" * 78)
    print(f"{'L(m)':>7} {'x=L/λ':>7} {'离散均值':>10} {'解析sin²/(πx)':>14} {'残差':>10}")
    max_resid = 0.0
    for L in Ls:
        for x in xs:
            lam = L / x
            disc = resp_discrete(L, lam, 400)  # 400 格足够细
            exact = resp_exact(x)
            resid = abs(disc - exact)
            max_resid = max(max_resid, resid)
            rows.append({"L_m": L, "x": x, "lambda_m": round(lam,2),
                         "discrete_mean": round(disc,6), "exact": round(exact,6),
                         "abs_resid": round(resid,6)})
            print(f"{L:>7.0f} {x:>7.2f} {disc:>10.5f} {exact:>14.5f} {resid:>10.6f}")

    print(f"\n最大残差 = {max_resid:.2e}  (判据1阈值 < 1e-3)")
    print(f"判据1：{'✅ 通过（恒等式精确成立）' if max_resid < 1e-3 else '❌ 失败'}")

    # 判据3：跨 L 响应只依赖 x（同一 x 不同 L 的响应应完全一致）
    print("\n判据3（标度）：同一 x 下，不同 L 的响应是否一致（应一致，因只依赖 x）")
    by_x = {}
    for r in rows:
        by_x.setdefault(r["x"], []).append(r["discrete_mean"])
    spread_max = 0.0
    for x, vals in sorted(by_x.items()):
        spread = (max(vals) - min(vals)) / (abs(sum(vals)/len(vals)) + 1e-12)
        spread_max = max(spread_max, spread)
    print(f"跨 L 最大相对离散 = {spread_max:.2e}  (应≈0，因为响应只依赖 x=L/λ)")
    print(f"判据3：{'✅ 平凡成立（L 是唯一长度尺度的必然）' if spread_max < 1e-6 else '⚠️ 需检查'}")

    # 展示 point 通道为何不是干净参考
    print("\n判据2 附带：'闸前单点' δh(L)=sin(2πL/λ)=sin(2πx) 会过零、变号，")
    print("故'池均/单点'不是测空间平均的干净对象（呼应诊断结论）。")

    out = {"purpose": "W5 判据1：空间平均算子静态核验",
           "max_abs_resid": max_resid, "ncell_sample": 400,
           "judgement1": "PASS" if max_resid < 1e-3 else "FAIL",
           "cross_L_spread": spread_max,
           "rows": rows, "wall_s": round(time.time()-t0,1)}
    json.dump(out, open("raw/bm01_w5_kc_static.json","w",encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n已写入 raw/bm01_w5_kc_static.json")

if __name__ == "__main__":
    main()
