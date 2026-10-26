# -*- coding: utf-8 -*-
"""
(W5) 空间低通截止波数 k_c 的显式表达 —— 第一步：用已有 E-FD 数据检验 sinc(L/λ) 闭式
================================================================================
预注册判据见 ../../WorkBuddy/2026-10-11-15-20-33/W5_预注册判据.md（先于本脚本写死）。

本脚本零新算例，只读已有 E-FD 频响数据，检验：
  判据1（形状吻合）：实测"池均/单点增益比"随 λ 变化，能否被 sinc(L/λ) 无自由参数拟合（R²≥0.95）；
  判据2（截止位置）：半功率点 λ 与 λ_c = 1.6587·L 偏差 ≤20%。

物理模型：
  水深扰动 δh(x) = δh0·sin(2πx/λ)。全池平均 = (1/L)∫0^L δh dx。
  单点（闸前 x≈L）读的是 δh(L) = δh0·sin(2πL/λ)。
  于是"池均振幅 / 单点振幅"的理论比 = |sinc(L/λ)| / |sin(2πL/λ)|。
  但单点相位随 λ 变化、且单点本身是 sin(2πL/λ) 会过零，直接用单点做参考不稳。

  更稳的参考：单点与闸前 100m 平均在 E-FD 里几乎重合（比值 0.97~1.0），
  且二者都接近"局部未滤波"的通道。因此用"池均 amp / 单点 amp"实测比 G_meas，
  与理论 |sinc(L/λ)| 比较（单点在 λ≫L 时 sin→sin(2πL/λ)→0 的过零需避开，
  但 E-FD 的频点 T≥1200 s 对应 λ = c·T ≥ 8.8·1200 ≈ 10560 m ≈ 1.76L，sin(2πL/λ) 尚在正值区）。

关键换算：λ = c·T，c = √(g·h) + Q0/(B·h) = 8.801 m/s（E-FD reference_time_scales 已给）。
  T=1200 s → λ≈10561 m → L/λ≈0.568 → sinc(0.568)=sin(1.785)/(1.785)=0.976/1.785=0.547
  实测 G_meas(T=1200) = mean/point = 2.303/8.971 = 0.257  ← 若理论成立应≈0.547，差了 2 倍？

  这说明"池均/单点"不是纯 sinc，因为单点通道本身在闸门处被边界条件放大。
  真正无歧义的理论对象是"全池平均算子对空间正弦的响应"，即 H(λ)=|sinc(L/λ)|，
  它应该对照"池均通道的绝对增益如何随 λ 减小"——但绝对增益还含整体水力学传递函数，
  不能直接比。

  因此第一步的判据1/2 需要更仔细的归一化。本脚本的职责是：
  (a) 把 E-FD 三通道数据按 λ 重排，输出"池均/单点"与"池均/avg100"两条实测比值曲线；
  (b) 与 sinc(L/λ) 叠图，量化偏差；
  (c) 若"池均/单点"比值本身不贴合 sinc，则如实报告：单点通道含边界放大，需改用
      "两个不同平均窗长"之间的比值（如 avg100 vs mean）来隔离纯空间平均效应。

诚实原则：本脚本只报告实测与理论的吻合/偏差，不改判据、不事后挑数据。
"""
import json, math, os

L = 6000.0
C = 8.801  # m/s，E-FD reference_time_scales.celerity_m_s

def sinc(x):
    # sinc(x) = sin(pi x)/(pi x)
    if abs(x) < 1e-9:
        return 1.0
    return math.sin(math.pi * x) / (math.pi * x)

def main():
    d = json.load(open("raw/bm01_fd_mech.json", encoding="utf-8"))
    fb = d["partB_frequency_response"]
    rows = []
    for key, blk in fb.items():
        T = blk["period_s"]
        lam = C * T
        ch = blk["channels"]
        mean = ch["mean"]["amp_mm"]; avg = ch["avg100"]["amp_mm"]; pt = ch["point"]["amp_mm"]
        if not (mean and avg and pt):
            continue
        rows.append({
            "T_s": T, "lambda_m": round(lam,1), "L_over_lambda": round(L/lam, 4),
            "mean_mm": mean, "avg100_mm": avg, "point_mm": pt,
            "ratio_mean_over_point": round(mean/pt, 4),
            "ratio_mean_over_avg100": round(mean/avg, 4),
            "sinc_L_over_lambda": round(sinc(L/lam), 4),
        })
    # 补充 E-FD-S 的极限环频段四频点
    d2 = json.load(open("raw/bm01_fd_probe_loopfreq.json", encoding="utf-8"))
    for key, blk in d2["results"].items():
        T = blk["channels"]["mean"]  # 注意结构不同，这里 period 从 key 取
    # 上面结构不对，单独重读
    rows2 = []
    for key, blk in d2["results"].items():
        T = float(key[1:])
        lam = C * T
        ch = blk["channels"]
        mean = ch["mean"]["amp_mm"]; avg = ch["avg100"]["amp_mm"]; pt = ch["point"]["amp_mm"]
        rows2.append({
            "T_s": T, "lambda_m": round(lam,1), "L_over_lambda": round(L/lam, 4),
            "mean_mm": mean, "avg100_mm": avg, "point_mm": pt,
            "ratio_mean_over_point": round(mean/pt, 4),
            "ratio_mean_over_avg100": round(mean/avg, 4),
            "sinc_L_over_lambda": round(sinc(L/lam), 4),
        })
    rows.extend(rows2)
    rows.sort(key=lambda r: r["T_s"])

    print("=" * 80)
    print("W5 第一步：E-FD 已有数据 vs sinc(L/λ) 闭式（零新算例）")
    print("=" * 80)
    print(f"{'T(s)':>7} {'λ(m)':>8} {'L/λ':>7} {'mean':>8} {'avg100':>8} {'point':>8} "
          f"{'mean/pt':>8} {'mean/avg':>9} {'sinc':>8}")
    for r in rows:
        print(f"{r['T_s']:>7.0f} {r['lambda_m']:>8.0f} {r['L_over_lambda']:>7.3f} "
              f"{r['mean_mm']:>8.2f} {r['avg100_mm']:>8.2f} {r['point_mm']:>8.2f} "
              f"{r['ratio_mean_over_point']:>8.3f} {r['ratio_mean_over_avg100']:>9.3f} "
              f"{r['sinc_L_over_lambda']:>8.3f}")

    # 关键诊断：mean/point 比值 vs sinc
    print("\n--- 诊断：实测 mean/point vs 理论 sinc(L/λ) ---")
    for r in rows:
        print(f"  L/λ={r['L_over_lambda']:.3f}  mean/pt={r['ratio_mean_over_point']:.3f}  "
              f"sinc={r['sinc_L_over_lambda']:.3f}  比值(meas/theory)="
              f"{r['ratio_mean_over_point']/max(r['sinc_L_over_lambda'],1e-9):.3f}")

    json.dump({"purpose": "W5 第一步：E-FD 已有数据 vs sinc 闭式",
               "rows": rows}, open("raw/bm01_w5_step1_diagnosis.json","w",encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print("\n已写入 raw/bm01_w5_step1_diagnosis.json")

if __name__ == "__main__":
    main()
