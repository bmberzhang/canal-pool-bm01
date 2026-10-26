"""BM-01 单池｜「前馈收益」随输运时滞 τ 的变化曲线。

回答一个问题：前馈到底在 τ 多大时才值得上？
单一变量扫描：只动 τ（= L/c，物理上等价于波速可调），池子面积 A 固定不变。
窗口 W 用绝对网格扫，不按 τ 的比例扫，避免把结论写进网格里。

自检：τ = 682 s（BM-01 动力波速口径）必须复现已知三个锚点。
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bm01_pool_step import Params, build_scenarios, judge          # noqa: E402
import bm01_inflow_integral_ff as FF                              # noqa: E402
import bm01_single_pool_ff as SP                                  # noqa: E402

TAU_GRID = [0.0, 60.0, 120.0, 240.0, 300.0, 375.0, 450.0, 600.0,
            682.0, 900.0, 1200.0, 1800.0, 2400.0, 3600.0]
W_GRID = [60.0, 150.0, 300.0, 450.0, 600.0, 900.0, 1200.0, 1800.0, 2400.0, 3600.0]


def sweep(base: Params, q: np.ndarray, sc) -> dict:
    out = {}
    for tau in TAU_GRID:
        row = {}
        r = FF.run(base, q, q, tau, "off")
        row["off"] = FF.stats(base, r["series"]["h1_dev"], r["series"]["a1"])

        r = FF.run(base, q, q, tau, "inst", gain=1.0)
        row["inst"] = FF.stats(base, r["series"]["h1_dev"], r["series"]["a1"])
        row["t90_inst"] = FF.t_reach_90(base, r["series"]["ff_open"], FF.STEP_DQ)

        best = None
        for W in W_GRID:
            r = FF.run(base, q, q, tau, "win", W=W, gain=1.0)
            s = FF.stats(base, r["series"]["h1_dev"], r["series"]["a1"])
            t90 = FF.t_reach_90(base, r["series"]["ff_open"], FF.STEP_DQ)
            if best is None or s["peak"] < best[1]["peak"]:
                best = (W, s, t90, r)
        row["W"] = best[0]
        row["win"] = best[1]
        row["t90_win"] = best[2]
        row["win_r"] = best[3]
        out[tau] = row
    return out


def main() -> int:
    # 阶跃工况的扰动时刻正好踩在控制节拍上，会掩盖相位效应；
    # real 用未知时刻的多阶跃，用来检验曲线上那个尖峰是不是采样巧合。
    kind = sys.argv[1] if len(sys.argv) > 1 else "step"
    base = Params()
    sc = next(s for s in build_scenarios(base) if s.key == "main")
    tau0 = SP.id_tau(base, "dynamic")
    q, _ = FF.make_inflow(base, kind, noise=(kind == "real"))
    Tc = base.ctrl_period_s
    t_act = (FF.STEP_DQ / base.kg) / base.rate_limit_per_s

    print("=" * 104)
    print("BM-01 单池｜前馈收益随输运时滞 τ 的变化（来水工况 = %s）" % kind)
    print("=" * 104)
    print("固定量：池面积 A = %.1e m2（%.0f m × %.0f m），控制周期 T_c = %.0f s，"
          "闸门开到位 Δa/rate = %.0f s"
          % (base.area_m2, base.pool_len_m, base.pool_width_m, Tc, t_act))
    print("单一变量：只改 τ = L/c（等价于波速可调）。τ 的参考值 %.0f s 来自动力波速口径。"
          % tau0)
    print()

    # ---------------- 自检（仅阶跃工况有已知锚点）----------------
    if kind == "step":
        print("自检（τ = %.0f s 必须复现已知锚点）：" % tau0)
        ref = FF.run(base, q, q, tau0, "off")
        p_off = FF.stats(base, ref["series"]["h1_dev"], ref["series"]["a1"])["peak"]
        ref = FF.run(base, q, q, tau0, "inst", gain=1.0)
        p_inst = FF.stats(base, ref["series"]["h1_dev"], ref["series"]["a1"])["peak"]
        ref = FF.run(base, q, q, tau0, "win", W=tau0, gain=1.0)
        p_win = FF.stats(base, ref["series"]["h1_dev"], ref["series"]["a1"])["peak"]
        print("  A1 纯反馈 |峰| = %.2f mm（锚点 66.47）→ %s"
              % (p_off, "OK" if abs(p_off - 66.47) < 0.02 else "**不一致**"))
        print("  A2 瞬时前馈 |峰| = %.2f mm（锚点 26.25）→ %s"
              % (p_inst, "OK" if abs(p_inst - 26.25) < 0.02 else "**不一致**"))
        print("  A3 窗口 W=τ |峰| = %.2f mm（锚点 17.84）→ %s"
              % (p_win, "OK" if abs(p_win - 17.84) < 0.05 else "**不一致**"))
        assert abs(p_off - 66.47) < 0.02 and abs(p_inst - 26.25) < 0.02
        print()
    else:
        print("(real 工况无固定锚点，自检见 step 模式)")
        print()

    rows = sweep(base, q, sc)

    print("=" * 104)
    print("表 1  前馈收益 vs τ")
    print("=" * 104)
    print("%7s %8s | %9s %9s %7s | %7s %9s %7s | %9s %9s"
          % ("τ(s)", "τ/T_c", "纯反馈mm", "瞬时FFmm", "降幅", "最优W", "窗口FFmm", "降幅",
             "瞬时净提前", "窗口净提前"))
    print("-" * 104)
    for tau in TAU_GRID:
        r = rows[tau]
        d_inst = 100 * (1 - r["inst"]["peak"] / r["off"]["peak"])
        d_win = 100 * (1 - r["win"]["peak"] / r["off"]["peak"])
        # 净提前量（正 = 命令到位早于水到池尾）
        lead_i = tau - (r["t90_inst"] - base.t_step_s)
        lead_w = tau - (r["t90_win"] - base.t_step_s)
        print("%7.0f %8.2f | %9.1f %9.1f %6.0f%% | %7.0f %9.1f %6.0f%% | %9.0f %9.0f"
              % (tau, tau / Tc, r["off"]["peak"], r["inst"]["peak"], d_inst,
                 r["W"], r["win"]["peak"], d_win, lead_i, lead_w))
    print("-" * 104)
    print("「降幅」= 相对同 τ 下纯反馈的降幅；「净提前」= τ −（命令涨到 90%% 所用时间），正=早于水到。")
    print()

    # ---------------- 表 2：门槛定位 ----------------
    print("=" * 104)
    print("表 2  门槛定位：前馈降幅跨过 20% / 50% 的 τ")
    print("=" * 104)
    for thr in (20.0, 50.0):
        hit = [tau for tau in TAU_GRID
               if 100 * (1 - rows[tau]["win"]["peak"] / rows[tau]["off"]["peak"]) >= thr]
        print("  窗口积分降幅 ≥ %.0f%%：%s"
              % (thr, ("τ ≥ %.0f s" % min(hit)) if hit else "整条网格都没有"))
    for thr in (20.0,):
        hit = [tau for tau in TAU_GRID
               if 100 * (1 - rows[tau]["inst"]["peak"] / rows[tau]["off"]["peak"]) >= thr]
        print("  瞬时前馈降幅 ≥ %.0f%%：%s"
              % (thr, ("τ ≥ %.0f s" % min(hit)) if hit else "整条网格都没有"))
    print()
    print("  对照门槛：控制周期 T_c = %.0f s；闸门开到位 t_act = %.0f s。" % (Tc, t_act))
    print()

    # ---------------- 表 3：判定等级有没有变 ----------------
    print("=" * 104)
    print("表 3  判定等级（正常条带 ±%.0f mm）" % (base.band_normal_m * 1000))
    print("=" * 104)
    print("%7s %12s %12s %12s" % ("τ(s)", "纯反馈", "瞬时FF", "窗口FF"))
    print("-" * 50)
    for tau in TAU_GRID:
        r = rows[tau]
        codes = []
        for key in ("off", "inst"):
            rr = FF.run(base, q, q, tau, {"off": "off", "inst": "inst"}[key],
                        gain=1.0 if key == "inst" else 1.0)
            codes.append(judge(base, sc, rr)["code"])
        codes.append(judge(base, sc, r["win_r"])["code"])
        print("%7.0f %12s %12s %12s" % (tau, codes[0], codes[1], codes[2]))
    print("-" * 50)
    print()

    # ---------------- 存曲线 ----------------
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "raw", "bm01_ff_tau_map.npz")
    np.savez(out,
             tau=np.array(TAU_GRID),
             off=np.array([rows[t]["off"]["peak"] for t in TAU_GRID]),
             inst=np.array([rows[t]["inst"]["peak"] for t in TAU_GRID]),
             win=np.array([rows[t]["win"]["peak"] for t in TAU_GRID]),
             W=np.array([rows[t]["W"] for t in TAU_GRID]),
             t90_inst=np.array([rows[t]["t90_inst"] for t in TAU_GRID]),
             t90_win=np.array([rows[t]["t90_win"] for t in TAU_GRID]))
    print("曲线数据已存：raw/bm01_ff_tau_map.npz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
