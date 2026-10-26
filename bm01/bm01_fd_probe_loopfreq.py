# -*- coding: utf-8 -*-
"""
(E-FD-S) 补充探针：在**实测极限环周期附近**（T = 5000 / 6000 / 7500 / 9000 s）
再测三通道频响，检验 E-FD Part B 留下的唯一疑点——
"环路实际振荡的那个频率上，取样尺度到底有没有造成通道差异？"

（Part B 只探到 1200–14400 s 的稀疏 5 点，且在 5400–8800 s 区间没有采样点。）

输出：raw/bm01_fd_probe_loopfreq.json
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bm01_fd_mech as FD

PERIODS = [5000.0, 6000.0, 7500.0, 9000.0]


def main():
    t0 = time.time()
    out = {"purpose": "E-FD-S 极限环周期附近的三通道频响补充探针",
           "periods_s": PERIODS, "ncell": 240, "results": {}}
    for T in PERIODS:
        t1 = time.time()
        r = FD.probe_open_loop(240, T)
        c = r["channels"]
        ra = c["avg100"]["amp_mm"] / max(c["mean"]["amp_mm"], 1e-9)
        rp = c["point"]["amp_mm"] / max(c["mean"]["amp_mm"], 1e-9)
        # 与池均的相位差（约定无关的相对量）
        pa = c["avg100"]["phase_deg"] - c["mean"]["phase_deg"]
        pp = c["point"]["phase_deg"] - c["mean"]["phase_deg"]
        out["results"]["T%.0f" % T] = {
            "channels": c,
            "ratio_avg100_over_mean": round(ra, 3),
            "ratio_point_over_mean": round(rp, 3),
            "phase_diff_avg100_minus_mean_deg": round(pa, 1),
            "phase_diff_point_minus_mean_deg": round(pp, 1),
            "wall_s": round(time.time() - t1, 1)}
        print("[FD-S] T=%6.0fs mean=%8.3f avg100=%8.3f point=%8.3f mm | r100=%5.2f rpt=%5.2f | dphi100=%6.1f dphipt=%6.1f (%.0fs)"
              % (T, c["mean"]["amp_mm"], c["avg100"]["amp_mm"], c["point"]["amp_mm"],
                 ra, rp, pa, pp, time.time() - t1), flush=True)
    out["wall_total_s"] = round(time.time() - t0, 1)
    os.makedirs("raw", exist_ok=True)
    jp = os.path.join("raw", "bm01_fd_probe_loopfreq.json")
    json.dump(out, open(jp, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("\n已写入 %s（%.0fs）" % (jp, out["wall_total_s"]))


if __name__ == "__main__":
    main()
