#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M3：±150 mm 判定阈值敏感性分析（零重跑，读已记录峰值重算判定码）。

要回答的审稿质疑：A1/A3 全部依赖 150mm 阈值，而 gate_avg 档谷底 154.23mm 只比
阈值高 2.8%，结论会不会是阈值人为制造的？

方法：峰值 peak_mean_mm 是 band 无关的物理量，判定码只是 peak≤band?A1:A3 的后处理。
故读 bm01_tune_valley.json（L=6000, ΔQ=45, ncell=240, mean/gate_avg 各 72 格）的
已记录峰值，在 {100,125,150,175,200} mm 五档阈值上重算判定码，看翻转是否稳健。

产物：raw/bm01_m3_band_sensitivity.json
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")

BANDS_MM = [100, 125, 150, 175, 200]
REF_BAND = 150

def load():
    with open(os.path.join(RAW, "bm01_tune_valley.json"), encoding="utf-8") as f:
        return json.load(f)

def collect_peaks(d):
    """返回 {arm: [peak_mean_mm,...]}，并附原始 judgement。"""
    out = {}
    for arm, m in d["maps"].items():
        peaks = []
        for kp, kis in m["grid"].items():
            for ki, cell in kis.items():
                pk = cell.get("peak_mean_mm")
                if pk is not None:
                    peaks.append((pk, cell.get("judgement")))
        out[arm] = peaks
    return out

def a1_count(peaks, band):
    """peaks: [(peak,judgement)] → A1 格数。"""
    return sum(1 for pk, _ in peaks if pk <= band)

def main():
    d = load()
    peaks = collect_peaks(d)

    rows = []
    for arm in ["mean", "gate_avg"]:
        ps = [pk for pk, _ in peaks[arm]]
        n = len(ps)
        mn, mx = min(ps), max(ps)
        # 各阈值下 A1 格数 + 谷底格（最小峰值）是否进带
        a1_by_band = {}
        for b in BANDS_MM:
            a1_by_band[b] = a1_count(peaks[arm], b)
        # 该档最小峰值对应的"翻转阈值"（最小峰值首次进 A1 的 band）
        flip_band = None
        for b in sorted(BANDS_MM):
            if mn <= b:
                flip_band = b
                break
        rows.append({
            "arm": arm,
            "n_cells": n,
            "min_peak_mm": round(mn, 2),
            "max_peak_mm": round(mx, 2),
            "A1_count_by_band": {str(b): a1_by_band[b] for b in BANDS_MM},
            "A1_ratio_by_band": {str(b): round(a1_by_band[b] / n, 3) for b in BANDS_MM},
            "min_peak_first_A1_at_band_mm": flip_band,
        })

    # 翻转稳健性：在每档阈值下，mean 与 gate_avg 的 A1 占比差
    flip_gap_by_band = {}
    m = {r["arm"]: r["A1_count_by_band"] for r in rows}
    for b in BANDS_MM:
        bs = str(b)
        flip_gap_by_band[bs] = round(m["mean"][bs] / rows[0]["n_cells"] - m["gate_avg"][bs] / rows[1]["n_cells"], 3)

    out = {
        "case_id": "M3",
        "purpose": "band_sensitivity_of_sampling_reversal",
        "question": "A1/A3 是否因 150mm 阈值人为制造？翻转在 100/125/150/175/200mm 下是否稳健？",
        "method": "零重跑。峰值 peak_mean_mm 是 band 无关物理量，判定码为 peak≤band?A1:A3 的后处理。读 bm01_tune_valley.json（L=6000/ΔQ=45/ncell=240，mean 与 gate_avg 各 72 格）已记录峰值重算。",
        "evidence_grade": "synthetic-1d-sve",
        "bands_mm": BANDS_MM,
        "reference_band_mm": REF_BAND,
        "arms": rows,
        "A1_ratio_gap_mean_minus_gateavg_by_band": flip_gap_by_band,
    }

    p = os.path.join(RAW, "bm01_m3_band_sensitivity.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    # 打印结果表
    print("M3 阈值敏感性（L=6000, ΔQ=45, ncell=240）")
    print()
    print("各反馈档：峰值范围 与 各阈值下 A1 格数")
    print(f"{'反馈档':<12}{'格数':>4}{'最小峰值':>10}{'最大峰值':>10}" +
          "".join(f"{'|'+str(b)+'mm':>9}" for b in BANDS_MM))
    print("-" * 75)
    for r in rows:
        print(f"{r['arm']:<12}{r['n_cells']:>4}{r['min_peak_mm']:>10.2f}{r['max_peak_mm']:>10.2f}" +
              "".join(f"{r['A1_count_by_band'][str(b)]:>9}" for b in BANDS_MM))
    print()
    print("翻转的 A1 占比差（mean − gate_avg），越大说明翻转越稳健：")
    for b in BANDS_MM:
        print(f"  {b:>4}mm: {flip_gap_by_band[str(b)]:.3f}")
    print()
    print("关键结论：")
    for r in rows:
        print(f"  {r['arm']}: 最小峰值 {r['min_peak_mm']}mm，首次进 A1 需要条带 ≥ {r['min_peak_first_A1_at_band_mm']}mm")
    print()
    print("写出:", p)

if __name__ == "__main__":
    main()
