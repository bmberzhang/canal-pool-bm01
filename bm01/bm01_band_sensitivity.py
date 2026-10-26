#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BM-01 判定条带灵敏度扫描。

问题：±0.15 m 这个"正常调节条带"是 assume（量级示例），不是官方给定值。
     那它到底影不影响判定结论？本脚本用已记录的场景时序重算 R1，
     扫一遍条带宽度，看每个场景的判定码在哪个条带处翻转。

不重跑水力模型：直接读 raw/*.json 里已有的 series_downsampled 与 metrics。
R2（饱和）与 band 无关，沿用记录值；R3（末段收敛）也用记录值。

用法：python bm01_band_sensitivity.py
产物：raw/bm01_band_sensitivity.json
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")

# 扫描网格（米）。覆盖从"桌面级严条带"到"中线量级宽条带"
BANDS_M = [0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30]

SCENARIOS = [
    ("bm01_main", "主案例 单池", ["pool1"]),
    ("bm01_chain_local_only", "对照一 串联+本地反馈", ["pool1", "pool2"]),
    ("bm01_chain_feedforward", "对照二 串联+前馈", ["pool1", "pool2"]),
    ("bm01_overload", "对照三 超能力", ["pool1"]),
]


def load(name):
    with open(os.path.join(RAW, name + ".json"), encoding="utf-8") as fh:
        return json.load(fh)


def peak_from_series(d, chan):
    """从降采样时序取该池偏差绝对值最大值（米）。"""
    key = f"{chan}_dev_m"
    s = d["series_downsampled"].get(key)
    if s is None:
        return None
    return max(abs(v) for v in s)


def code_for_band(d, chans, band_m):
    """按 band 重算 R1，与记录的 R2/R3 合成判定码。返回 (code, r1_ok, worst_peak)。"""
    r1_ok = True
    worst = 0.0
    for c in chans:
        pk = peak_from_series(d, c)
        if pk is None:
            pk = d["metrics"][c]["peak_dev_m"]      # 回退到记录值
        worst = max(worst, pk)
        if pk > band_m:
            r1_ok = False
    # R2 / R3 与条带无关，直接取记录
    r2_ok = all(r["pass"] for r in d["judgement"]["rules"]
                if r["rule"] == "R2_no_sustained_saturation")
    r3_ok = all(r["pass"] for r in d["judgement"]["rules"]
                if r["rule"] == "R3_settled")
    if r1_ok and r2_ok and r3_ok:
        return "A1", True, worst
    sc_step = d["scenario"]["step_dq"]
    cap = d["judgement"]["absorb_capacity_m3s"]
    # 越带即 A3（记录里的两条 A3 分支只在措辞上不同）
    _ = sc_step, cap
    return "A3", False, worst


def main():
    rows = []
    for name, label, chans in SCENARIOS:
        d = load(name)
        rec_code = d["judgement"]["code"]
        rec_band = d["params"]["band_normal_m"]
        entry = {"key": name, "label": label, "recorded_code": rec_code,
                 "recorded_band_m": rec_band, "channels": chans, "sweep": []}
        for b in BANDS_M:
            code, r1_ok, worst = code_for_band(d, chans, b)
            entry["sweep"].append({
                "band_m": b, "band_mm": round(b * 1000, 1),
                "code": code, "r1_pass": r1_ok,
                "worst_peak_mm": round(worst * 1000, 2),
            })
        # 翻转点：最小能让该场景变 A1 的条带
        flip = next((s["band_mm"] for s in entry["sweep"] if s["code"] == "A1"), None)
        entry["min_band_mm_for_A1"] = flip
        entry["recorded_worst_peak_mm"] = round(
            max(peak_from_series(d, c) or d["metrics"][c]["peak_dev_m"] for c in chans) * 1000, 2)
        rows.append(entry)

    out = {
        "case_id": "BM-01",
        "purpose": "band_sensitivity",
        "question": "±0.15 m 正常调节条带是 assume；判定码对它有多敏感？",
        "method": "读已记录时序重算 R1；R2/R3 与条带无关，沿用记录值。未重跑水力模型。",
        "evidence_grade": "synthetic-reduced-model",
        "note": "本扫描只回答'条带取值敏感性'，不改变 BM-01 的 A1 判定与证据等级；参数仍为量级示例。",
        "bands_mm": [round(b * 1000, 1) for b in BANDS_M],
        "scenarios": rows,
    }
    p = os.path.join(RAW, "bm01_band_sensitivity.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    # 打印
    print("BM-01 判定条带灵敏度扫描  （记录条带 = %.0f mm）" % (rows[0]["recorded_band_m"] * 1000))
    print()
    hdr = "场景".ljust(22) + "".join(f"{b*1000:>6.0f}" for b in BANDS_M)
    print(hdr + "   翻转点")
    print("-" * len(hdr) + "---------")
    for e in rows:
        line = e["label"].ljust(20) + "  "
        for s in e["sweep"]:
            mark = "A1" if s["code"] == "A1" else "A3"
            line += f"{mark:>6}"
        flip = e["min_band_mm_for_A1"]
        line += f"   {'无' if flip is None else f'{flip:.0f} mm'}"
        print(line)
    print()
    print("各场景最坏峰值（mm）：")
    for e in rows:
        print(f"  {e['label']:<22} {e['recorded_worst_peak_mm']:>8.2f}   "
              f"（记录判定 {e['recorded_code']}）")
    print()
    print("写出：", p)


if __name__ == "__main__":
    main()
