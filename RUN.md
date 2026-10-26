# Reproduction commands and expected results

All commands assume the repository root as the working directory and the environment from
`requirements.txt` (Python 3.11–3.13). Every driver writes its result to `bm01/raw/` (or `bm01r/raw/`)
and prints the decision codes it obtains, so a re-run can be checked against the numbers below.

Decision codes: **A1** = absorbable within the ±150 mm band, **A3** = non-absorbable / degrade.

---

## 1. The controlled counterexample (Table 2, Figure 2, Section 3)

```bash
cd bm01
python bm01_sve1d_v2.py            # ~2–3 min, 18 equations-of-state cases
```

Expected, at dx = 25 m (ncell = 240), L = 6000 m, ΔQ = +45 m³/s, T_c = 300 s:

| Feedback | Pool-mean peak | Verdict |
|---|---:|---|
| pool-mean level | **75.05 mm** | A1 |
| gate-adjacent single point | **470.78 mm** | A3 |

Grid check across dx = 50 / 25 / 12.5 m (Supplementary Tables S4, S5):

```bash
python bm01_w8ap_grid_table.py     # → 74.7 / 75.0 / 74.8 mm for the pool-mean arm
python bm01_sve1d_v2_sensor.py     # grid independence of the two feedback channels
```

## 2. The lumped 0D benchmark case (Section 2.2)

```bash
python bm01_pool_step.py --scenario all --out . --verify
```

Four scenarios, main case +45 m³/s → A1; two-pool chain with local feedback → A3 (pool 2 leaves the band
for 23 min); the same chain with feedforward → A1; +200 m³/s (beyond capacity) → A3.

## 3. Tuning independence (Section 4)

```bash
python bm01_tune_map.py            # gain-plane scan at fixed 0.44 ratio
python bm01_tune_valley.py         # refinement of the valley floor
python bm01_e3_gainset_window_s3.py
```

## 4. Observability (Section 5)

```bash
python bm01_obs_estimator.py       # ideal pool mean vs volume-balance estimate
```

## 5. Mechanism A: control period vs wave travel time (Section 6, Table 3)

```bash
python bm01_pred_ctrlperiod.py     # ratio 0.22 / 0.44 / 0.66 / 1.32
python bm01_w3pp_6000.py
python bm01_l4000_trigger.py
```

## 6. The real reach, 25.6 km (Section 6.3)

```bash
cd ../bm01r
python bm01_pred_bm01r.py          # 49.7 mm (A1) vs 338.8 mm (A3) → the flip survives on a real geometry
python bm01r_masscheck.py          # mass balance, 1.56% of net inflow / 9.4e-7 of throughput
python bm01r_0d_1d_check.py
```

## 7. Mechanism B: spatial averaging (Section 7)

```bash
cd ../bm01
python bm01_w5_kc_static.py        # direct numerical check of the sinc(L/λ) low-pass
python bm01_w5_kc_identity.py      # identity check
python bm01_p11_sampling_window_scan.py
```

## 8. Independent solver: EPA SWMM 5.2 (Section 8, Figure 4)

```bash
python bm01_ext_swmm.py --mode pool --ncell 240 --routing 2 --horizon 43200 --out raw/bm01_ext_swmm_pool_n240.json
python bm01_ext_swmm.py --mode gate --ncell 240 --routing 2 --horizon 43200 --out raw/bm01_ext_swmm_gate_n240.json
python bm01_ext_compare.py         # table of all 16 runs → raw/bm01_ext_compare.json
```

The archived Figure 4 runs use `--horizon 43200` (12 h) on top of the default 7200 s spin-up, i.e. a
50460 s simulation. Dropping `--horizon` silently falls back to the 14400 s default and truncates the tail
of the run: the peaks are unchanged (they occur early) but `a_final` and the late limit cycle are not,
so always pass the flag when comparing against the archive.

Expected SWMM peaks over ncell = 20/40/60/120/240 (all archived in `raw/`):

| Arm | Pool-mean peak range | Gate-adjacent peak range | Verdict |
|---|---|---:|---|
| SWMM, pool feedback | 67.4 – 72.6 mm | 257 – 1025 mm | A1 at every grid |
| SWMM, gate feedback | 202 – 467 mm | 126 – 438 mm | A3 at every grid |

The verdict pair (A1, A3) is therefore reproduced by a code that shares no line with the in-house solver.
A quick smoke run without a grid choice is simply `python bm01_ext_swmm.py` (ncell = 20, ~1 s).

## 9. Decision-code sensitivity to the band (Supplementary Tables)

```bash
python bm01_band_sensitivity.py    # reads archived peaks, re-derives the codes at ±100…±300 mm
python bm01_m3_band_sensitivity.py
```

## 10. Figures

| Figure | Script |
|---|---|
| Figure 1 (schematic) | `bm01/_make_fig0_schematic.py` |
| Figure 2 (main-case time series) | `bm01/bm01_sve1d_v2.py`, then `bm01_sve1d_v2_plot.py` |
| Figure 3 (observability) | `bm01/bm01_obs_estimator.py` |
| Figure 4 (SWMM replication) | `bm01/bm01_ext_figure.py` (needs `raw/bm01_ext_compare.json`) |
| Figure S6 (mechanism time series) | `bm01/bm01_review_mech_timeseries.py`, `bm01_review_mech_plot.py` |

## 11. Determinism

Re-running `bm01_sve1d_v2.py` on the archived input reproduces `raw/bm01_sve1d_v2.json` **byte for byte
except the `runtime_s` field** (checked 2026-10-26 on the packaging machine). No random seeds are used in
the deterministic drivers; the multi-case and Monte-Carlo style drivers
(`bm01_mc_multicase.py`, `bm01_a_24h_antiwindup.py`) fix their seeds at the top of the file.
