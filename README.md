# Feedback quantity flips the absorbability verdict for upstream steps in canal pools

Code and data for the manuscript

> Zhang, X. (2026) *Feedback quantity can flip the absorbability verdict for upstream steps in canal pools:
> a controlled numerical counterexample*. Submitted to **Journal of Hydroinformatics** (IWA Publishing).

This repository contains the 1D Saint-Venant solver, the lumped 0D storage model, the controller, the
independent EPA SWMM replication, every raw result file behind the figures and tables, and the plotting
scripts. It is the material referred to in the manuscript's *Data and code availability* statement and in
Supplementary Material note **S11**.

---

## 1. What is in here

```
.
├── bm01/                     synthetic canal pool (the controlled counterexample), L = 2000–12000 m
│   ├── bm01_sve1d_core.py        1D Saint-Venant solver, minimal core (~120 lines)
│   ├── bm01_sve1d_v2.py          main counterexample driver (pool-mean vs gate-adjacent feedback)
│   ├── bm01_sve1d_v2_sensor.py   grid independence of the two feedback channels
│   ├── bm01_pool_step.py         lumped 0D storage model + A1/A3 decision rules (the benchmark case)
│   ├── bm01_ext_swmm.py          EPA SWMM 5.2 replication of the main case (via pyswmm)
│   ├── bm01_ext_compare.py       side-by-side of the in-house and SWMM results
│   ├── raw/                      every archived JSON/CSV result and the SWMM input deck
│   ├── figures/                  archived figure files (PNG/SVG)
│   ├── multiday/                 multi-day campaign with 1-min time series
│   ├── normalized/               normalised record table (records.jsonl)
│   └── review/                   reproducibility note and internal review record
└── bm01r/                    real reach (Tanghe–Fangshuihe, 25.6 km) parameterised from published data
    ├── bm01r_*_check.py          geometry, gate, mass-balance, convergence and 0D/1D checks
    ├── bm01_pred_bm01r.py        the two feedback arms on the real reach
    └── raw/, figures/            archived results
```

## 2. Mapping to the manuscript

| Manuscript item | Script | Archived result |
|---|---|---|
| §2.1 governing equations, numerical scheme | `bm01/bm01_sve1d_core.py` | — (solver) |
| §2.2 lumped 0D model, volume balance | `bm01/bm01_pool_step.py`, `bm01/bm01_p0p_0d_1d_equivalence.py` | `raw/bm01_p0p_0d_1d_equivalence.json` |
| §2.3 controller, feedback quantities, metrics | `bm01/bm01_sve1d_compare.py` (`Ctrl`), `bm01/bm01_sve1d_v2.py` | `raw/bm01_sve1d_v2.json` |
| §2.4 solver verification, grid sensitivity | `bm01/bm01_ver_solver.py`, `bm01/bm01_w8ap_grid_table.py`, `bm01/bm01_g_grid_robustness.py`, `bm01/bm01_window_convergence.py` | `raw/bm01_ver_*.json`, `raw/_gridconv_960.json` |
| Table 2, Figure 2 (the counterexample) | `bm01/bm01_sve1d_v2.py` | `raw/bm01_sve1d_v2.json` |
| §3.5 assessment window, peak capture | `bm01/bm01_p11_sampling_window_scan.py`, `bm01/bm01_window_convergence.py` | `raw/bm01_p11_sampling_window_scan.json` |
| §4 gain-plane scans (tuning independence) | `bm01/bm01_tune_map.py`, `bm01/bm01_tune_valley.py`, `bm01/bm01_gain_sensitivity.py`, `bm01/bm01_e3_gainset_window_s3.py` | `raw/bm01_tune_map.json`, `raw/bm01_gain_sensitivity.json` |
| §5 observability, volume-balance estimate | `bm01/bm01_obs_estimator.py` | `raw/bm01_obs_estimator.json` |
| §6 Mechanism A, Table 3 (control period vs travel time) | `bm01/bm01_pred_ctrlperiod.py`, `bm01/bm01_w3pp_6000.py`, `bm01/bm01_l4000_trigger.py`, `bm01/bm01_w4pp_two_type.py` | `raw/bm01_pred_ctrlperiod.json`, `raw/bm01_w3pp_6000.json` |
| §6.3 the real reach (25.6 km) | `bm01r/bm01_pred_bm01r.py`, `bm01r/bm01r_*.py` | `bm01r/raw/bm01_pred_bm01r.json` |
| §7 Mechanism B, closed form and dispersion chain | `bm01/bm01_w5_kc.py`, `bm01/bm01_w5_kc_static.py`, `bm01/bm01_w5_kc_identity.py`, `bm01/bm01_sve1d_v2_sensor.py` | `raw/bm01_w5_kc*.json`, `raw/bm01_sve1d_v2_sensor.json` |
| Figure S6 mechanism time series | `bm01/bm01_review_mech_timeseries.py`, `bm01/bm01_review_mech_plot.py` | `bm01/figS6_mech_ts.svg` |
| §8 independent-solver verification, Figure 4 | `bm01/bm01_ext_swmm.py`, `bm01/bm01_ext_compare.py`, `bm01/bm01_ext_figure.py` | `raw/_ext_swmm*.json`, `raw/_ext_swmm_pool.inp` |
| Figure 1 (schematic) | `bm01/_make_fig0_schematic.py` | `figures/` |
| Decision-code sensitivity to the ±150 mm band | `bm01/bm01_band_sensitivity.py`, `bm01/bm01_m3_band_sensitivity.py` | `raw/bm01_band_sensitivity.json` |

The remaining `bm01_*.py` files are the exploratory and diagnostic runs that led to these results. They are
kept so that the audit trail is complete; not all of them are cited in the paper.

## 3. Quick start

```bash
python -m pip install -r requirements.txt

# 1D counterexample (main case, ~5 min on a laptop)
cd bm01 && python bm01_sve1d_v2.py

# lumped 0D benchmark case, all four scenarios, with a re-run consistency check
python bm01_pool_step.py --scenario all --out . --verify

# band-sensitivity audit (no re-run of physics: reads archived peaks)
python bm01_band_sensitivity.py

# independent solver (needs EPA SWMM 5.2 through pyswmm)
python bm01_ext_swmm.py && python bm01_ext_compare.py
```

Each driver writes its JSON to `raw/` and prints the decision codes it obtains. See `RUN.md` for the
expected numbers, so you can check that a re-run lands on the published values.

## 4. Expected headline numbers

| Quantity | Value in the paper | Where it comes from |
|---|---|---|
| Pool-mean arm peak, main case | 75.0 mm (A1, absorbable) | `bm01_sve1d_v2.py` |
| Gate-adjacent arm peak, main case | 470.7 mm (A3, non-absorbable) | `bm01_sve1d_v2.py` |
| Pool-mean peak at dx = 50 / 25 / 12.5 m | 74.7 / 75.0 / 74.8 mm | `bm01_w8ap_grid_table.py` |
| Real reach, pool-mean vs gate-adjacent | 49.7 vs 338.8 mm (flip A1 → A3) | `bm01r/bm01_pred_bm01r.py` |

## 5. Scope and caveats

* BM-01 is a **synthetic reduced-model benchmark case**, not a calibrated model of a specific canal. The
  evidence grade is stated in `bm01/README.md`.
* The band ±150 mm is an assumed value for the benchmark, with its own sensitivity study
  (`bm01_band_sensitivity.py`), not a design standard.
* The SWMM arm is a *comparability* design: the orifice coefficient is tuned so that the SWMM gate law
  equals the in-house gate law pointwise (derivation in the header of `bm01_ext_swmm.py`). It is an
  independent-code check, not an independent-physics check.
* Nothing here is an engineering recommendation for a real canal.

## 6. License and citation

* Code: MIT (`LICENSE`). Data and figures: CC BY 4.0 (`LICENSE-DATA`).
* Citation metadata: `CITATION.cff`. If you use this material, please cite the Journal of Hydroinformatics
  article above, and the archived release (DOI) when it becomes available.

## 7. Environment used

Windows 11, Python 3.13, numpy 2.4.6, matplotlib 3.11.2, pyswmm (EPA SWMM 5.2 engine). Exact versions and
file checksums are in `PROVENANCE.md` and `checksums_sha256.txt`.
