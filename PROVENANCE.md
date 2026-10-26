# Provenance

Recorded so that a reader can reproduce the archived numbers and so that the "random seeds and file
versions are recorded in the repository" statement of Supplementary Material note S11 is literally true.

## 1. Packaging

| Item | Value |
|---|---|
| Package assembled | 2026-10-26 |
| Source of the material | working directories of the BM-01 and BM-01R cases |
| Files | 282 (91 in `bm01/`, 12 in `bm01r/`, the rest in `raw/`, `figures/`, `multiday/`, `normalized/`, `review/`) |
| Total size | ≈21 MB |
| Excluded on purpose | third-party book notes, Word reports, `__pycache__`, exploratory scratch directories |

## 2. Environment actually used

| Component | Version |
|---|---|
| OS | Windows 11 |
| Python | 3.13 |
| numpy | 2.4.6 |
| matplotlib | 3.11.2 |
| pyswmm | 1.0.1 |
| swmm-toolkit (EPA SWMM 5.2 engine) | 0.17.0 |

`numpy` and `matplotlib` are the only hard requirements; `pyswmm` is needed for the independent-solver arm
only (`bm01_ext_swmm.py`, `bm01_ext_compare.py`, `bm01_ext_figure.py`).

## 3. Reproduction check performed at packaging time

`bm01/bm01_sve1d_v2.py` was re-run inside this package directory on 2026-10-26 and compared with the copy
of `raw/bm01_sve1d_v2.json` archived before packaging:

* result: **identical byte for byte, except the `runtime_s` field** (129.56 s archived vs 140.59 s on
  re-run);
* every peak, decision code, grid level and conservation figure is unchanged.

`bm01/bm01_band_sensitivity.py` was re-run and produced a file **byte-identical** to the archived one.

The SWMM arm was re-run at the default resolution (`bm01_ext_swmm.py`, ncell = 20) and completed normally,
reproducing the archived pool-mean peak of 71.55 mm.

## 4. Determinism statement

No random numbers enter the deterministic drivers: the 1D solver, the 0D model, the controller, the grid
and window scans, and the SWMM driver are all seed-free. The drivers that do sample
(`bm01_mc_multicase.py`, `bm01_a_24h_antiwindup.py`) set their seed explicitly at the top of the file.

## 5. Integrity

SHA-256 checksums of every file are listed in `checksums_sha256.txt` (generated at packaging time):

```bash
sha256sum -c checksums_sha256.txt      # Linux / macOS
certutil -hashfile <file> SHA256       # Windows, single file
```
