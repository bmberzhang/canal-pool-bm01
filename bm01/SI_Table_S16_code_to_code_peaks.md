### Supplementary Table S16. Case-by-case code-to-code comparison of pool-mean peaks

BM-01 main case, upstream step +45 m3/s, 12 h window, band +/-150 mm. The two codes use line-by-line equivalent PI laws (Kp = 1.6, Ki = 1.6e-3, Tc = 300 s, deadband 10 mm, rate limit 1/3000 s-1). Only the grid resolutions at which both codes were run are comparable row by row; the remaining rows document the within-code grid study.

| Feedback | dx (m) | SVE-1D ncell | SVE-1D peak (mm) | SWMM ncell | SWMM peak (mm) | Diff. (%, SVE denom.) | Verdict SVE / SWMM |
|---|---|---|---|---|---|---|---|
| pool mean | 300 | - | - | 20 | 71.55 | - | - / A1 |
|  | 150 | - | - | 40 | 72.55 | - | - / A1 |
|  | 100 | 60 | 74.87 | 60 | 72.18 | 3.6 | A1 / A1 |
|  | 50 | 120 | 74.69 | 120 | 72.40 | 3.1 | A1 / A1 |
|  | 25 | 240 | 75.05 | 240 | 67.37 | 10.2 | A1 / A1 |
| gate single point | 300 | - | - | 20 | 202.40 | - | - / A3 |
|  | 150 | - | - | 40 | 396.05 | - | - / A3 |
|  | 100 | 60 | 314.79 | 60 | 326.82 | -3.8 | A3 / A3 |
|  | 50 | 120 | 813.76 | 120 | 467.15 | 42.6 | A3 / A3 |
|  | 25 | 240 | 470.78 | 240 | 354.38 | 24.7 | A3 / A3 |

Notes.

1. Difference is defined with the in-house value as denominator, 100 x (peak_SVE - peak_SWMM) / peak_SVE, as stated in Section 8.2.
2. The three rows where both codes were run at the same dx give 3.1% (dx = 50 m), 3.6% (dx = 100 m) and 10.2% (dx = 25 m).
3. The 10.2% row is not scatter: between dx = 50 and 25 m the SWMM pool-mean peak falls from 72.40 to 67.37 mm whereas the in-house value rises slightly, 74.69 to 75.05 mm. The difference is therefore generated inside SWMM by a resolution-dependent loss (240 conduits in series, an artificial junction loss at every connection, and the orifice carried by a single link with a centroid head), not by a disagreement on the physics.
4. The verdict is identical in every comparable row and at every grid in both codes: pool-mean feedback A1, gate-adjacent single point A3. The 0D model gives 72.74 mm (A1), so three model types span 67-75 mm, all A1.
5. The gate single-point arm is compared only qualitatively: its amplitude is grid-dependent in both codes (158.5% and 130.8% dispersion) and the two codes place the single point differently (network node versus last cell).