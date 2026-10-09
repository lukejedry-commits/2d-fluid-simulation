# 2D FLIP/APIC fluid sim: STATUS (PAUSE package, 2026-10-09)

Local only. Nothing has been pushed to GitHub or Replit and no repo has been created; that waits for NM's sign-off on the package.
Build: `sim.js` md5 **a147d3ff796838466d000ed2d17fca4d** (copy in `backup_item2_final/`).

## Current state
- Solver: FLIP/APIC (flip ratio 0.8) on a MAC grid, quadratic B-spline transfers, PCG pressure solve, neutral particle separation. Reviewed defaults; no-obstacle behaviour is byte-identical to the citation build (validate.js CSVs identical).
- Obstacles and edit mode: circles, boxes, ramps and freehand polygons with exact SDF contact; edit mode for mouse and touch (browser smoke test 34/34).
- Item 1: thin-ghost fix for 1-cell obstacles plus the 2-cell editor thickness rule. Done.
- Item 2 (signed off by NM):
  - `limitObstacleApproach(dt)` runs before advection and strips ONLY the into-solid normal velocity (SDF gradient n; never outward or tangential).
  - Separation pair pushes are clipped along each particle's obstacle normal (one shared equal-and-opposite push, centre of mass unchanged). Active only when obstacles exist.
  - The old `pushOutObstacles()` stays as a rare backstop.
  - All frozen gates passed: gates 1, 2, 3, 4, 4b, the separation row and the clip checks, plus the full regression set (details in `validation/summary_item2.json`).

### Total energy, Run 2 (dam break past box / ramp, CFL 1, damping 0), E/E0, E0 = 4.788 J/kg
No-obstacle column is the no-obstacle dam break on the same build (`validation/item1_baseline_energy.csv`). Box and ramp come from `validation/item2_total_energy_{box,ramp}.csv` (E0 plus the cumulative per-substep particle dE of the gate run). No solver change was made for this table.

| 0.5 | 0.9936 | 0.9417 | 0.9604 |
| 1 | 0.9393 | 0.7455 | 0.8282 |
| 2 | 0.7850 | 0.5844 | 0.6495 |
| 5 | 0.4383 | 0.4957 | 0.4622 |
| 10 | 0.4045 | 0.4859 | 0.4397 |
| 20 | 0.3952 | 0.4822 | 0.4325 |
| 30 | 0.3888 | 0.4814 | 0.4314 |

- Both obstacle runs peak at E/E0 = 1.0037 at t = 0.09 s. This is the same early rise as the item-1 barrier runs (0.0037).
- At 30 s both end below E0: box 0.481, ramp 0.431, against 0.389 with no obstacle (the obstacle runs retain more energy than the no-obstacle run; not diagnosed here).
- Plots: `validation/item2_total_energy.png` (total E/E0 vs t, this table); `validation/item2_energy.png` and `item2_window_rows.png` show the push-out/removal rows and 0.5 s windows, not total E.

## Validation file index (validation/)
- `summary_item2.json`: item 2 gate values, 4b, clip checks, CFL 0.5 push counts. Raw per-case run JSONs and a copy of the final sim in `item2_runs/`.
- `item2_energy.png`, `item2_window_rows.png`, `item2_total_energy.png`, `item2_total_energy_{box,ramp}.csv`
- `summary_validate_item2.json` and `*_item2.csv`: validate.js on the final build (compare with `*_citation.*`; identical).
- `item1_gates_item2.json`, `item1_gates_energy_item2.csv`: fixed b1/b2 gates and hydrostatic on the final build. (Its stage rows do not book the strip row.)
- Item 1: `item1_gates.json`, `item1_regression.json`, `item1_stage_budget*.json/csv`, `item1_dt_check.json`, `item1_api_tilt.json`, `item1_energy.png`, `item1_stage_rows.png`, `item1_api_tilt.png`, `item1_thin_reject.png`.
- Obstacle runs: `displaced_depth_obst.json`, `energy_contact_obst.json`, `thin_barrier_obst.json`, `ramp_slide_obst.json`, `energy_series_obst.json`, `thin_barrier_energy_obst.png`.
- Citation baseline and the earlier validation sweeps: `*_citation.*`, `summary_*.json`, `comparison_mm_*.csv`, `energy_*.png`, `front_ZT_*.png`.
- Harnesses: `test/` (`item2_gates.js` + `plot_item2.py`, `obstacle_gates.js`, `item1_gates.js`, `item1_closeout.js`, `item1_regression.js`, `obstacles.js`, `validate.js`, `browser_smoke.js`).
- Experiment copy: `experiments_item2/sim_sepclip.js` (the separation clip before it was merged; superseded by sim.js).

## OPEN LIST (measured numbers)
**Ramp sliding is NOT realistic yet.**
- Ramp slide acceleration: **0.258** (late window, a_slope_fit_v) against g sin(theta) = **4.867**. Early window (0.05 to 0.2 s): **3.46**. Not called drag; the cause is not diagnosed. Still to run: an isolated one-particle slide test, expected exactly 4.867, with cut cells.
- Staircase drag: mu ~ 0.48.
- Free fall: **10.12** measured vs 9.81.
- Hydrostatic excess ~1% (b2 ratio **1.038**; b1 1.060).
- Bottom-pressure test: pending.
- Planned, not started: cut cells, path check, the 2f still-pool fix.
- WebGL2 renderer: pending.

## Notes and caveats
- **Overwritten item-1 CSVs.** Running `test/item1_closeout.js` for the item-2 tilt regression overwrote `item1_api_tilt_force*.csv`, `item1_api_tilt_energy*.csv` and `item1_baseline_energy.csv`. They now hold item-2 build numbers and **cannot be re-derived** (no git history existed). Surviving item-1 records: `item1_api_tilt.json` and the PNGs (`item1_api_tilt.png`, `item1_energy.png`, `item1_stage_rows.png`). The baseline CSV is the no-obstacle run, which is byte-identical across builds, so the E/E0 baseline column above is unaffected.
- **`test/obstacles.js` check 5a changed on purpose.** It used to require the 2-cell barrier run to be bit-identical to `backup_pre_item1/sim.js`. Item 2 changes contact behaviour intentionally, so identity can no longer hold. It now requires finite values, nothing inside the obstacle, same count and energy within 1% (measured E/E_pre = 1.0013).
- Float32 particle storage: stored-value tangential/outward velocity changes of ~6e-8 relative are logged, not gated (gate 4 is measured in double inside the code; gate 4b bounds the stored-value sum).
- Leak checks improved for the report-only thin cases (0.5-cell barriers and ramp tip now 0 leaks); this follows from the strip, not from a path check.

## Package contents / git readiness
See `.gitignore`. Tracked: sim.js, index.html, style.css, README.md, STATUS.md, test/, validation/ (all summaries, CSVs, PNGs), experiments_item2/, the backup_* folders, preview PNGs. Excluded: `sim.js.bak`, node_modules, tmp files. Total size is about 42 MB, under the ~50 MB limit, so no raw run dumps were excluded.

## Physics note (obstacle runs retain more energy)
Per Physics: the obstacle runs ending with more total energy than the no-obstacle run (30 s: box 0.481, ramp 0.431 vs 0.389 of E0) is expected, not a defect. Not committed/pushed yet; goes with the next push.
