# 2D water: FLIP/APIC fluid

Vanilla HTML/CSS/JS water with no build step and no libraries. Open `index.html` in any static server (for example Replit's web Preview). Only water particles are simulated. Empty canvas is air. A surface mesh (metaballs) comes later.

## Solver
`sim.js` is a **FLIP/APIC** hybrid on a **MAC (staggered) grid**:

- Particles carry position, velocity and an APIC affine matrix **C** (2x2).
- P2G uses APIC transfers with quadratic B-spline weights. G2P rebuilds the velocity and C, then blends in a FLIP delta (`flipRatio`, default 0.8).
- Grid steps, in order: mark solid, fluid and air cells; gravity; light viscosity and damping; box walls (**separating** free-slip walls by default, see below); pressure projection by **preconditioned conjugate gradient with a MIC(0) preconditioner** (modified incomplete Cholesky, τ = 0.97). PCG stops on the residual: relative residual ≤ 1e-6, or max residual·dt ≤ 1e-4 (in |div u|·dt units). `params.pcgMaxIters` (default 500) is only a safety cap; a hit is logged to `FlipSolver.pcgCapLog` with a rate-limited console warning. The old red-black Gauss-Seidel (SOR 1.9, fixed `params.pressureIters`) is still available with `FlipSolver.params.pressureSolver = "gs"`. A **drift correction** compares particle density per cell with the rest density, so volume does not creep. Velocities are then extrapolated 2 faces into the air.
- Timestep: each substep uses `dt = CFL·h / (u_max + sqrt(g·h))` with CFL = 1, up to 8 substeps per 1/60 s frame. The raw frame gap is never used. If a frame reaches the cap, it is logged (`FlipSolver.capLog`, plus a rate-limited console warning) and the sim falls behind real time instead of exceeding CFL 1. Particles get a light push-apart pass and are clamped to the walls. The push-apart pass is **momentum-neutral** by default (`params.separation = "neutral"`):
  - Every pair push is equal and opposite.
  - Each axis is cut to what both particles can move without leaving the wall bounds, so the wall clamp can never undo one side of a push.
  - As a result the pass cannot move the centre of mass (measured |Δy_com| ≤ 4e-9 m per substep, float32 rounding).
  - `params.separation = "legacy"` restores the old pass. There the clamp could cancel the half of a push that went into the floor, which lifted the fluid.
  - Known side effect: a pair stacked almost vertically on the floor can no longer be pushed apart vertically, so the floor row ends up slightly crowded (see Validation).
- **Walls** (`params.boundary = "separating"`, default, step 2g):
  - Starts from the 2d free-slip mirror walls (normal = 0 both ways; ghost faces massless).
  - Separating (Signorini) contact via an active set on top of PCG: release any wall face whose fluid cell has p < 0, then put back any released face the fluid moves into. Cap 5 solves per substep. Never clamps pressure to max(p,0).
  - `params.boundary = "mirror"` restores the 2d two-sided walls. `params.boundary = "legacy"` restores the old one-sided walls. Comparison only.
  - Harness: `resetColumn()` reseeds the RNG, and `test/validate.js` runs each column in a fresh process.

- `params.damping` (default 0.15 /s) is a light linear damping that lets the browser pool settle. Set it to 0 for physics validation.
- World: the canvas width is a **3 m tank**, with gravity 9.81 m/s². The grid is about 10 px per cell (64–128 cells wide; 128x72 on a 16:9 desktop) with about 3.2 particles per cell. The default dam break uses about 6,400 particles.

The code has these sections:
0. **`ObstacleShapes`**: pure geometry for fixed obstacles (signed distance functions, shape builders, freehand simplification). No DOM.
1. **`FlipSolver`**: the physics core. It never touches the DOM or canvas. It exposes typed arrays (`px`, `py`, `vx`, `vy`, `count`, `nx`, `ny`, `cellSize`, `radius`) and the obstacle API (see below).
2. **`Canvas2DRenderer`**: reads those arrays and draws discs coloured by speed (blue to white foam), then the obstacles. Its interface is `{ resize(w, h, dpr), draw(solver, view[, overlay]) }`; `overlay` is the optional Object Edit Mode state (ghost shape, selection, hint line). A future WebGL2 renderer only has to implement the same interface.
3. **App**: canvas sizing, HUD, sliders, mouse/touch, Object Edit Mode and the main loop.

## Controls
| control | id | default | what it does |
|---|---|---|---|
| water amount | `sl-n` | 22 % | Fraction of the tank area in the dam column. Changing it restarts the dam break. The HUD shows the particle count. |
| gravity | `sl-g` | 9.81 m/s² | Gravity. |
| FLIP blend | `sl-flip` | 0.80 | 0 = pure APIC (smooth, a little damped). 1 = full FLIP delta (most lively, noisier). |
| pressure solver max iterations (PCG cap) | `sl-iters` | 500 | Safety cap on PCG iterations per substep. PCG normally stops on its residual after about 7–17 iterations, so this rarely changes anything. A very low cap (10) can stop the solve early in violent frames; each cap hit is logged. In GS mode (`params.pressureSolver = "gs"`) the slider sets the fixed Gauss-Seidel iteration count instead. |
| viscosity | `sl-mu` | 0.05 | Grid velocity diffusion. 0 = water, 1 = syrupy. |

Buttons: **Pause**, **Reset**, **Dam break** (a column 0.58 of the tank height on the left, then released), **Edit objects**. **Hold the left mouse button** (or a finger) to hold a round drop of water at the cursor. Release it to throw it with the mouse velocity. Keys: **Space** pauses, **R** resets (obstacles are kept), **E** toggles Object Edit Mode.

### Object Edit Mode (fixed obstacles)
Press **E** or the **Edit objects** button. The sim pauses (your previous pause state comes back when you leave). While editing, the mouse/touch never spawns water. On phones the sliders are hidden while editing so more of the tank can be touched.

| action | mouse / keyboard | touch / panel |
|---|---|---|
| pick shape | **1** Circle · **2** Box · **3** Ramp · **4** Triangle · **5** Freehand | tool buttons |
| place + size | drag on empty canvas (circle: centre → radius; box/triangle: corner → corner; ramp: drag from its high end to its toe). A click without dragging places a default-size shape. | same with a finger |
| freehand | hold and draw a closed loop, release to close it | same |
| select / move | click an object; drag it to move it | same |
| rotate | **[** / **]** (15°, Shift = 3°), mouse wheel, or drag the round handle above the shape (Shift snaps to 15°) | ⟲ / ⟳ buttons or the handle |
| confirm | **Enter** | **Confirm** |
| cancel the ghost | **Esc** (a second Esc clears the selection) | **Cancel** |
| delete selected | **Delete** / **Backspace** | **Delete** (**Clear all** removes every object) |
| leave | **E** | **Done** |

A new or moved shape is a translucent ghost until you confirm it. Confirming writes it into the solver; any particle whose centre is inside it is deleted (the status line shows how many; `stats().obstacleDeleted` counts them since the last reset). Leaving edit mode with a ghost still on screen confirms it; press Esc first to discard it. Edit mode also shows the rasterized solid cells the grid actually sees (faint squares).

**Thickness rule.** On Confirm (Enter, the Confirm button, or Done/E with a ghost on screen) every shape must be at least **2 grid cells (2h, 4.7 cm on the 128-column grid)** thick everywhere. This applies to every tool, including freehand shapes and ramp/triangle tips. The check (`FlipSolver.checkObstacleThickness(shape)`, built on `ObstacleShapes.thicknessCheck`) uses the solver's own raster (a cell is solid if its centre is inside the shape). Tank walls count as solid, so a part lying on the floor or against a wall merges with it. Other obstacles do not count, because they can be removed later. A cell fails when:
- it is more than h/2 from every disc of diameter 2h that fits inside the shape plus the walls. This is a morphological opening, searched on an h/4 lattice. A strip thinner than 2h has no such disc at all; 90° corners pass and corners sharper than about 84° fail. Or:
- it belongs to a one-cell-thick row or column of the raster (two adjacent solid cells with open cells on both sides). That is the configuration behind the one-cell energy blow-up described below.

What happens on Confirm:
- **Polygon passes after clipping** (ramp, triangle, freehand): if the only problem is a sharp corner, the editor clips every convex corner sharper than 88°. The first cut is 2 cells long and grows in half-cell steps up to 6 cells until the shape passes; a cut never removes more than 45 % of an edge. The clipped polygon is placed, keeps its kind, and the status line says "sharp tip(s) clipped so every part is at least 2 grid cells thick". For example, a ramp in mid-air loses its toe point; a ramp whose toe lies on the floor does not.
- **Anything else that fails is not placed** (a thin neck in a freehand shape, or a shape smaller than one cell). The failing cells are drawn in red, the status line says "Not placed: parts of this freehand (red) are thinner than 2 grid cells (4.7 cm). Make it thicker, or Esc to cancel.", and the ghost stays so you can resize or move it. Done/E will not leave edit mode until the shape is fixed or discarded with Esc.

Boxes and circles drawn with the editor are already at least 3 cells across (circle r ≥ 1.5 cells), so they always pass. The API `addObstacle` does not apply the rule: scripted shapes may still be thinner, and the solver handles one-cell-thick solids itself (see below).

### How obstacles work (physics)
- Shapes (solver units, metres, y down): `{ type: "circle", cx, cy, r }`, `{ type: "rect", cx, cy, w, h, angle }` (oriented box), `{ type: "polygon", cx, cy, angle, points: [[x, y], ...], kind }` (ramp, triangle and freehand are polygons; points are local around the centroid). Each has an exact signed distance function; the polygon inside test is even-odd. Builders: `ObstacleShapes.circle / box / ramp / triangle / polygon / freehand`.
- **Grid**: obstacles are rasterized into the solid-cell mask (a cell is solid if its centre is inside a shape) and joined with the tank walls whenever they change. Every grid step that handles wall cells treats obstacle cells the same way: wall-normal faces are zero (or released by the separating active set exactly like box-wall faces), faces inside an obstacle carry no P2G mass and are refilled from the fluid side (normal neighbours count as 0, tangential ones are copied, which reproduces the walls' free-slip mirror on flat sides). A ghost face inside a solid exactly one cell thick has fluid-side faces on both sides, and both sides' stencils read it. Such a face takes the P2G-mass-weighted average of those faces instead of the plain average (see "One-cell-thick obstacles" below). The density sampler, drift term, PCG matrix and extrapolation all see them as SOLID.
- **Particles** collide with the exact shape, like the wall clamp: after advection, a particle closer than its radius is projected out along the SDF gradient and its velocity into the obstacle is removed (tangential velocity kept). A particle wedged between an obstacle and a wall jumps to the nearest free spot (counted in `stats().pushStuck`).
- `reset()` keeps obstacles (pass `reset(mode, { clearObstacles: true })` to drop them); `reset` and `resetColumn` never spawn particles inside or touching an obstacle. A held drop released over an obstacle loses the particles that overlap it (`stats().dropClipped`).
- Freehand strokes are closed, simplified with Ramer-Douglas-Peucker (eps = half a cell, at most 64 vertices). A stroke that crosses itself becomes its convex hull; strokes that are a line or thinner than about 1.5 cells are rejected with a message.
- API: `addObstacle(shape[, { removeParticles }]) → id`, `updateObstacle(id, shape)`, `removeObstacle(id)`, `clearObstacles()`, `obstacles` (read-only copies with ids), `hitTest(x, y[, tol]) → id | -1`, `obstacleSdf(x, y)`, `lastRemovedParticles`, `obstacleMask`. Obstacles are stored in metres, so they survive a window resize.
- With no obstacles the solver is bit-identical to the version before obstacles (checked by `test/obstacles.js`; `test/validate.js` output files are byte-identical apart from wall-clock timings).

**Known limitations**
- The pressure solve sees obstacles as whole grid cells (staircase edges at about 10 px per cell). Particles follow the smooth shape, so a few particles (about 0.1 % in the demo) can sit in a cell whose centre is inside the shape; the grid treats their cell as solid.
- Features thinner than about one cell may not appear on the grid at all; particles still bounce off the exact shape there. The editor's thickness rule keeps such features out of edited shapes (sharp tips are clipped, thin necks rejected). Shapes added through the API are not checked.
- Obstacles are static (no moving or rotating solids while the sim runs) and have no friction setting (free-slip, like the walls). In the comparison-only `pressureMass = "weighted"` mode the mass inside obstacles is dropped instead of folded like the wall ghosts.
- Obstacles make the pressure problem harder: in the 10 s demo dam break the separating active set needs about 2.1 solves per substep instead of 1.6, and PCG needs more iterations, so each substep costs about 40 % more (about 10–20 % more per frame, because the slowed flow takes fewer substeps).

### One-cell-thick obstacles (item 1)
**Bug.** A full-height barrier exactly one grid column thick (b1) blew up during a dam break. With the speed clamp off and damping 0, energy left E0 0.06 s after contact and the run hit NaN at t = 0.236 s (CFL 1) or 0.220 s (CFL 0.5). A 2-column barrier (b2) was fine.

**Cause.** The ghost v faces inside a one-column barrier have a copy source on both sides: the wet face upstream and the dry face downstream. `fillGhosts` averaged the two equally, so the value changed between two points of the same substep:
- At P2G time, which is the FLIP reference (`uOld/vOld`), the dry face is 0, so the ghost = v_wet/2.
- After the two-layer extrapolation, the dry face holds a copy of v_wet, so the ghost = v_wet.

The FLIP delta therefore added about +v/2 at the ghost every substep (×0.8 FLIP blend × stencil weight). That compounded into exponential growth of the tangential velocity along the barrier. The per-stage budget (`test/item1_budget.js` → `validation/item1_stage_budget.json`) shows it. In the growth window [1.01·E0, 10·E0] at CFL 1, extrapolation contributes +128 % of the measured energy change and the obstacle ghost fill +21 %; at CFL 0.5 the figures are +129 % and +20 %.

**Fix.** This is the only solver change. Thin ghosts (level-1 ghosts with ≥ 2 copy sources, flagged in `buildGhosts`) weight their copy sources by P2G face mass (`uW/vW`). A dry side has zero mass and never contributes, so the ghost holds the same value at P2G time and after extrapolation. Every other ghost uses the old code path. Obstacles at least 2 cells thick, and runs with no obstacles, are bit-identical to `backup_pre_item1/sim.js` (checked in `test/obstacles.js`), and the `test/validate.js` tables are byte-identical to the citation dump.

**Gates.** `test/item1_gates.js` → `validation/item1_gates.json` runs b1 and b2 for 30 s at CFL 1, damping 0, clamp off. Both runs:
- have 0 substep-cap hits and 0 PCG cap hits;
- reach max |u| 9.98 m/s;
- end with E(30)/E0 = 0.776 (b1) and 0.773 (b2);
- have 0 leaks and min SDF = 1.0 r;
- push-out PE never fires, and normal-velocity removal is ≤ 0.

The hydrostatic force on the wet face over the last 5 s is 1.045 × ½ρgd² for b1 and 1.040 × for b2, inside the ±2h/d = ±6.3 % band.

Surface tension was removed. In a 3 m tank, real capillary forces are smaller than one grid cell. Rest density and stiffness were SPH concepts and were replaced by FLIP blend and the pressure solve.

Developer options: `index.html?warm=3` runs 3 s of simulation before the first frame (for screenshots). `index.html?demoObstacles=1` preloads a circle and a ramp in the dam-break path (`obstacles_preview.png` is `?demoObstacles=1&warm=0.6`).

## Test
`node test/obstacles.js [seconds width height]` tests the obstacles headless. It checks that the no-obstacle path is bit-identical to `backup_pre_objects/sim.js`, runs a 10 s dam break with a circle and a ramp against the same run without them (leak check every frame: no particle centre inside an obstacle; NaNs; particle count; fluid cells and interior density; ms per frame and per substep), drops a box into the settled pool (particles deleted and counted, no leaks after), and checks the API and freehand builder. Item 1 checks:
- a 2-cell barrier is bit-identical to `backup_pre_item1/sim.js`;
- a 1-cell barrier keeps its energy bounded like the 2-cell one (the pre-fix build reaches E/E0 = 2.86 with the speed clamp pinned at 14 m/s);
- the thickness rule cases: thin box and bar fail; 3-cell boxes at any angle, circles, floor shelves, floor ramps and the demo shapes pass; a mid-air ramp is clipped and passes; a freehand shape with a thin neck fails.

Item 1 harnesses: `test/item1_budget.js` (per-substep stage energy budget at CFL 1 and 0.5, `--sim` to point at another build), `test/item1_gates.js` (30 s frozen gates and hydrostatic force) and `test/item1_regression.js` (citation byte-compare, Run 1 Δd, plain-build b1/b2).

`node --experimental-websocket test/browser_smoke.js` (Node 20; plain `node` on Node 22+) drives `index.html` in headless Chrome over the DevTools protocol with real mouse, keyboard and touch events. It covers edit mode, every tool, move, rotate (keys, wheel, handle), delete, Esc, Enter, the Done/Confirm buttons, drops being disabled while editing, and an iPhone-sized touch viewport. It also covers the thickness rule: a freehand dumbbell with a 1.2-cell neck is rejected on Enter with the "Not placed" message and stays pending, Esc discards it, and a mid-air ramp is placed with its tips clipped. It writes `validation/item1_thin_reject.png` (the red failing cells), `obstacles_preview.png`, `obstacles_editmode.png` and `obstacles_editmode_phone.png`. It needs `google-chrome` or `chromium` and no npm packages.

`node test/headless.js [width height flipRatio viscosity]` runs the solver core in Node with no DOM. It runs a 30 s dam break, then a dropped drop, then 30 thrown drops. It reports NaN and out-of-box counts, volume drift, settling, and ms per frame. It also reports substep-cap frames: how much sim time each one covered, and the maximum CFL in them and over the whole run.

## Validation (Martin & Moyce 1952 column collapse)
`node test/validate.js [--solver pcg|gs] [--sep neutral|legacy] [--boundary separating|mirror|legacy] [--suffix _x]` uses a 128x72 grid (h = 3/128 m) with damping 0 and an exact column from `FlipSolver.resetColumn(a, H)` (2x2 particles per cell, edges on cell edges). It runs H/a = 2 (24x48 cells) and H/a = 1 (24x24 cells) for 30 s each and writes, to `validation/`:
- `front_*`, `energy_*`, `volume_*` and `div_*` CSVs (the div CSV also has PCG iterations and residual per substep);
- `comparison_mm*.csv`, on the release clock and on the shifted clock (T_shifted = T_sim + T_ref − T_1.44, with T_1.44 interpolated);
- `summary*.json` with every metric.

The default suffix is `_<solver>`. `python3 test/plot_validation.py _pcg` (needs matplotlib) makes `front_ZT_pcg.png` (both clocks) and `energy_pcg.png`; `_2c` does the same for step 2c. The reference points in `validation/mm_reference.csv` are the measured "mean" columns of Tables 1 and 2 of the original paper, with sources.

File sets in `validation/`:
- **`*_pcg`**: round 3 (PCG + MIC(0), legacy separation).
- **No suffix** (`front_h2.csv`, `summary.json`, `front_ZT.png`, ...): the round-2 Gauss-Seidel run, kept unchanged.
- **`*_gs`**: the same GS solver re-run with the round-3 harness. It reproduces round 2 and adds the new metrics.
- **`*_2c`**: step 2c (PCG + momentum-neutral separation, legacy walls).
- **`*_2c_legacy`**: the same harness with `--sep legacy`. It reproduces `*_pcg` bit for bit and adds the separation centre-of-mass columns.
- **`*_2d`**: step 2d (mirror walls + neutral separation + PCG).
- **`*_2d_base2c`** (`--boundary legacy`): the 2c state.
- **`*_2d_baseR3`** (`--sep legacy --boundary legacy`): round 3. Both baselines reproduce their earlier runs bit for bit.
- **`*_2g`**: step 2g, the current default (separating walls + neutral separation + PCG).
- **`*_2g_base2d`** / **`*_2g_baseR3`**: 2d and round-3 baselines under the new harness.
- **`*_harness_old`** / **`*_harness_new`**: H/a=1 re-baseline for the RNG reseed (old = same process after h2; new = fresh process + reseed).
- **`energy_budget_2g.json`**, **`test/ceiling_drop.js`**: 2g diagnostics.
- **`partA_post_impact_2g.json`**, **`still_pool_2f.json`**: Part A post-impact budget + 2f still-pool gate (FAIL → stopped).
- **`energy_budget_2d.json`**: stage-by-stage energy budgets over the pre-impact window, from `node test/energy_budget.js [--cells 48|24] [--sep ..] [--boundary ..] [--sim ..]`. The tool writes an instrumented copy of `sim.js` to the OS temp dir.

Round 3 (PCG), against the frozen thresholds:

| metric | H/a = 2 | H/a = 1 | threshold |
|---|---|---|---|
| solver-only max \|div u\|·dt | 9.99e-5 (GS 1.08e-2) | 1.00e-4 (GS 5.5e-3) | < 1e-2, pass |
| PCG iterations mean / p95 / max, cap hits | 9.9 / 12 / 17, 0 | 7.1 / 9 / 13, 0 | |
| interior density / rest, min over samples | 0.937 (212 of 1800 out of band) | 0.937 (294 out of band) | 0.98–1.02, **fail** |
| interior density / rest at 30 s | 0.9985 | 0.9986 | 0.99–1.01, pass |
| fluid cells, max \|N−N0\|/N0 | 17.1 % | 21.5 % | ≤ 2 %, **fail** |
| (N30 − N10)/N0 | −1.22 % | +1.39 % | ≤ 1 %, fail (absolute) |
| max dZ/dT | 1.713 | 2.046 (4 samples > 2) | ≤ 2: pass for H/a = 2, **fail** for H/a = 1 |
| settled depth, last 5 s | 0.2119 m | 0.1061 m | 0.2109 / 0.1055 ± 0.0234, pass |
| Martin & Moyce, release and shifted clock | all pass (shifted minimum 0.928) | all pass (shifted minimum 0.988) | ≥ 0.90 × reference |

Lowering the pressure residual about 100x left density, volume and dZ/dT almost unchanged. Those failures are therefore not caused by the pressure solve. Full numbers are in `validation/summary_pcg.json`.

Step 2c (momentum-neutral separation, `--sep neutral`), against the frozen thresholds. `|N30−N10|` is read as an absolute value.

| metric | H/a = 2 | H/a = 1 | threshold |
|---|---|---|---|
| separation pass, max \|Δy_com\| per substep, after the wall clamp | 2.2e-9 m (legacy 2.0e-5) | 4.0e-9 m (legacy 3.8e-5) | ~0, pass |
| pre-impact energy rise (unsmoothed) | +0.612 % (0.37–0.63 s) | +0.692 % (0.33–0.68 s) | ≤ 0.20 %, **fail** |
| settled depth, last 5 s | 0.2122 m | 0.1052 m | in band, pass |
| interior density / rest, min | 0.954 (149 of 1800 out of band) | 0.932 (281 out of band) | 0.98–1.02, **fail** |
| fluid cells max / \|N30−N10\| | 15.7 % / 1.22 % | 21.4 % / 4.51 % | 2 % / 1 %, **fail** |
| max dZ/dT | 1.696 | 2.041 | ≤ 2: pass for H/a = 2, **fail** for H/a = 1 |
| Martin & Moyce, shifted clock, lowest ratio | 0.927 (a = 1⅛ in.), 0.953 | 0.990, 0.985 | ≥ 0.90, pass |

The neutral pass leaves the pre-impact rise almost unchanged. A per-stage energy breakdown on an instrumented copy shows that, in that window:
- The pass itself adds exactly 0.
- The legacy clamp lift was only 0.0085 J/kg (H/a = 2).
- The surplus comes from two sources:
  - Gravity work against advection: +0.043 J/kg, which halves at CFL 0.5.
  - The `enforceBoundary()` call after projection: +0.029 J/kg, independent of dt.

Side effects of the neutral pass:
- The H/a = 2 pool keeps sloshing (KE 0.021 J/kg at 30 s, versus 0.001 with legacy).
- The floor row holds 12–15 % extra particles, with 170–220 overlapping pairs there. Legacy has none.

Full numbers are in `validation/summary_2c.json`.

Step 2d (mirror walls), against the 2d pass list. H/a = 2 / H/a = 1; the 2c baseline is in brackets.

| item | H/a = 2 | H/a = 1 | result |
|---|---|---|---|
| 1. pre-impact rise (≤ 0.20 %) | 0.586 % (0.612) | 0.662 % (0.692) | **fail** |
| 2. `enforceBoundary` rows over the window (≤ +1e-4 J/kg) | post 0; pre + in-P2G −0.0024 (+0.029) | post 0; −0.0071 (+0.023) | pass |
| 3. post-impact spikes (≤ 0.20 % or halved at CFL 0.5) | 0.221 %; 0.219 % at CFL 0.5 | 0.991 % | **fail** |
| 4. upswings 8–30 s (≤ 0.05 %) | max 0.034 % | max 0.030 % | pass |
| 5. settled depth | 0.2116 m | 0.1056 m | pass |
| 6. Martin & Moyce, shifted clock | 1⅛ in. 0.924, 2¼ in. 0.953 | 0.989, 0.986 | pass |
| density, every sample | min 0.927 (99 of 1800 out of band) | min 0.941 (130 out of band) | **fail** |
| fluid cells max / \|N30−N10\| | 21.6 % / 0.78 % | 19.6 % / 3.65 % | **fail** / pass, **fail** / **fail** |
| dZ/dT | 1.702 | 2.040 | pass / **fail** |

The wall now does no positive work in the budget, but the pre-impact rise barely moves. The energy the old copy showed reappears in the pressure row. What is left is:
- Gravity on the grid against particle advection: +0.041 / +0.072 J/kg. This is the explicit gravity/advection staggering, and it roughly halves at CFL 0.5.
- The pressure projection: −0.008 / −0.047 J/kg at CFL 1.0. At CFL 0.5 it turns positive (+0.019 / +0.001 without the drift term). The projection uses unit face masses, so it does not reduce the particle-mass-weighted energy.

See `validation/energy_budget_2d.json` and `summary_2d.json`.

Step 2g (separating walls), against the 2g pass list. H/a = 2 / H/a = 1.

Harness first: `resetColumn()` reseeds the RNG; each column runs in a fresh process. H/a=1 post-impact spike OLD harness **0.991 %** → NEW **0.824 %** (see `validation/harness_reseed_h1.json`). H/a=2 is bit-identical across harnesses.

Physics: active-set separating contact on top of 2d mirror walls. Release any wall face whose fluid cell has p < 0; put back any released face the fluid moves into (re-project from the current field). Cap 5 solves. Never clamps pressure. `params.boundary = "mirror"` / `"legacy"` keep the old walls.

| item | H/a = 2 | H/a = 1 | result |
|---|---|---|---|
| Wall energy rows (pre-impact) | −0.0024 J/kg | −0.0071 J/kg | PASS |
| Late upswings 8–30 s ≤ 0.05 % | max **0.070 %** at 8.73–8.97 s | max 0.034 % | **FAIL** / PASS |
| Ceiling-drop t_release ≤ legacy+0.10 s and ≤ 0.55 s | center 0.38 (leg 0.32) PASS; left 0.40 (0.06) **FAIL**; right 0.50 (0.08) **FAIL** | — | **FAIL** |
| Dam-break lid t_release | 0.44 s (legacy 1.82) | — | PASS |
| Active-set: worstNegP / worstUInto / capNeg / compFail | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | PASS |
| Active-set passes mean / p95 / max, cap hits | 1.46 / 3 / 5, 0 | 1.14 / 2 / 5, 0 | — |
| Solver-only max \|div u\|·dt | 1.00e-4 | 9.97e-5 | PASS |
| Pre-impact rise ≤ 0.20 % | 0.586 % | 0.662 % | **FAIL** |
| M&M shifted, 1⅛ in. ratio | **0.9235** | — | PASS (2.35 % margin) |

Full numbers in `validation/summary_2g.json`.

### Part A — post-impact budget (2g vs 2d) and ceiling corner band
H/a=1 window **1.5–1.7 s**, fresh process + reseed. Stage rows (J/kg): gravity dominates the positive side; wall `enforceBoundary` totals are **−0.0080 (2g)** / **−0.0081 (2d)** — walls **do not** feed the spike. See `validation/partA_post_impact_2g.json`.

Ceiling dual-face |v_side| under 2g and Validation corner band `t_phys = √(2h/g) + R/|v_side|` (h=0.0234 m, R=0.117 m, floor 0.1 m/s; pass if t_release ≤ t_phys+0.10 and ≤ 0.55):

| drop | |v_side| (m/s) | t_phys | t_release | band |
|------|-----------|--------|-----------|------|
| left | 0.842 | 0.208 | 0.40 | FAIL |
| right | 0.734 | 0.229 | 0.50 | FAIL |
| center | 4.97 | 0.093 | 0.38 | FAIL |

### Step 2f — mass-weighted pressure (**PARKED FAIL**)
Default `params.pressureMass` reverted to **`"equal"`**. `"weighted"` remains behind the flag.

One dig cycle (still-pool after each item; see `validation/still_pool_2f_dig.json`):
1. Shared `uScale/vScale = ρ0/ρ_f` for both Laplacian and `u -= scale·Δp` — still grows (max|u| 0.367).
2. Fold normal ghost face mass onto interface faces — unchanged (0.367).
3. Explicit `ρ_f = m_f/V_f` with `V_f = h²`; clamp frac ~0.021% — unchanged (0.367).

Equal baseline still fails the frozen assert (max|u| ≈ 0.013). Weighted growth is a separate bug. Dam-break for 2f **not run**.

### Step 2e — floor-row sideways separation redirect
When the floor fully blocks a vertical neutral push, redirect into ±x (equal-opposite, side-wall clipped, |δ| capped). Artifacts `*_2e*`.
Pairs &lt;0.8r at t=30: H/a=2 **18** (2g 17), H/a=1 **14** (2g 15). Early t=5 slightly worse. Pre-impact still ~0.57%/0.70%. Wall complementarity fail on H/a=2 (`worstUInto` &gt; 0). M&amp;M 1⅛ in. shifted min ratio **0.925** (margin 7.5%).

### Step 2a — soften drift
`DRIFT_RATE` / `params.driftRate` **4.0 → 1.0**. Artifacts `*_2a*`. Pre-impact 0.55%/0.78%; post spikes still high; late upswings appear on H/a=1; pairs t=30 rose (77/52). Softening alone did not fix energy.

### Step 2b — solid-cell density sampler
Physics default `densitySampler: "renorm"` (skip SOLID corners, renormalize). Always also fills `cellDensityLegacy`. Validation **grades on legacy**; both posted under `volume.densitySamplers`. Artifacts `*_2b*`.
- Renorm ρ30 ≈ **1.001 / 0.998**; legacy grade ρ30 **0.987 / 0.972** (many out-of-band — distribution changed under renorm physics).
- Pairs &lt;0.8r at t=30 improved to **10 / 6**.
- Pre-impact **0.58% / 0.92%**; post spikes still fail; late upswings ≥0.02% present (gate ≤0.05% mostly held).
- M&amp;M 1⅛ in. shifted min ratio **0.926**.

Chain complete (2f parked → 2e → 2a → 2b). Waiting on Lead for verdict; do not start further items.

### Obstacle contact at the source (item 2)
Two solver changes, both only when obstacles exist (the no-obstacle run is byte-identical): `limitObstacleApproach(dt)` runs before advection and strips ONLY the into-solid normal velocity (v·n < 0 within reach, SDF gradient n; never outward or tangential velocity) so a particle stops at r + 4e-5 h; and the pair push of `separateParticles` is clipped along each particle's obstacle normal (one shared equal-and-opposite push, centre of mass unchanged). `pushOutObstacles()` stays as a rare backstop. Harness: `node test/item2_gates.js [--sim path] [--tag T]` then `python3 test/plot_item2.py` (writes `validation/summary_item2.json`, `item2_energy.png`, `item2_window_rows.png`). `test/obstacles.js` check 5a (2-cell barrier bit-identical to pre-item-1) was replaced by a finite / none-inside / energy-within-1 % check because item 2 intentionally changes contact. Backup of the previous build: `backup_pre_item2/`.
