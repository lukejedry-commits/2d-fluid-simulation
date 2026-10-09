// Frozen-threshold validation runs for the obstacle build (physics review group).
// Usage: node test/obstacle_gates.js            → runs every case, each in a fresh Node process
//        node test/obstacle_gates.js --case X   → one case, prints JSON (used by the parent)
// Writes validation/*_obst.{json,csv}. sim.js is NOT modified: cases that need CFL 0.5
// or the stage-by-stage energy budget load an instrumented copy written to the OS temp
// dir (same technique as test/energy_budget.js). Measurement-only instrumentation: the
// budget copy is checked to give a bit-identical final state to the plain build.
// Common settings: 128x72 grid, h = 3/128 m, damping 0, reviewed defaults, 30 s sim time.
const fs = require("fs");
const os = require("os");
const path = require("path");
const crypto = require("crypto");
const { spawn } = require("child_process");

const ROOT = path.join(__dirname, "..");
const SIM_PATH = path.join(ROOT, "sim.js");
const OUT = path.join(ROOT, "validation");
const NX = 128, NY = 72, H_CELL = 3 / 128, T_END = 30, G = 9.81;
const args = process.argv.slice(2);
const CASE = args.includes("--case") ? args[args.indexOf("--case") + 1] : null;

// ---------------------------------------------------------------------------------
// solver loading (plain, CFL-patched and/or energy-budget instrumented temp copies)
// ---------------------------------------------------------------------------------
function loadSim(opts) {
  const cfl = opts.cfl || 1, budget = !!opts.budget;
  if (cfl === 1 && !budget) return require(SIM_PATH);
  let src = fs.readFileSync(SIM_PATH, "utf8");
  const patch = (a, b) => {
    if (src.split(a).length !== 2) throw new Error("obstacle_gates: patch anchor not found exactly once: " + a.slice(0, 70));
    src = src.replace(a, b);
  };
  if (cfl !== 1) patch("  const CFL = 1.0;", "  const CFL = " + cfl + ";");
  if (budget) {
    patch("    if (mirror) enforceBoundary();\n", "    if (mirror) { const k0 = kg(); enforceBoundary(); bP2G = kg() - k0; }\n");
    patch("      if (d >= pr) continue;\n", "      if (d >= pr) continue;\n      PB.pushed++; PB.depthSum += pr - d; if (pr - d > PB.depthMax) PB.depthMax = pr - d;\n");
    patch("        if (vn < 0) { vx -= vn * o.nx; vy -= vn * o.ny; }\n",
      "        if (vn < 0) { const e0 = vx * vx + vy * vy; vx -= vn * o.nx; vy -= vn * o.ny; PB.velRem += 0.5 * (vx * vx + vy * vy - e0); PB.velRemEvents++; }\n");
    const sub = src.match(/  function substep\(dt\) \{\n[\s\S]*?\n  \}\n/);
    if (!sub) throw new Error("substep() not found");
    const expected = [
      "pinHeld();", "advectParticles(dt);", "separateParticles(SEPARATION_ITERS);", "collideWalls();", "if (obsActive) pushOutObstacles();",
      "particlesToGrid();", "markCells();", "applyGravity(dt);", "applyViscosity();", "enforceBoundary();", "solvePressure(dt);",
      "extrapolateVelocities();", "enforceBoundary();", "if (substepHook) substepHook(divergenceStats(dt), dt);", "gridToParticles();",
    ];
    const body = sub[0].split("\n").slice(1, -2).map((l) => l.trim());
    if (JSON.stringify(body) !== JSON.stringify(expected)) throw new Error("substep() changed; update the instrumented copy:\n" + body.join("\n"));
    patch(sub[0], `  const B = {}; let bP2G = 0; const SER = [];
  const PB = { pushed: 0, depthSum: 0, depthMax: 0, velRem: 0, velRemEvents: 0, substeps: 0, pushedSubsteps: 0 };
  const add = (k, x) => { B[k] = (B[k] || 0) + x; };
  function kp() { let k = 0, q = 0; const fy = (ny - 1) * h, g = params.gravity; for (let i = 0; i < count; i++) { k += 0.5 * (pvx[i] * pvx[i] + pvy[i] * pvy[i]); q += g * (fy - py[i]); } return [k / count, q / count]; }
  function kg() { let k = 0; for (let i = 0; i < u.length; i++) k += uW[i] * u[i] * u[i]; for (let i = 0; i < v.length; i++) k += vW[i] * v[i] * v[i]; return (0.5 * k) / count; }
  function substep(dt) {
    pinHeld();
    let a = kp();
    advectParticles(dt); let b = kp(); add("advection (particle PE)", b[0] + b[1] - a[0] - a[1]);
    separateParticles(SEPARATION_ITERS); a = kp(); add("separation (particles)", a[0] + a[1] - b[0] - b[1]);
    collideWalls(); b = kp(); add("wall clamp (particles)", b[0] + b[1] - a[0] - a[1]);
    const vr0 = PB.velRem, np0 = PB.pushed;
    if (obsActive) pushOutObstacles();
    a = kp();
    const vr = (PB.velRem - vr0) / count, pPE = a[1] - b[1], pKE = a[0] - b[0] - vr;
    add("push-out PE", pPE); add("push-out KE (excl. normal-velocity removal)", pKE); add("normal-velocity removal", vr);
    SER.push([simTime + dt, pPE, pKE, vr, PB.pushed - np0]);
    PB.substeps++; if (PB.pushed > np0) PB.pushedSubsteps++;
    const k0 = a[0]; bP2G = 0;
    particlesToGrid(); markCells(); const g0 = kg(); add("transfer P2G (incl. ghost fold / obstacle ghost mass)", g0 - k0 - bP2G); add("enforceBoundary (inside P2G)", bP2G);
    applyGravity(dt); let g = kg(); add("gravity (grid; damping = 0)", g - g0);
    applyViscosity(); let g2 = kg(); add("viscosity", g2 - g);
    enforceBoundary(); g = kg(); add("enforceBoundary (pre-projection)", g - g2);
    solvePressure(dt); g2 = kg(); add("pressure projection (incl. drift, separating active set)", g2 - g);
    extrapolateVelocities(); g = kg(); add("extrapolation", g - g2);
    enforceBoundary(); g2 = kg(); add("enforceBoundary (post-projection, obstacle ghost fill)", g2 - g);
    if (substepHook) substepHook(divergenceStats(dt), dt);
    gridToParticles(); a = kp(); add("transfer G2P (incl. FLIP, speed clamp)", a[0] - g2);
  }
`);
    patch("    setSubstepHook(fn) { substepHook = fn || null; },", "    setSubstepHook(fn) { substepHook = fn || null; },\n    budget: B, pushStats: PB, pushSeries: SER, energyNow: kp,");
  }
  const tmp = path.join(os.tmpdir(), `sim_obst_gates_${process.pid}_${cfl}_${budget ? "b" : "p"}.js`);
  fs.writeFileSync(tmp, src);
  const mod = require(tmp);
  fs.unlinkSync(tmp);
  return mod;
}

function setup(mod) {
  const S = mod.FlipSolver;
  S.init(NX, NY, H_CELL);
  S.params.damping = 0; // common setting; everything else = reviewed defaults
  return S;
}
function energy(S) {
  const fy = (S.ny - 1) * S.cellSize;
  let k = 0, q = 0;
  for (let i = 0; i < S.count; i++) { k += 0.5 * (S.vx[i] * S.vx[i] + S.vy[i] * S.vy[i]); q += G * (fy - S.py[i]); }
  return { ke: k / S.count, pe: q / S.count, e: (k + q) / S.count };
}
function stateHash(S) {
  const hsh = crypto.createHash("sha256");
  for (const A of [S.px, S.py, S.vx, S.vy]) hsh.update(Buffer.from(A.buffer, A.byteOffset, S.count * 4));
  return hsh.digest("hex").slice(0, 16);
}
function config(S, mod) {
  const p = S.params;
  return { cellSize: S.cellSize, grid: [S.nx - 2, S.ny - 2], radius: S.radius, cfl: mod.FlipSolver.CFL, maxSubsteps: S.MAX_SUBSTEPS, damping: p.damping, flipRatio: p.flipRatio, viscosity: p.viscosity, solver: p.pressureSolver, boundary: p.boundary, separation: p.separation, pressureMass: p.pressureMass, driftRate: p.driftRate, densitySampler: p.densitySampler, gravity: p.gravity, fill: p.fill };
}

// geometry helpers (cells → metres; interior cells are 1..128 x 1..72, floor face at y = 73 h)
const H = H_CELL, FLOOR = (NY + 1) * H;
const OS = () => require(SIM_PATH).ObstacleShapes;
function floorBox(i0, cellsW, cellsH) { // cell-aligned box on the floor: cells i0..i0+cellsW-1, top cellsH rows
  return { type: "rect", cx: (i0 + cellsW / 2) * H, cy: FLOOR - (cellsH / 2) * H, w: cellsW * H, h: cellsH * H, angle: 0 };
}
const BOX1 = { i0: 70, w: 42, h: 8 }; // cells 70..111 on the floor (dam column ends at ~cell 49)
const RAMP2 = () => OS().ramp(112 * H, FLOOR - 8 * H, 70 * H, FLOOR); // same footprint, rises to the right

// ---------------------------------------------------------------------------------
// Run 1 — displaced depth
// ---------------------------------------------------------------------------------
function run1(withBox) {
  const mod = loadSim({});
  const S = setup(mod);
  let removedAtAdd = 0, boxCells = 0;
  if (withBox) {
    S.addObstacle(floorBox(BOX1.i0, BOX1.w, BOX1.h));
    removedAtAdd = S.lastRemovedParticles;
  }
  S.reset("dam"); // default dam-break column (fill 0.22), seeded AFTER the box exists
  const n0 = S.count;
  const st0 = S.stats();
  boxCells = st0.obstacleCells;
  const M = S.obstacleMask, nx = S.nx, ny = S.ny, h = S.cellSize;
  const boxCols = new Set();
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) if (M[i + j * nx]) boxCols.add(i);
  // the same exclusion set (the Run-1 box columns) is used for BOTH runs
  const excl = new Set(); for (let i = BOX1.i0; i < BOX1.i0 + BOX1.w; i++) excl.add(i);
  const cols = []; for (let i = 1; i < nx - 1; i++) if (!excl.has(i)) cols.push(i);
  function depth() {
    // validate.js columnDepth(): floor-connected FLUID cells walked up from the floor;
    // highest particle centre in the top cell + 0.25 h. Averaged over `cols` only.
    const T = S.cellType, topY = new Float64Array(nx).fill(Infinity), topRow = new Int32Array(nx).fill(-1);
    const fluidPerCol = new Int32Array(nx);
    for (let i = 1; i < nx - 1; i++) {
      for (let j = 1; j < ny - 1; j++) if (T[i + j * nx] === 1) fluidPerCol[i]++;
      let j = ny - 2;
      if (T[i + j * nx] !== 1) continue;
      while (j > 1 && T[i + (j - 1) * nx] === 1) j--;
      topRow[i] = j;
    }
    for (let k = 0; k < S.count; k++) {
      const i = Math.floor(S.px[k] / h), j = Math.floor(S.py[k] / h);
      if (topRow[i] === j && S.py[k] < topY[i]) topY[i] = S.py[k];
    }
    let s = 0, sf = 0;
    for (const i of cols) { s += topY[i] === Infinity ? 0 : FLOOR - topY[i] + 0.25 * h; sf += fluidPerCol[i]; }
    return { d: s / cols.length, dCells: (sf / cols.length) * h };
  }
  let nan = 0, inside = 0, minSdf = Infinity, countChanged = false;
  const rows = [], win = [];
  while (S.simTime < T_END - 1e-9) {
    S.step();
    const st = S.stats();
    nan = Math.max(nan, st.nan);
    inside = Math.max(inside, st.insideObstacle);
    if (st.minObstacleSdf < minSdf) minSdf = st.minObstacleSdf;
    if (S.count !== n0) countChanged = true;
    const dd = depth();
    rows.push([S.simTime.toFixed(5), dd.d.toFixed(6), dd.dCells.toFixed(6), st.fluidCells]);
    if (S.simTime >= T_END - 5 - 1e-9) win.push({ d: dd.d, dc: dd.dCells, fc: st.fluidCells });
  }
  const mean = (k) => win.reduce((a, b) => a + b[k], 0) / win.length;
  const st = S.stats();
  return {
    case: withBox ? "box" : "none", config: config(S, mod), particles0: n0, particles30: S.count, countChanged, removedAtAdd,
    obstacleDeleted: st.obstacleDeleted, nanMax: nan, insideMax: inside, minSdfOverR: withBox ? minSdf / S.radius : null,
    boxCells, boxColumnsInMask: [...boxCols].sort((a, b) => a - b), excludedColumns: [BOX1.i0, BOX1.i0 + BOX1.w - 1], columnsAveraged: cols.length,
    d_last5s: mean("d"), d_last5s_min: Math.min(...win.map((w) => w.d)), d_last5s_max: Math.max(...win.map((w) => w.d)),
    dFluidCells_last5s: mean("dc"), fluidCells_last5s: mean("fc"), fluidCells30: st.fluidCells, fluidCells0: rows.length ? +rows[0][3] : null,
    samples: win.length, series: rows,
  };
}

// ---------------------------------------------------------------------------------
// Run 2 — energy budget at contact
// ---------------------------------------------------------------------------------
function run2(kind, cfl, plain) {
  const mod = loadSim({ cfl, budget: !plain });
  const S = setup(mod);
  if (kind === "box") S.addObstacle(floorBox(BOX1.i0, BOX1.w, BOX1.h));
  if (kind === "ramp") S.addObstacle(RAMP2());
  S.reset("dam");
  const n0 = S.count, e0 = energy(S);
  let nan = 0, inside = 0, minSdf = Infinity;
  while (S.simTime < T_END - 1e-9) {
    S.step();
    const st = S.stats();
    nan = Math.max(nan, st.nan); inside = Math.max(inside, st.insideObstacle);
    if (st.minObstacleSdf < minSdf) minSdf = st.minObstacleSdf;
  }
  const e1 = energy(S), st = S.stats();
  const res = { case: kind, cfl, config: config(S, mod), particles: [n0, S.count], nanMax: nan, insideMax: inside, minSdfOverR: kind === "none" ? null : minSdf / S.radius, capHits: S.capLog.length, finalStateHash: stateHash(S), simTime: S.simTime, pushStuck: st.pushStuck };
  if (plain) return res;
  const rows = {}; let sum = 0;
  for (const [k, v] of Object.entries(S.budget)) { rows[k] = v; sum += v; }
  const SER = S.pushSeries, PB = S.pushStats;
  // gate quantities
  const pushPE = rows["push-out PE"] || 0, pushKE = rows["push-out KE (excl. normal-velocity removal)"] || 0, velRem = rows["normal-velocity removal"] || 0;
  let maxWin = -Infinity, maxWinT = null, minWin = Infinity, i0 = 0, acc = 0, maxVelRemSub = -Infinity;
  const pref = [0];
  for (const r of SER) pref.push(pref[pref.length - 1] + r[1] + r[2] + r[3]);
  for (let j = 0; j < SER.length; j++) {
    while (SER[i0][0] <= SER[j][0] - 0.5) i0++;
    const w = pref[j + 1] - pref[i0];
    if (w > maxWin) { maxWin = w; maxWinT = [SER[i0][0], SER[j][0]]; }
    if (w < minWin) minWin = w;
    if (SER[j][3] > maxVelRemSub) maxVelRemSub = SER[j][3];
    acc += SER[j][4];
  }
  Object.assign(res, {
    E0: e0.e, E30: e1.e, dE_measured: e1.e - e0.e, sumOfRows: sum, closure: sum - (e1.e - e0.e), rows,
    pushPE, pushKE, velRem, pushTotal: pushPE + pushKE + velRem, maxVelRemPerSubstep: maxVelRemSub,
    maxWindow05s: maxWin, maxWindow05sAt: maxWinT, minWindow05s: minWin,
    substeps: PB.substeps, substepsWithPush: PB.pushedSubsteps, pushedTotal: PB.pushed,
    pushedPerSubstep: PB.pushed / PB.substeps, pushedPerPushSubstep: PB.pushedSubsteps ? PB.pushed / PB.pushedSubsteps : 0,
    meanPushDepth_m: PB.pushed ? PB.depthSum / PB.pushed : 0, meanPushDepth_r: PB.pushed ? PB.depthSum / PB.pushed / S.radius : 0, maxPushDepth_r: PB.depthMax / S.radius,
    velRemEvents: PB.velRemEvents,
    series: SER.map((r) => [r[0].toFixed(6), r[1].toExponential(6), r[2].toExponential(6), r[3].toExponential(6), r[4]]),
  });
  return res;
}

// ---------------------------------------------------------------------------------
// Run 3 — thin-feature leaks
// ---------------------------------------------------------------------------------
function barrierShape(x0c, x1c) { // full height, extended through lid and floor wall cells
  return { type: "rect", cx: 0.5 * (x0c + x1c) * H, cy: 0.5 * (NY + 2) * H, w: (x1c - x0c) * H, h: (NY + 2) * H, angle: 0 };
}
const BARRIERS = {
  b1: [64, 65], // 1 cell, cell-aligned
  b2: [64, 66], // 2 cells, cell-aligned
  b05c: [64.25, 64.75], // 0.5 cell centred on a cell centre (marked)
  b05e: [64.6, 65.1], // 0.5 cell containing no cell centre (not marked)
};
function run3(kind, cfl) {
  const mod = loadSim({ cfl });
  const S = setup(mod);
  const res = { case: kind, cfl };
  const r = S.radius;
  if (kind === "tip") {
    // elevated thin ramp ("ski jump"), tip pointing upstream: underside horizontal 12 cells
    // above the floor, tip at x = 70 h (zero thickness), back at x = 100 h (4 cells thick)
    const tipX = 70 * H, backX = 100 * H, yb = FLOOR - 12 * H, yt = FLOOR - 16 * H;
    S.addObstacle(OS().ramp(backX, yt, tipX, yb));
    S.reset("dam");
    const n = S.count, side = new Int8Array(n); // 1 above top face, -1 below underside, 0 elsewhere
    const yTop = (x) => yb + ((yt - yb) * (x - tipX)) / (backX - tipX);
    const classify = (x, y) => (x < tipX || x > backX ? 0 : y < yTop(x) ? 1 : y > yb ? -1 : 0);
    for (let i = 0; i < n; i++) side[i] = classify(S.px[i], S.py[i]);
    let t2u = 0, u2t = 0, minSdf = Infinity, insideMax = 0, nan = 0;
    const crossT = [];
    S.setSubstepHook((d) => {
      let ins = 0;
      for (let i = 0; i < S.count; i++) {
        const c = classify(S.px[i], S.py[i]);
        if (side[i] === 1 && c === -1) { t2u++; if (crossT.length < 50) crossT.push(+d.t.toFixed(4)); }
        if (side[i] === -1 && c === 1) u2t++;
        side[i] = c;
        const sd = S.obstacleSdf(S.px[i], S.py[i]);
        if (sd < minSdf) minSdf = sd;
        if (sd < 0) ins++;
      }
      insideMax = Math.max(insideMax, ins);
    });
    while (S.simTime < T_END - 1e-9) { S.step(); nan = Math.max(nan, S.stats().nan); }
    const M = S.obstacleMask; let cells = 0, firstCol = Infinity;
    for (let j = 0; j < S.ny; j++) for (let i = 0; i < S.nx; i++) if (M[i + j * S.nx]) { cells++; firstCol = Math.min(firstCol, i); }
    return Object.assign(res, { config: config(S, mod), particles: [n, S.count], topToUnderside: t2u, undersideToTop: u2t, firstCrossTimes: crossT, minSdfOverR: minSdf / r, insideMaxPerSubstep: insideMax, nanMax: nan, maskCells: cells, maskFirstColumn: firstCol, tipColumn: 70, missingTipLength_cells: firstCol - 70, capHits: S.capLog.length });
  }
  const [a, b] = kind === "none" ? [64, 64] : BARRIERS[kind];
  if (kind !== "none") S.addObstacle(barrierShape(a, b));
  S.reset("dam");
  const xR = (kind === "none" ? 64 : b) * H; // downstream side = right of the barrier's right face
  const n = S.count, ever = new Uint8Array(n);
  let maxDown = 0, maxDownT = null, firstLeakT = null, minSdf = Infinity, minSdfT = null, insideMax = 0, subs = 0, nan = 0;
  let initDown = 0; for (let i = 0; i < n; i++) if (S.px[i] > xR) initDown++;
  S.setSubstepHook((d) => {
    subs++;
    let down = 0, ins = 0;
    for (let i = 0; i < S.count; i++) {
      if (S.px[i] > xR) { down++; ever[i] = 1; }
      if (kind !== "none") {
        const sd = S.obstacleSdf(S.px[i], S.py[i]);
        if (sd < minSdf) { minSdf = sd; minSdfT = d.t; }
        if (sd < 0) ins++;
      }
    }
    if (down > maxDown) { maxDown = down; maxDownT = d.t; }
    if (down > 0 && firstLeakT === null) firstLeakT = d.t;
    insideMax = Math.max(insideMax, ins);
  });
  while (S.simTime < T_END - 1e-9) { S.step(); nan = Math.max(nan, S.stats().nan); }
  let finalDown = 0; for (let i = 0; i < S.count; i++) if (S.px[i] > xR) finalDown++;
  let everN = 0; for (let i = 0; i < n; i++) everN += ever[i];
  const M = S.obstacleMask; let cells = 0; const colsMarked = new Set();
  for (let j = 0; j < S.ny; j++) for (let i = 0; i < S.nx; i++) if (M[i + j * S.nx]) { cells++; colsMarked.add(i); }
  return Object.assign(res, {
    config: config(S, mod), barrier_cells: kind === "none" ? null : [a, b], thickness_cells: kind === "none" ? 0 : b - a, particles: [n, S.count],
    initialDownstream: initDown, maxDownstream: maxDown, maxDownstreamT: maxDownT, finalDownstream: finalDown, everDownstream: everN, firstLeakT,
    minSdfOverR: kind === "none" ? null : minSdf / r, minSdfT, insideMaxPerSubstep: insideMax, substeps: subs, nanMax: nan, capHits: S.capLog.length,
    maskCells: cells, maskColumns: [...colsMarked].sort((x, y) => x - y), maskRowsExpected: NY,
  });
}

// ---------------------------------------------------------------------------------
// Run 4 — sheet sliding on a frictionless ramp (and free-fall baseline)
// ---------------------------------------------------------------------------------
function run4(withRamp) {
  const mod = loadSim({});
  const S = setup(mod);
  const A = [1 * H, FLOOR - 40 * H], C = [71 * H, FLOOR]; // ramp top vertex at the left wall, toe at x = 71 h
  const L = Math.hypot(C[0] - A[0], C[1] - A[1]), th = Math.atan2(40, 70);
  const t = [(C[0] - A[0]) / L, (C[1] - A[1]) / L], nrm = [Math.sin(th), -Math.cos(th)];
  S.resetColumn(5, 5); // 100 particles (no obstacle yet), then moved onto the sheet below
  const r = S.radius, rows = 4, per = 25, sp = 0.6 * H;
  if (S.count !== rows * per) throw new Error("unexpected seed count " + S.count);
  let k = 0;
  for (let rr = 0; rr < rows; rr++) for (let m = 0; m < per; m++, k++) {
    const s = 8 * H + m * sp, dn = 1.05 * r + rr * sp;
    S.px[k] = A[0] + s * t[0] + dn * nrm[0]; S.py[k] = A[1] + s * t[1] + dn * nrm[1]; S.vx[k] = 0; S.vy[k] = 0;
  }
  let removed = 0;
  if (withRamp) { S.addObstacle(OS().ramp(A[0], A[1], C[0], C[1])); removed = S.lastRemovedParticles; }
  const n = S.count;
  const samples = [];
  let minSdf = Infinity, nan = 0;
  const rec = (tt) => {
    let ms = 0, vs = 0, vy = 0, maxS = -Infinity, maxY = -Infinity, dn = 0;
    const vrow = [0, 0, 0, 0];
    for (let i = 0; i < n; i++) {
      const dx = S.px[i] - A[0], dy = S.py[i] - A[1], s = dx * t[0] + dy * t[1];
      const vsi = S.vx[i] * t[0] + S.vy[i] * t[1];
      ms += s; vs += vsi; vy += S.vy[i]; dn += dx * nrm[0] + dy * nrm[1]; vrow[(i / per) | 0] += vsi / per;
      if (s > maxS) maxS = s; if (S.py[i] > maxY) maxY = S.py[i];
      if (withRamp) { const sd = S.obstacleSdf(S.px[i], S.py[i]); if (sd < minSdf) minSdf = sd; }
    }
    samples.push({ t: tt, s: ms / n, vs: vs / n, vy: vy / n, maxS, maxY, dn: dn / n, v0: vrow[0], v1: vrow[1], v2: vrow[2], v3: vrow[3] });
  };
  rec(0);
  S.setSubstepHook((d) => rec(d.t));
  while (S.simTime < T_END - 1e-9) { S.step(); nan = Math.max(nan, S.stats().nan); }
  // clean window: 0.05 s → leading particle 6 cells before the toe (ramp) / lowest particle 5 cells above the floor (free fall)
  const end = samples.find((q) => q.t > 0.05 && (withRamp ? q.maxS > L - 6 * H : q.maxY > FLOOR - 5 * H));
  const t1 = end ? end.t : 1;
  let W = samples.filter((q) => q.t >= 0.2 && q.t <= t1);
  const fitIn = (key, ta, tb) => { const keep = W; W = samples.filter((q) => q.t >= ta && q.t <= tb); const r = fit(key); W = keep; return r; };
  W = samples.filter((q) => q.t >= 0.05 && q.t <= t1);
  const fit = (key) => { // least-squares slope of key vs t
    const m = W.length, mt = W.reduce((a, q) => a + q.t, 0) / m, mv = W.reduce((a, q) => a + q[key], 0) / m;
    let num = 0, den = 0; for (const q of W) { num += (q.t - mt) * (q[key] - mv); den += (q.t - mt) ** 2; }
    return num / den;
  };
  const quad = () => { // s = s0 + v0 t + a t^2 / 2 (least squares)
    const X = W.map((q) => [1, q.t, 0.5 * q.t * q.t]), y = W.map((q) => q.s);
    const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], bb = [0, 0, 0];
    for (let i = 0; i < X.length; i++) for (let a = 0; a < 3; a++) { bb[a] += X[i][a] * y[i]; for (let c = 0; c < 3; c++) M[a][c] += X[i][a] * X[i][c]; }
    for (let a = 0; a < 3; a++) for (let c = a + 1; c < 3; c++) { const f = M[c][a] / M[a][a]; for (let e = a; e < 3; e++) M[c][e] -= f * M[a][e]; bb[c] -= f * bb[a]; }
    const x = [0, 0, 0]; for (let a = 2; a >= 0; a--) { let s = bb[a]; for (let c = a + 1; c < 3; c++) s -= M[a][c] * x[c]; x[a] = s / M[a][a]; }
    return x[2];
  };
  return {
    case: withRamp ? "ramp" : "freefall", config: config(S, mod), theta_deg: (th * 180) / Math.PI, gSinTheta: G * Math.sin(th), particles: [n, S.count], removedAtAdd: removed,
    window: [0.05, t1], windowSamples: W.length, a_slope_fit_v: fit("vs"), a_slope_fit_s: quad(), a_vertical_fit: fit("vy"),
    a_slope_early_0p05_0p2: fitIn("vs", 0.05, 0.2), a_vertical_early_0p05_0p2: fitIn("vy", 0.05, 0.2),
    a_slope_byRow_window: [fit("v0"), fit("v1"), fit("v2"), fit("v3")], // row 0 touches the ramp, row 3 outermost
    a_slope_byRow_early: [fitIn("v0", 0.05, 0.2), fitIn("v1", 0.05, 0.2), fitIn("v2", 0.05, 0.2), fitIn("v3", 0.05, 0.2)],
    v_slope_max: Math.max(...samples.filter((q) => q.t <= t1).map((q) => q.vs)), v_slope_at_window_end: W.length ? W[W.length - 1].vs : null,
    v_slope_byRow_at_window_end: W.length ? [W[W.length - 1].v0, W[W.length - 1].v1, W[W.length - 1].v2, W[W.length - 1].v3] : null,
    meanNormalOffset_r_start: W.length ? W[0].dn / r : null, meanNormalOffset_r_end: W.length ? W[W.length - 1].dn / r : null,
    minSdfOverR: withRamp ? minSdf / r : null, nanMax: nan,
    series: samples.filter((q, i) => i % 2 === 0 && q.t <= Math.max(1.5, t1 + 0.2)).map((q) => [q.t.toFixed(5), q.s.toFixed(6), q.vs.toFixed(6), q.vy.toFixed(6), q.dn.toFixed(6)]),
  };
}

// ---------------------------------------------------------------------------------
const CASES = {
  r1_none: () => run1(false), r1_box: () => run1(true),
  r2_none: () => run2("none", 1), r2_box: () => run2("box", 1), r2_ramp: () => run2("ramp", 1),
  r2_box_plain: () => run2("box", 1, true), r2_ramp_plain: () => run2("ramp", 1, true),
  r2_none_cfl05: () => run2("none", 0.5), r2_box_cfl05: () => run2("box", 0.5), r2_ramp_cfl05: () => run2("ramp", 0.5),
  r3_none: () => run3("none", 1), r3_b05c: () => run3("b05c", 1), r3_b05e: () => run3("b05e", 1), r3_b1: () => run3("b1", 1), r3_b2: () => run3("b2", 1),
  r3_b1_cfl05: () => run3("b1", 0.5), r3_tip: () => run3("tip", 1),
  r4_ramp: () => run4(true), r4_freefall: () => run4(false),
};

if (CASE && !CASES[CASE]) { console.error("unknown case " + CASE); process.exit(2); }
if (CASE) {
  const t0 = Date.now();
  const res = CASES[CASE]();
  res.wallSeconds = (Date.now() - t0) / 1000;
  process.stdout.write(JSON.stringify(res), () => process.exit(0)); // exit only after the pipe drained
}

// ---- parent: one fresh process per case, 6 at a time --------------------------------
const CACHE = path.join(os.tmpdir(), "obst_gates_cache"); // per-case results, so an interrupted suite can resume
function child(name) {
  const cf = path.join(CACHE, name + ".json");
  if (args.includes("--resume") && fs.existsSync(cf)) return Promise.resolve(JSON.parse(fs.readFileSync(cf, "utf8")));
  return new Promise((resolve, reject) => {
    const p = spawn(process.execPath, [__filename, "--case", name], { stdio: ["ignore", "pipe", "pipe"] });
    let out = "", err = "";
    p.stdout.on("data", (d) => (out += d));
    p.stderr.on("data", (d) => (err += d));
    p.on("close", (code) => (code === 0 ? (fs.mkdirSync(CACHE, { recursive: true }), fs.writeFileSync(cf, out), resolve(JSON.parse(out))) : reject(new Error(name + " exited " + code + "\n" + err))));
  });
}
if (!CASE) (async () => {
  const names = Object.keys(CASES), R = {};
  let next = 0;
  const workers = Array.from({ length: +(process.env.GATES_JOBS || 6) }, async () => {
    while (next < names.length) { const nm = names[next++]; R[nm] = await child(nm); console.error("done " + nm + " (" + R[nm].wallSeconds.toFixed(0) + " s)"); }
  });
  await Promise.all(workers);
  fs.mkdirSync(OUT, { recursive: true });
  const strip = (o) => { const c = Object.assign({}, o); delete c.series; return c; };
  const csv = (file, header, rows) => fs.writeFileSync(path.join(OUT, file), header + "\n" + rows.map((r) => r.join(",")).join("\n") + "\n");
  const rd = (x) => (x == null ? null : +x.toFixed(6));

  // Run 1
  const d0 = R.r1_none, d1 = R.r1_box, dd = d1.d_last5s - d0.d_last5s;
  const valid1 = d1.removedAtAdd === 0 && d1.obstacleDeleted === 0 && !d1.countChanged && !d0.countChanged && d1.nanMax === 0 && d0.nanMax === 0 && d1.insideMax === 0;
  const run1 = {
    gate: "|Δd − 0.0615| ≤ 0.0234 m (Δd in [0.0381, 0.0849])", valid: valid1, deltaD: dd, deltaD_fluidCellMethod: d1.dFluidCells_last5s - d0.dFluidCells_last5s,
    pass: valid1 && Math.abs(dd - 0.0615) <= 0.0234, theory_A_over_W: (336 * H * H) / (128 * H), none: strip(d0), box: strip(d1),
  };
  fs.writeFileSync(path.join(OUT, "displaced_depth_obst.json"), JSON.stringify(run1, null, 1));
  csv("displaced_depth_none_obst.csv", "t,d_surface_m,d_fluidcells_m,fluidCells", d0.series);
  csv("displaced_depth_box_obst.csv", "t,d_surface_m,d_fluidcells_m,fluidCells", d1.series);

  // Run 2
  const g2 = (r) => {
    const closureOk = Math.abs(r.closure) <= 1e-6;
    return {
      closureOk, velRemLE0: r.velRem <= 0, cumulativeLE1e4: r.pushTotal <= 1e-4, windowLE1e4: r.maxWindow05s <= 1e-4,
      pass: closureOk ? r.velRem <= 0 && r.pushTotal <= 1e-4 && r.maxWindow05s <= 1e-4 : "INVALID (closure)",
    };
  };
  const run2 = {
    gates: { velocityRemoval: "row alone ≤ 0", pushPlusRemoval: "cumulative ≤ +1e-4 J/kg AND no 0.5 s window above +1e-4", closure: "rows sum to measured ΔE within 1e-6 J/kg (else invalid)" },
    instrumentationIdentical: { box: R.r2_box.finalStateHash === R.r2_box_plain.finalStateHash, ramp: R.r2_ramp.finalStateHash === R.r2_ramp_plain.finalStateHash, hashes: [R.r2_box.finalStateHash, R.r2_box_plain.finalStateHash, R.r2_ramp.finalStateHash, R.r2_ramp_plain.finalStateHash] },
    cfl1: { none: strip(R.r2_none), box: Object.assign(strip(R.r2_box), { gate: g2(R.r2_box) }), ramp: Object.assign(strip(R.r2_ramp), { gate: g2(R.r2_ramp) }) },
    cfl05: { none: strip(R.r2_none_cfl05), box: Object.assign(strip(R.r2_box_cfl05), { gate: g2(R.r2_box_cfl05) }), ramp: Object.assign(strip(R.r2_ramp_cfl05), { gate: g2(R.r2_ramp_cfl05) }) },
  };
  for (const k of ["box", "ramp"]) {
    const a = R["r2_" + k], b = R["r2_" + k + "_cfl05"];
    run2["halving_" + k] = { pushPEKE_cfl1: a.pushPE + a.pushKE, pushPEKE_cfl05: b.pushPE + b.pushKE, ratioPEKE: (b.pushPE + b.pushKE) / (a.pushPE + a.pushKE), pushTotal_ratio: b.pushTotal / a.pushTotal, velRem_ratio: b.velRem / a.velRem, pushed_ratio: b.pushedTotal / a.pushedTotal };
  }
  fs.writeFileSync(path.join(OUT, "energy_contact_obst.json"), JSON.stringify(run2, null, 1));
  for (const nm of ["r2_box", "r2_ramp", "r2_box_cfl05", "r2_ramp_cfl05"]) csv(nm.replace("r2_", "energy_contact_series_") + "_obst.csv", "t_end,pushPE,pushKE_other,velRemoval,pushed", R[nm].series);

  // Run 3
  const g3 = (r) => r.maxDownstream === 0 && r.everDownstream === 0 && r.minSdfOverR >= 0;
  const run3 = {
    gate: "CFL 1, no path check: 1-cell and 2-cell barriers 0 leaks and min SDF ≥ 0; 0.5-cell, CFL 0.5 and ramp tip report-only",
    b1: Object.assign(R.r3_b1, { pass: g3(R.r3_b1) }), b2: Object.assign(R.r3_b2, { pass: g3(R.r3_b2) }),
    b05_centred: R.r3_b05c, b05_noCellCentre: R.r3_b05e, b1_cfl05: R.r3_b1_cfl05, noBarrierBaseline: R.r3_none, rampTip: R.r3_tip,
  };
  fs.writeFileSync(path.join(OUT, "thin_barrier_obst.json"), JSON.stringify(run3, null, 1));

  // Run 4
  const run4 = { reportOnly: true, ramp: strip(R.r4_ramp), freefall: strip(R.r4_freefall) };
  fs.writeFileSync(path.join(OUT, "ramp_slide_obst.json"), JSON.stringify(run4, null, 1));
  csv("ramp_slide_obst.csv", "t,s_mean_m,v_slope_mps,vy_mps,normal_offset_m", R.r4_ramp.series);
  csv("ramp_slide_freefall_obst.csv", "t,s_mean_m,v_slope_mps,vy_mps,normal_offset_m", R.r4_freefall.series);

  const summary = {
    run1: { pass: run1.pass, valid: run1.valid, deltaD: rd(dd), d_none: rd(d0.d_last5s), d_box: rd(d1.d_last5s) },
    run2: { box: run2.cfl1.box.gate, ramp: run2.cfl1.ramp.gate, instrumentationIdentical: run2.instrumentationIdentical },
    run3: { b1: run3.b1.pass, b2: run3.b2.pass },
  };
  fs.writeFileSync(path.join(OUT, "summary_obst.json"), JSON.stringify(summary, null, 1));
  console.log(JSON.stringify(summary, null, 1));
})().catch((e) => { console.error(e); process.exit(1); });
