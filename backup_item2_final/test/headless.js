// Headless check for sim.js (FLIP/APIC). Not referenced by index.html.
// Usage: node test/headless.js [widthPx heightPx]
const path = require("path");
const { FlipSolver: solver, computeLayout } = require(path.join(__dirname, "..", "sim.js"));
// thin adapter over the solver API
const sim = {
  params: solver.params,
  initWorld: (w, h) => { const L = computeLayout(w, h); solver.init(L.cellsX, L.cellsY, L.cellSize); },
  resetSim: (m) => solver.reset(m),
  simulateFrame: () => solver.step(),
  getStats: () => { const s = solver.stats(); s.grid = { nx: solver.nx - 2, ny: solver.ny - 2, h: solver.cellSize }; return s; },
  getParticles: () => ({ count: solver.count, px: solver.px, py: solver.py, pvx: solver.vx, pvy: solver.vy, h: solver.cellSize, nx: solver.nx, ny: solver.ny, pr: solver.radius }),
  setMouse: (x, y, vx, vy) => solver.setPointer(x, y, vx, vy),
  spawnHeldDrop: (x, y) => solver.spawnDrop(x, y),
  releaseHeldDrop: () => solver.releaseDrop(),
};

const W = parseInt(process.argv[2] || "1600", 10);
const H = parseInt(process.argv[3] || "900", 10);
if (process.argv[4]) sim.params.flipRatio = parseFloat(process.argv[4]);
if (process.argv[5]) sim.params.viscosity = parseFloat(process.argv[5]);
sim.initWorld(W, H);
sim.resetSim("dam");

const fps = 60;
function surfaceProfile() {
  // top-most particle per grid column, only counting columns with >= 3 particles in a cell stack
  const P = sim.getParticles();
  const cols = P.nx - 2;
  const colCount = new Int32Array(cols);
  const colTop = new Float32Array(cols).fill(1e9);
  // per column, count particles then estimate water depth = count / (rest per cell) * h
  for (let i = 0; i < P.count; i++) {
    const c = Math.min(cols - 1, Math.max(0, Math.floor(P.px[i] / P.h) - 1));
    colCount[c]++;
  }
  return colCount;
}

function fmt(n, d = 3) {
  return Number(n).toFixed(d);
}

let s0 = null;
const t0Stats = [];
let frontHitT = -1;
let maxRunup = 0;
let comXs = [];
let totalMs = 0;
let maxMs = 0;
let frames = 0;
let worstNan = 0;
const frameTimes = [];
let substepTotal = 0;
let worstOut = 0;
let globalVmax = 0;

function frame() {
  const t = performance.now();
  sim.simulateFrame();
  const dt = performance.now() - t;
  totalMs += dt;
  if (dt > maxMs) maxMs = dt;
  frameTimes.push(dt);
  substepTotal += solver.stats().substeps;
  frames++;
}

// warm-up first frame measures rest density
frame();
s0 = sim.getStats();
const P = sim.getParticles();
const tankW = (P.nx - 2) * P.h;
const tankH = (P.ny - 2) * P.h;
const startVolume = s0.fluidCells;
console.log(`canvas ${W}x${H}  grid ${s0.grid.nx}x${s0.grid.ny} (h=${fmt(s0.grid.h, 4)} m)  tank ${fmt(tankW, 2)}x${fmt(tankH, 2)} m`);
console.log(`particles ${s0.count}  restDensity ${fmt(s0.restDensity, 2)} /cell  flipRatio ${sim.params.flipRatio}  iters ${sim.params.pressureIters}`);

const depthAtRest = (s0.fluidCells * P.h * P.h) / tankW;
const samples = [];
for (let f = 1; f <= 30 * fps; f++) {
  frame();
  const st = sim.getStats();
  worstNan = Math.max(worstNan, st.nan);
  worstOut = Math.max(worstOut, st.outside);
  globalVmax = Math.max(globalVmax, st.vmax);
  const Pp = sim.getParticles();
  // wave front & run-up at right wall
  let maxX = 0;
  let comX = 0;
  let runTop = 1e9;
  for (let i = 0; i < Pp.count; i++) {
    if (Pp.px[i] > maxX) maxX = Pp.px[i];
    comX += Pp.px[i];
    if (Pp.px[i] > (Pp.nx - 1) * Pp.h - 3 * Pp.h && Pp.py[i] < runTop) runTop = Pp.py[i];
  }
  comX /= Pp.count;
  comXs.push(comX);
  if (frontHitT < 0 && maxX > (Pp.nx - 1) * Pp.h - Pp.pr - 0.01) frontHitT = f / fps;
  const runup = (Pp.ny - 1) * Pp.h - runTop;
  if (f < 6 * fps && runup > maxRunup && runTop < 1e8) maxRunup = runup;
  if (f % (fps * 2) === 0 || f === fps * 10) {
    samples.push({ t: f / fps, ...st, comX });
  }
}

console.log("\n t(s)  fluidCells  vol%   interiorRho/rest  vmax(m/s)  KE/particle  substeps  comX(m)");
for (const s of samples) {
  console.log(
    `${fmt(s.t, 0).padStart(5)}  ${String(s.fluidCells).padStart(10)}  ${fmt((100 * s.fluidCells) / startVolume, 1).padStart(5)}  ${fmt(s.interiorDensity / s.restDensity, 3).padStart(16)}  ${fmt(s.vmax, 2).padStart(9)}  ${fmt(s.ke, 4).padStart(11)}  ${String(s.substeps).padStart(8)}  ${fmt(s.comX, 3)}`
  );
}
// slosh: count center-of-mass reversals after the front hit
let reversals = 0;
let lastSign = 0;
for (let k = 2; k < comXs.length; k++) {
  const d = comXs[k] - comXs[k - 2];
  if (Math.abs(d) < 1e-4) continue;
  const sgn = Math.sign(d);
  if (lastSign && sgn !== lastSign) reversals++;
  lastSign = sgn;
}
// flatness at 30 s: per-column particle counts -> depth
const cc = surfaceProfile();
let mean = 0;
for (const c of cc) mean += c;
mean /= cc.length;
let varr = 0;
for (const c of cc) varr += (c - mean) * (c - mean);
const std = Math.sqrt(varr / cc.length);
const st30 = sim.getStats();
console.log(`\nwave front hit far wall at t=${fmt(frontHitT, 2)} s; max run-up at right wall ${fmt(maxRunup, 2)} m (initial column ${fmt(tankH * 0.58, 2)} m, rest depth ~${fmt(depthAtRest * 0.0 + (st30.fluidCells * P.h * P.h) / tankW, 2)} m)`);
console.log(`center-of-mass direction reversals over 30 s (sloshing): ${reversals}`);
console.log(`flatness at 30 s: per-column particle count mean ${fmt(mean, 1)}, std ${fmt(std, 2)} (${fmt((100 * std) / mean, 1)}%)`);
console.log(`volume (fluid cells) start ${startVolume} -> 10 s ${samples.find((s) => s.t === 10).fluidCells} -> 30 s ${st30.fluidCells}  (${fmt((100 * (st30.fluidCells - startVolume)) / startVolume, 1)}% at 30 s)`);
console.log(`interior density / rest at 30 s: ${fmt(st30.interiorDensity / st30.restDensity, 3)}`);
console.log(`NaN max ${worstNan}, outside-box max ${worstOut}, max speed ${fmt(globalVmax, 2)} m/s`);
{
  const ft = frameTimes.slice(60).sort((a, b) => a - b); // skip 1 s of JIT warm-up
  const q = (x) => ft[Math.min(ft.length - 1, Math.floor(x * ft.length))];
  console.log(`dam phase timing (30 s): ${fmt(totalMs / frames, 2)} ms/frame avg, p50 ${fmt(q(0.5), 2)}, p95 ${fmt(q(0.95), 2)}, max after warm-up ${fmt(ft[ft.length - 1], 1)} ms; ${fmt(totalMs / substepTotal, 2)} ms/substep; avg ${fmt(substepTotal / frames, 2)} substeps/frame`);
  frameTimes.length = 0;
}

// --- drop test: one held drop released from up high into the settled pool
totalMs = 0; maxMs = 0; frames = 0;
sim.setMouse(tankW * 0.5 + P.h, tankH * 0.2 + P.h, 0, 0);
sim.spawnHeldDrop(tankW * 0.5 + P.h, tankH * 0.2 + P.h);
for (let f = 0; f < 20; f++) frame();
const nBeforeDrop = sim.getStats().count;
sim.setMouse(tankW * 0.5 + P.h, tankH * 0.2 + P.h, 0, 3);
sim.releaseHeldDrop();
const surfY = (P.ny - 1) * P.h - (st30.fluidCells * P.h * P.h) / tankW;
let maxAirborne = 0;
let maxSplashH = 0;
for (let f = 0; f < 3 * fps; f++) {
  frame();
  const Pp = sim.getParticles();
  if (f < 25) continue; // let it hit first
  let air = 0;
  let top = 1e9;
  for (let i = 0; i < Pp.count; i++) {
    if (Pp.py[i] < surfY - 3 * Pp.h) {
      air++;
      if (Pp.py[i] < top) top = Pp.py[i];
    }
  }
  if (air > maxAirborne) maxAirborne = air;
  if (top < 1e8 && surfY - top > maxSplashH) maxSplashH = surfY - top;
}
const sd = sim.getStats();
console.log(`\ndrop: ${nBeforeDrop - st30.count} particles; after impact up to ${maxAirborne} particles >3 cells above the pool surface, splash height ${fmt(maxSplashH, 2)} m; NaN ${sd.nan}, outside ${sd.outside}, vmax ${fmt(sd.vmax, 2)}`);

// --- spam test: 30 drops thrown hard in quick succession
let spamWorstNan = 0, spamWorstOut = 0, spamVmax = 0;
for (let d = 0; d < 30; d++) {
  const cx = P.h + tankW * (0.1 + 0.8 * Math.random());
  const cy = P.h + tankH * (0.1 + 0.3 * Math.random());
  sim.setMouse(cx, cy, 0, 0);
  sim.spawnHeldDrop(cx, cy);
  for (let f = 0; f < 3; f++) frame();
  sim.setMouse(cx, cy, (Math.random() - 0.5) * 20, (Math.random() - 0.3) * 20);
  sim.releaseHeldDrop();
  for (let f = 0; f < 6; f++) {
    frame();
    const st = sim.getStats();
    spamWorstNan = Math.max(spamWorstNan, st.nan);
    spamWorstOut = Math.max(spamWorstOut, st.outside);
    spamVmax = Math.max(spamVmax, st.vmax);
  }
}
for (let f = 0; f < 5 * fps; f++) {
  frame();
  const st = sim.getStats();
  spamWorstNan = Math.max(spamWorstNan, st.nan);
  spamWorstOut = Math.max(spamWorstOut, st.outside);
  spamVmax = Math.max(spamVmax, st.vmax);
}
const ss = sim.getStats();
console.log(`spam 30 thrown drops: particles ${ss.count}, NaN max ${spamWorstNan}, outside max ${spamWorstOut}, max speed ${fmt(spamVmax, 2)} m/s, interior rho/rest ${fmt(ss.interiorDensity / ss.restDensity, 3)}`);
console.log(`substep-cap hits (8/frame): ${solver.capLog.length}` + (solver.capLog.length ? ` — e.g. ${solver.capLog.slice(0, 5).map((c) => `frame ${c.frame} u_max ${c.umax.toFixed(2)}`).join("; ")}` : ""));
if (solver.capLog.length) {
  const cov = solver.capLog.map((c) => c.covered);
  const allShort = solver.capLog.every((c) => c.covered < solver.FIXED_DT - 1e-12 && c.substeps === solver.MAX_SUBSTEPS);
  console.log(`  cap frames: sim time covered ${(Math.min(...cov) * 1000).toFixed(3)}–${(Math.max(...cov) * 1000).toFixed(3)} ms (< ${(1000 / 60).toFixed(3)} ms in every capped frame: ${allShort}); max CFL inside capped frames ${Math.max(...solver.capLog.map((c) => c.maxCfl)).toFixed(4)}`);
}
console.log(`max advective CFL (u_max*dt/h) over every substep of the whole run: ${solver.maxCflObserved.toFixed(4)}; PCG safety-cap hits: ${solver.pcgCapLog.length}`);
console.log(`drop/spam timing: ${fmt(totalMs / frames, 2)} ms/frame avg, worst ${fmt(maxMs, 1)} ms`);
const ok = worstNan === 0 && worstOut === 0 && spamWorstNan === 0 && spamWorstOut === 0 && Math.abs(st30.fluidCells - startVolume) / startVolume < 0.1;
console.log(ok ? "\nPASS" : "\nFAIL");
process.exit(ok ? 0 : 1);
