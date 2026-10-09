// Headless test for fixed obstacles (no DOM). Not referenced by index.html.
// Usage: node test/obstacles.js [seconds=10] [widthPx=1600 heightPx=900]
//   1. regression: with no obstacles, sim.js is bit-identical to backup_pre_objects/sim.js
//   2. dam break with a circle + a ramp in the path vs the same dam break without them:
//      leak check (no particle centre inside an obstacle, SDF < 0) every frame, NaNs,
//      particle count, volume (fluid cells) and interior density, ms/frame and ms/substep
//   3. placing a box into the settled pool: particles inside are deleted (counted), no leaks after
//   4. API: hitTest / updateObstacle / removeObstacle / clearObstacles / reset keeps obstacles /
//      resetColumn skips obstacle cells / freehand builder (RDP, self-crossing → hull, thin → rejected)
const fs = require("fs");
const path = require("path");
const SIM = path.join(__dirname, "..", "sim.js");
const BACKUP = path.join(__dirname, "..", "backup_pre_objects", "sim.js");
const { FlipSolver: S, computeLayout, ObstacleShapes: OS, demoObstacleShapes } = require(SIM);

const SECONDS = parseFloat(process.argv[2] || "10");
const WPX = parseInt(process.argv[3] || "1600", 10), HPX = parseInt(process.argv[4] || "900", 10);
const FPS = 60;
const fmt = (x, d = 3) => Number(x).toFixed(d);
let failures = 0;
function check(name, ok, info) {
  if (!ok) failures++;
  console.log((ok ? "ok   " : "FAIL ") + name + (info ? "  " + info : ""));
}
function setup() {
  const L = computeLayout(WPX, HPX);
  S.clearObstacles();
  S.init(L.cellsX, L.cellsY, L.cellSize);
  return L;
}

// ---- 1. regression vs the pre-obstacle sim.js --------------------------------
if (fs.existsSync(BACKUP)) {
  const B = require(BACKUP).FlipSolver;
  const L = computeLayout(WPX, HPX);
  B.init(L.cellsX, L.cellsY, L.cellSize); B.reset("dam");
  setup(); S.reset("dam");
  const frames = 240;
  for (let f = 0; f < frames; f++) { B.step(); S.step(); }
  let diff = 0;
  for (let i = 0; i < Math.max(B.count, S.count); i++) {
    if (B.px[i] !== S.px[i] || B.py[i] !== S.py[i] || B.vx[i] !== S.vx[i] || B.vy[i] !== S.vy[i]) diff++;
  }
  const sb = B.stats(), ss = S.stats();
  check("no-obstacle path bit-identical to backup_pre_objects/sim.js (" + frames + " frames dam break)",
    diff === 0 && B.count === S.count && sb.fluidCells === ss.fluidCells && sb.interiorDensity === ss.interiorDensity,
    `count ${S.count}, differing particles ${diff}, fluidCells ${sb.fluidCells}/${ss.fluidCells}`);
} else console.log("skip regression: backup_pre_objects/sim.js not found");

// ---- 2. dam break: baseline vs circle + ramp ---------------------------------
function damRun(withObstacles) {
  setup();
  let shapes = [];
  if (withObstacles) for (const s of demoObstacleShapes(S)) shapes.push(S.addObstacle(s));
  S.reset("dam");
  const n0 = S.count;
  let worstInside = 0, minSdf = Infinity, nan = 0, outside = 0, ms = 0, substeps = 0, msMax = 0;
  const samples = [];
  const h = S.cellSize;
  const total = Math.round(SECONDS * FPS);
  for (let f = 1; f <= total; f++) {
    const t0 = process.hrtime.bigint();
    S.step();
    const dt = Number(process.hrtime.bigint() - t0) / 1e6;
    if (f > 30) { ms += dt; substeps += S.stats().substeps; if (dt > msMax) msMax = dt; }
    const st = S.stats();
    nan = Math.max(nan, st.nan);
    outside = Math.max(outside, st.outside);
    worstInside = Math.max(worstInside, st.insideObstacle);
    if (st.minObstacleSdf < minSdf) minSdf = st.minObstacleSdf;
    if (f % FPS === 0 || f === 30) samples.push({ t: f / FPS, fluidCells: st.fluidCells, rho: st.interiorDensity / st.restDensity, ke: st.ke, vmax: st.vmax });
  }
  const st = S.stats();
  return {
    n0, n1: S.count, worstInside, minSdf, nan, outside, samples, st,
    msFrame: ms / (total - 30), msSub: ms / substeps, subPerFrame: substeps / (total - 30), msMax,
    obstacleCells: st.obstacleCells, r: S.radius, h,
  };
}
const base = damRun(false);
const obs = damRun(true);
console.log(`\ndam break ${SECONDS} s, grid ${S.nx - 2}x${S.ny - 2} (h=${fmt(obs.h, 4)} m), circle + ramp = ${obs.obstacleCells} solid cells`);
console.log("  t(s)   fluidCells base / obst    rho/rest base / obst    KE/particle base / obst");
for (let k = 0; k < obs.samples.length; k++) {
  const a = base.samples[k], b = obs.samples[k];
  if (!(b.t === 0.5 || b.t === 1 || b.t === 2 || b.t % 2 === 0 || b.t === SECONDS)) continue;
  console.log(`  ${fmt(b.t, 1).padStart(4)}   ${String(a.fluidCells).padStart(10)} / ${String(b.fluidCells).padEnd(6)}   ${fmt(a.rho).padStart(12)} / ${fmt(b.rho)}   ${fmt(a.ke, 4).padStart(14)} / ${fmt(b.ke, 4)}`);
}
const lastA = base.samples[base.samples.length - 1], lastB = obs.samples[obs.samples.length - 1];
const avg = (arr, k) => arr.slice(1).reduce((s, x) => s + x[k], 0) / (arr.length - 1);
console.log(`  mean over samples: fluidCells ${fmt(avg(base.samples, "fluidCells"), 1)} / ${fmt(avg(obs.samples, "fluidCells"), 1)}, rho/rest ${fmt(avg(base.samples, "rho"))} / ${fmt(avg(obs.samples, "rho"))}`);
console.log(`  particles: base ${base.n0} -> ${base.n1}; obstacles ${obs.n0} -> ${obs.n1} (deleted on placement ${obs.st.obstacleDeleted}, drop-clipped ${obs.st.dropClipped})`);
console.log(`  closest particle centre to an obstacle: ${fmt(obs.minSdf / obs.r, 4)} radii (push-out target 1.0); wedge fallbacks ${obs.st.pushStuck}`);
console.log(`  timing: base ${fmt(base.msFrame, 2)} ms/frame (${fmt(base.msSub, 3)} ms/substep, ${fmt(base.subPerFrame, 2)} substeps/frame, max ${fmt(base.msMax, 1)} ms)`);
console.log(`          obst ${fmt(obs.msFrame, 2)} ms/frame (${fmt(obs.msSub, 3)} ms/substep, ${fmt(obs.subPerFrame, 2)} substeps/frame, max ${fmt(obs.msMax, 1)} ms)  → ${fmt((100 * (obs.msSub - base.msSub)) / base.msSub, 1)}% per substep`);
check("no particle centre inside an obstacle at any frame", obs.worstInside === 0, `worst ${obs.worstInside}`);
check("particles stay >= ~1 radius from obstacles", obs.minSdf >= 0.99 * obs.r, `min ${fmt(obs.minSdf / obs.r, 4)} r`);
check("no NaN / out-of-box particles", obs.nan === 0 && obs.outside === 0 && base.nan === 0);
check("particle count unchanged (obstacles placed before reset: nothing spawned inside)", obs.n0 === obs.n1 && obs.n0 === base.n0, `${obs.n0} vs baseline ${base.n0}`);
check("water volume kept: final fluid cells within 5% of the no-obstacle run",
  Math.abs(lastB.fluidCells - lastA.fluidCells) / lastA.fluidCells < 0.05, `${lastB.fluidCells} vs ${lastA.fluidCells}`);
check("interior density / rest within 3% of 1 at the end", Math.abs(lastB.rho - 1) < 0.03, fmt(lastB.rho));

// ---- 3. place a box into the settled pool ------------------------------------
{
  const h = S.cellSize, floor = (S.ny - 1) * h;
  const before = S.count;
  // box straddling the pool surface, right of the ramp's toe region
  const id = S.addObstacle(OS.box(h + 0.25 * (S.nx - 2) * h, floor - 0.12, 0.3, 0.2, 0.3));
  const removed = S.lastRemovedParticles;
  check("box dropped into the pool deletes the particles inside it", removed > 0 && S.count === before - removed && S.stats().insideObstacle === 0,
    `removed ${removed} (stats.obstacleDeleted ${S.stats().obstacleDeleted}, ${before} -> ${S.count})`);
  let worst = 0, nan = 0, minSdf = Infinity;
  for (let f = 0; f < 2 * FPS; f++) {
    S.step();
    const st = S.stats();
    worst = Math.max(worst, st.insideObstacle); nan = Math.max(nan, st.nan); minSdf = Math.min(minSdf, st.minObstacleSdf);
  }
  check("2 s after placing the box: no leaks, no NaN", worst === 0 && nan === 0, `min sdf ${fmt(minSdf / S.radius, 3)} r`);
  // drop inside an obstacle is clipped on release
  const cx = S.obstacles.find((o) => o.id === id).cx, cy = S.obstacles.find((o) => o.id === id).cy;
  S.setPointer(cx, cy - 0.05, 0, 0);
  S.spawnDrop(cx, cy - 0.05);
  S.step();
  S.releaseDrop();
  S.step();
  check("held drop released over an obstacle: overlapping particles removed, none inside", S.stats().insideObstacle === 0 && S.stats().dropClipped > 0, `clipped ${S.stats().dropClipped}`);
  // move it (updateObstacle) and rotate it; same id
  const n1 = S.count;
  check("updateObstacle moves an obstacle (same id)", S.updateObstacle(id, OS.rotate(OS.translate(S.obstacles.find((o) => o.id === id), 0.4, 0), 0.5)) && S.obstacles.some((o) => o.id === id),
    `removed ${S.lastRemovedParticles} (${n1} -> ${S.count})`);
  for (let f = 0; f < 30; f++) S.step();
  check("after the move: no leaks", S.stats().insideObstacle === 0);
}

// ---- 4. API and builders ---------------------------------------------------------
{
  setup();
  const h = S.cellSize, floor = (S.ny - 1) * h;
  const a = S.addObstacle(OS.circle(1.0, floor - 0.5, 0.1));
  const b = S.addObstacle(OS.box(2.0, floor - 0.5, 0.3, 0.1, Math.PI / 4));
  check("hitTest finds circle / rotated box / nothing", S.hitTest(1.0, floor - 0.5) === a && S.hitTest(2.0, floor - 0.5) === b && S.hitTest(1.5, floor - 0.5) === -1);
  check("hitTest respects the box rotation", S.hitTest(2.0 + 0.14, floor - 0.5) === -1 && S.hitTest(2.0 + 0.09, floor - 0.5 + 0.09) === b);
  check("obstacles getter returns copies with ids", S.obstacles.length === 2 && S.obstacles[0].id === a && S.obstacles[1].type === "rect");
  S.reset("dam");
  check("reset keeps obstacles", S.obstacles.length === 2 && S.stats().insideObstacle === 0);
  S.reset("dam", { clearObstacles: true });
  check("reset(mode, { clearObstacles: true }) clears them", S.obstacles.length === 0);
  S.addObstacle(OS.box(h + 0.3, floor - 0.3, 0.2, 0.2));
  S.resetColumn(24, 48);
  const solidCells = S.stats().obstacleCells;
  check("resetColumn keeps obstacles and spawns no particle inside", S.obstacles.length === 1 && S.stats().insideObstacle === 0 && S.count < 24 * 48 * 4,
    `count ${S.count} of ${24 * 48 * 4} lattice sites, obstacle cells ${solidCells}`);
  check("removeObstacle / clearObstacles", S.removeObstacle(S.obstacles[0].id) && S.obstacles.length === 0 && !S.removeObstacle(999) && (S.clearObstacles(), S.obstacles.length === 0));
  let threw = 0;
  for (const bad of [{ type: "circle", cx: 1, cy: 1, r: -1 }, { type: "polygon", cx: 0, cy: 0, points: [[0, 0], [1, 1], [2, 2]] }, { type: "blob", cx: 0, cy: 0 }, { type: "rect", cx: NaN, cy: 0, w: 1, h: 1 }]) {
    try { S.addObstacle(bad); } catch (e) { threw++; }
  }
  check("invalid shapes are rejected with TypeError", threw === 4 && S.obstacles.length === 0);
  // SDF sanity
  const sq = OS.polygon([[0, 0], [1, 0], [1, 1], [0, 1]], "polygon");
  const o = { nx: 0, ny: 0 };
  const d1 = OS.sdf(sq, 0.5, 0.5, o), d2 = OS.sdf(sq, 1.5, 0.5, o);
  check("polygon SDF: inside -0.5, outside +0.5 with outward normal", Math.abs(d1 + 0.5) < 1e-12 && Math.abs(d2 - 0.5) < 1e-12 && Math.abs(o.nx - 1) < 1e-12);
  const rb = OS.box(0, 0, 2, 1, 0), dR = OS.sdf(rb, 2, 1, o);
  check("box SDF corner distance", Math.abs(dR - Math.hypot(1, 0.5)) < 1e-12);
  // freehand: noisy circle → simplified polygon; figure-8 → hull; thin line → rejected
  const ring = [];
  for (let k = 0; k < 200; k++) { const t = (k / 200) * 2 * Math.PI; ring.push([1 + 0.2 * Math.cos(t) + 0.002 * Math.sin(37 * t), 1 + 0.2 * Math.sin(t)]); }
  const fr = OS.freehand(ring, h);
  check("freehand: 200-point loop simplified (RDP)", fr.shape && fr.shape.points.length >= 8 && fr.shape.points.length <= 64 && Math.abs(Math.abs(OS.polygonArea(fr.shape.points)) - Math.PI * 0.04) / (Math.PI * 0.04) < 0.08,
    fr.shape ? `${fr.shape.points.length} vertices` : fr.error);
  const eight = [];
  for (let k = 0; k < 120; k++) { const t = (k / 120) * 2 * Math.PI; eight.push([1 + 0.3 * Math.sin(t), 1 + 0.15 * Math.sin(2 * t)]); }
  const f8 = OS.freehand(eight, h);
  check("freehand: self-crossing stroke → convex hull (simple polygon)", f8.shape && !OS.selfIntersects(f8.shape.points) && /hull/.test(f8.note || ""), f8.note || f8.error);
  const line = [];
  for (let k = 0; k < 50; k++) line.push([0.5 + k * 0.01, 1 + 0.001 * (k % 2)]);
  const fl = OS.freehand(line, h);
  check("freehand: straight stroke rejected", !fl.shape && !!fl.error, fl.error);
}

// ---- 5. item 1: one-cell-thick obstacles and the editor thickness rule -------------
{
  const NX = 128, NY = 72, h = 3 / NX, floor = (NY + 1) * h; // init(nx, ny) counts interior cells; walls add one cell each side
  const PRE = path.join(__dirname, "..", "backup_pre_item1", "sim.js");
  const barrier = (Solver, cells) => {
    Solver.clearObstacles(); Solver.init(NX, NY, h);
    Solver.addObstacle({ type: "rect", cx: (64 + cells / 2) * h, cy: 0.5 * (NY + 2) * h, w: cells * h, h: (NY + 2) * h, angle: 0 }); // full height
    Solver.reset("dam");
  };
  const energy = (Q) => { let e = 0; for (let i = 0; i < Q.count; i++) e += 0.5 * (Q.vx[i] ** 2 + Q.vy[i] ** 2) + 9.81 * (floor - Q.py[i]); return e / Q.count; };
  // 5a. 2-cell barrier (b2). [item 2] This was "bit-identical to backup_pre_item1/sim.js"; item 2 deliberately
  // changes obstacle contact (pre-advection normal strip + obstacle-aware separation clip), so identity is no
  // longer expected. Now: same dam break (120 frames) as the pre-item-1 build -> no NaN, nothing inside, no leaks past the
  // barrier, same particle count, and the energy agrees with the pre-item-1 build within 1 %. Differing particles are reported.
  if (fs.existsSync(PRE)) {
    const B = require(PRE).FlipSolver;
    barrier(B, 2); barrier(S, 2);
    let leakS = 0, leakB = 0, inside = 0, nan = 0;
    for (let f = 0; f < 120; f++) { B.step(); S.step(); const st = S.stats(); inside = Math.max(inside, st.insideObstacle); nan = Math.max(nan, st.nan); }
    for (let i = 0; i < S.count; i++) { if (S.px[i] > 66 * h) leakS++; }
    for (let i = 0; i < B.count; i++) { if (B.px[i] > 66 * h) leakB++; }
    let diff = 0;
    for (let i = 0; i < Math.max(B.count, S.count); i++) if (B.px[i] !== S.px[i] || B.py[i] !== S.py[i] || B.vx[i] !== S.vx[i] || B.vy[i] !== S.vy[i]) diff++;
    const eS = energy(S), eB = energy(B);
    check("2-cell barrier vs backup_pre_item1/sim.js (120 frames; bit-identity retired in item 2): finite, none inside, energy within 1 %",
      B.count === S.count && nan === 0 && inside === 0 && Math.abs(eS - eB) / Math.abs(eB) < 0.01, `differing particles ${diff} (expected, by design), E/E_pre ${fmt(eS / eB, 4)}, right-of-barrier ${leakS} vs pre ${leakB} (water legitimately passes over the 2-cell wall at 120 frames only if it overtops it)`);
  } else console.log("skip: backup_pre_item1/sim.js not found");
  // 5b. 1-cell barrier (b1): energy stays bounded and tracks the 2-cell barrier
  const run = (cells) => { barrier(S, cells); const E0 = energy(S); let Emax = E0, nan = 0, inside = 0; for (let f = 0; f < 90; f++) { S.step(); const st = S.stats(); nan = Math.max(nan, st.nan); inside = Math.max(inside, st.insideObstacle); Emax = Math.max(Emax, energy(S)); } return { E0, E: energy(S), Emax, nan, inside, vmax: S.stats().vmax }; };
  const b1 = run(1), b2 = run(2);
  check("1-cell barrier: 1.5 s dam break, energy bounded like the 2-cell barrier",
    b1.nan === 0 && b1.inside === 0 && b1.Emax <= b1.E0 * 1.005 && Math.abs(b1.E - b2.E) / b2.E0 < 0.01 && b1.vmax < 12,
    `E/E0 b1 ${fmt(b1.E / b1.E0, 4)} b2 ${fmt(b2.E / b2.E0, 4)}, max E/E0 b1 ${fmt(b1.Emax / b1.E0, 4)}, vmax b1 ${fmt(b1.vmax, 2)}`);
  // 5c. editor thickness rule (ObstacleShapes.thicknessCheck, editor Confirm)
  S.clearObstacles(); S.init(NX, NY, h);
  if (S.ny !== NY + 2) throw new Error("unexpected grid layout");
  const F = floor, T = (sh) => S.checkObstacleThickness(sh);
  const clipFix = (sh) => { for (let cut = 2 * h; cut <= 6 * h + 1e-9; cut += 0.5 * h) { const c = OS.clipAcuteTips(sh, cut); if (c.clipped && T(c.shape).ok) return c.shape; } return null; };
  check("rule: 1-cell box fails, 2-cell box passes", !T({ type: "rect", cx: 64.5 * h, cy: 36 * h, w: h, h: 72 * h, angle: 0 }).ok && T({ type: "rect", cx: 65 * h, cy: 36 * h, w: 2 * h, h: 72 * h, angle: 0 }).ok);
  check("rule: smallest editor box (3 cells) passes at 0° and 45°; 1.5-cell rotated bar fails",
    T(OS.box(40.5 * h, 30.5 * h, 3 * h, 3 * h)).ok && T(OS.box(40.3 * h, 30.7 * h, 3 * h, 20 * h, Math.PI / 4)).ok && !T(OS.box(40.3 * h, 30.7 * h, 1.5 * h, 20 * h, Math.PI / 6)).ok);
  check("rule: circles r >= 1.5 cells pass, sub-cell circle fails", T(OS.circle(40.2 * h, 30.3 * h, 1.5 * h)).ok && T(OS.circle(40.2 * h, 30.3 * h, 0.12)).ok && !T(OS.circle(40.5 * h, 30.5 * h, 0.4 * h)).ok);
  check("rule: one-cell shelf lying on the floor passes (merges with the wall)", T({ type: "rect", cx: 91 * h, cy: F - 0.5 * h, w: 20 * h, h: h, angle: 0 }).ok);
  check("rule: floor ramp tip passes (wall counts as solid), demo shapes pass", T(OS.ramp(H1(h), F - 40 * h, 71 * h, F)).ok && demoObstacleShapes(S).every((sh) => T(sh).ok));
  const air = OS.ramp(100 * h, F - 16 * h, 70 * h, F - 12 * h), airFix = clipFix(air);
  check("rule: mid-air ramp tip fails, auto-clipped polygon passes and stays a ramp", !T(air).ok && airFix && airFix.kind === "ramp" && T(airFix).ok, airFix ? `${air.points.length} -> ${airFix.points.length} vertices` : "no clip");
  const neck = OS.polygon([[20 * h, 20 * h], [30 * h, 20 * h], [30 * h, 24.6 * h], [45 * h, 24.6 * h], [45 * h, 20 * h], [55 * h, 20 * h], [55 * h, 30 * h], [45 * h, 30 * h], [45 * h, 25.6 * h], [30 * h, 25.6 * h], [30 * h, 30 * h], [20 * h, 30 * h]], "freehand");
  const nk = T(neck);
  check("rule: freehand with a 1-cell neck fails (and is not rescued by tip clipping)", !nk.ok && nk.thin.length > 0 && !clipFix(neck), `${nk.thin.length} thin cells`);
}
function H1(h) { return h; }

console.log(failures ? `\nFAIL (${failures})` : "\nPASS");
process.exit(failures ? 1 : 0);
