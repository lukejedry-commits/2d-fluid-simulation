// Item 1 close-out runs (no solver changes; sim.js is loaded from an OS-temp copy with harness patches).
// Usage: node test/item1_closeout.js            -> runs every case below in parallel child processes
//        node test/item1_closeout.js --case X   -> one case, JSON on stdout (fresh process)
// Common: 128x72 grid (h = 3/128 m), damping 0, speed clamp OFF (MAX_SPEED = Infinity), NaN guard at
// the start of every substep, CFL 1, reviewed defaults, 30 s.
// Cases:
//   base_dam      no-obstacle dam break (baseline for validation/item1_energy.png)
//   tilt_b1       API run: 1-cell full-height solid at column 64 added with addObstacle (editor rule not
//                 applied). Water on both sides, each side tilted about its OWN centre with the SAME
//                 slope (same direction), equal mean depth d = 16 cells on both sides; amplitude 2 cells
//                 at the outer ends of the wider (right, 64-cell) side, 1.97 cells on the left (63 cells).
//   tilt_none     the same initial particle positions without the barrier (comparison for max |u|)
//   tilt_b1_eqN   extra: same slope, right-side mean depth lowered so the particle counts match exactly
//                 (the sides are 63 and 64 cells wide, so equal depth means 1.6 % more particles right)
// Seeding (harness patch: resetColumn(.., .., pts) adds the given particles at rest): per side, a
// lattice with ~2x2 particles per cell: x evenly spaced in [side start + r, side end - r] (r = particle
// radius, so no particle touches a solid and none is deleted), rows y = floor - r - k h/2; a point is
// fluid if its height above the floor is below the tilted surface eta(x) = d + s (x - x_centre).
// Measured: E(t) per frame; per substep: net horizontal force on the barrier (pressure of the final
// active-set pass, p = P rho h / dt, rho = 1000, on the wet faces of columns 63 and 65 unless the face
// is released), leaks (particles that ever end up on the other side), min SDF / r, max particle |u|,
// thin-ghost statistics (v ghosts inside the barrier whose two copy-source faces both have P2G mass,
// and the sign of the two source velocities at the last fill of each substep).
const fs = require("fs"), os = require("os"), path = require("path"), { spawn } = require("child_process");
const ROOT = path.join(__dirname, ".."), SIM = path.join(ROOT, "sim.js"), V = path.join(ROOT, "validation");
const NX = 128, NY = 72, H = 3 / 128, G = 9.81, RHO = 1000, T_END = +(process.env.TEND || 30), FLOOR = (NY + 1) * H;
const CASE = (() => { const i = process.argv.indexOf("--case"); return i >= 0 ? process.argv[i + 1] : null; })();
const D_CELLS = 16, SLOPE = 2 / 32; // mean depth 16 cells (0.375 m); slope: 2 cells over a 32-cell half width

function load() {
  let src = fs.readFileSync(SIM, "utf8");
  const patch = (a, b) => { if (src.split(a).length !== 2) throw new Error("anchor: " + a.slice(0, 70)); src = src.replace(a, b); };
  patch("  const MAX_SPEED = 14;", "  const MAX_SPEED = Infinity; // [closeout] clamp OFF");
  patch("  function substep(dt) {\n    pinHeld();\n", "  function substep(dt) {\n    for (let i = 0; i < count; i++) if (!(Number.isFinite(px[i]) && Number.isFinite(py[i]) && Number.isFinite(pvx[i]) && Number.isFinite(pvy[i]))) throw new Error(\"NaN guard: particle \" + i + \" at t=\" + simTime);\n    pinHeld();\n");
  patch("  function resetColumn(cellsA, cellsH) {", "  function resetColumn(cellsA, cellsH, pts) {");
  patch("    for (let j = 0; j < cellsH; j++) {\n      for (let i = 0; i < cellsA; i++) {", "    if (pts) { for (const q of pts) addParticle(q[0], q[1], 0, 0); return; } // [closeout] custom seed\n    for (let j = 0; j < cellsH; j++) {\n      for (let i = 0; i < cellsA; i++) {");
  patch("  function substep(dt) {\n", "  const IG = { pacc: null, solves: 0, maxIt: 0, thinCalls: 0, thinFaces: 0, bothWet: 0, opp: 0, lastBoth: 0, lastOpp: 0, lastA: 0, lastB: 0, lastFaces: 0, lastWsZero: 0 };\n  function substep(dt) {\n");
  patch("    lastPcg.iters = it;\n", "    lastPcg.iters = it;\n    IG.solves++; if (it > IG.maxIt) IG.maxIt = it;\n");
  patch("    if (m === 0 || bNorm2 === 0 || maxR * resScale <= PCG_ABS_TOL) return;\n", "    if (m === 0 || bNorm2 === 0 || maxR * resScale <= PCG_ABS_TOL) { IG.solves++; return; }\n");
  patch("  function solvePressurePCG(dt) {\n", "  function solvePressurePCG(dt) {\n    solvePressurePCG_raw(dt);\n    if (!IG.pacc || IG.pacc.length !== pcgP.length) IG.pacc = new Float64Array(pcgP.length);\n    for (let c = 0; c < pcgP.length; c++) IG.pacc[c] += pcgP[c];\n  }\n  function solvePressurePCG_raw(dt) {\n");
  patch("      u.set(uSnap); v.set(vSnap); enforceBoundary(); solvePressurePCG(dt);", "      u.set(uSnap); v.set(vSnap); if (IG.pacc) IG.pacc.fill(0); enforceBoundary(); solvePressurePCG(dt);");
  // thin-ghost statistics (read-only): per fill call, faces with both copy sources carrying P2G mass
  patch("    const F = g.faces, S = g.src, TH = g.thin;\n", "    const F = g.faces, S = g.src, TH = g.thin;\n    let igB = 0, igO = 0, igA = 0, igBB = 0, igF = 0, igZ = 0;\n");
  patch("        if (ws > 0) { G[F[q]] = (wv / ws) * (nc / c); continue; }\n",
    "        { igF++; let nw = 0, a = 0, b = 0; for (let t = 4 * q; t < 4 * q + 4; t++) { const k = S[t]; if (k >= 0 && W[k] > 0) { nw++; if (nw === 1) a = G[k]; else b = G[k]; } } if (nw >= 2) { igB++; igA += Math.abs(a); igBB += Math.abs(b); if (a * b < 0) igO++; } if (!(ws > 0)) igZ++; }\n        if (ws > 0) { G[F[q]] = (wv / ws) * (nc / c); continue; }\n");
  // book at the end of fillGhosts (only calls that have thin faces)
  const fgEnd = "    fillGhosts(u, ghostU, uW);\n    fillGhosts(v, ghostV, vW);\n";
  patch(fgEnd, fgEnd.replace("fillGhosts(v, ghostV, vW);", "fillGhosts(v, ghostV, vW);"));
  patch("  function fillGhosts(G, g, W) {\n", "  function fillGhosts(G, g, W) { fillGhostsRaw(G, g, W); }\n  function fillGhostsRaw(G, g, W) {\n");
  src = src.replace(/(function fillGhostsRaw\(G, g, W\) \{[\s\S]*?\n    \}\n)(  \})/, (m0, body, close) => body + "    if (igF) { IG.thinCalls++; IG.thinFaces += igF; IG.bothWet += igB; IG.opp += igO; IG.lastFaces = igF; IG.lastBoth = igB; IG.lastOpp = igO; IG.lastA = igB ? igA / igB : 0; IG.lastB = igB ? igBB / igB : 0; IG.lastWsZero = igZ; }\n" + close);
  if (!/IG\.thinCalls\+\+/.test(src)) throw new Error("fillGhosts end anchor");
  patch("    setSubstepHook(fn) { substepHook = fn || null; },", "    setSubstepHook(fn) { substepHook = fn || null; },\n    get IG() { return IG; }, get wallRelease() { return wallRelease; },");
  const tmp = path.join(os.tmpdir(), `sim_item1co_${process.pid}.js`);
  fs.writeFileSync(tmp, src); const m = require(tmp); fs.unlinkSync(tmp); return m;
}

function tiltSeed(S, dRightCells) {
  const r = S.radius, pts = [], sides = [];
  const mk = (x0, x1, dCells) => { // side spans [x0, x1] (solid faces), mean depth dCells
    const W = x1 - x0, ncol = Math.round(W / H), npt = 2 * ncol, dx = (W - 2 * r) / (npt - 1), xc = 0.5 * (x0 + x1), d = dCells * H;
    let n = 0;
    for (let c = 0; c < npt; c++) {
      const x = x0 + r + c * dx, eta = d + SLOPE * (x - xc);
      for (let k = 0; ; k++) { const hgt = r + k * 0.5 * H; if (hgt >= eta) break; pts.push([x, FLOOR - hgt]); n++; }
    }
    sides.push({ x0, x1, width_m: W, columns: ncol, ptsPerRow: npt, dx_over_h: dx / H, d_target_m: d, n, d_fromCount_m: (n * dx * 0.5 * H) / W, amplitude_cells_at_ends: (SLOPE * (W / 2 - r)) / H });
  };
  mk(H, 64 * H, D_CELLS);
  mk(65 * H, (NX + 1) * H, dRightCells);
  return { pts, sides };
}

function runCase(name) {
  const S = load().FlipSolver;
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  const barrier = name === "tilt_b1" || name === "tilt_b1_eqN";
  if (barrier) S.addObstacle({ type: "rect", cx: 64.5 * H, cy: 0.5 * (NY + 2) * H, w: H, h: (NY + 2) * H, angle: 0 });
  let seed = null;
  if (name === "base_dam") S.reset("dam");
  else {
    let dR = D_CELLS;
    if (name === "tilt_b1_eqN") { // lower the right depth until the counts match
      const nL = tiltSeed(S, D_CELLS).sides[0].n; let lo = D_CELLS * 0.9, hi = D_CELLS;
      for (let it = 0; it < 60; it++) { const mid = 0.5 * (lo + hi); if (tiltSeed(S, mid).sides[1].n > nL) hi = mid; else lo = mid; }
      dR = lo; while (tiltSeed(S, dR).sides[1].n < nL) dR += 1e-6;
    }
    seed = tiltSeed(S, dR);
    S.resetColumn(0, 0, seed.pts);
  }
  const n0 = S.count, fy = (S.ny - 1) * H, nx = S.nx, ny = S.ny;
  const E = () => { let k = 0, q = 0, vm = 0; for (let i = 0; i < S.count; i++) { const v2 = S.vx[i] ** 2 + S.vy[i] ** 2; k += 0.5 * v2; q += G * (fy - S.py[i]); if (v2 > vm) vm = v2; } return { e: (k + q) / S.count, vmax: Math.sqrt(vm) }; };
  // columnDepth (validate.js method) per side at t=0
  const colDepth = (i0, i1) => { const top = new Float64Array(nx).fill(Infinity); for (let k = 0; k < S.count; k++) { const i = Math.floor(S.px[k] / H); if (S.py[k] < top[i]) top[i] = S.py[k]; } let s = 0, n = 0; for (let i = i0; i <= i1; i++) { if (top[i] < Infinity) s += FLOOR - top[i] + 0.25 * H; n++; } return s / n; };
  const side0 = new Int8Array(n0); for (let i = 0; i < n0; i++) side0[i] = S.px[i] < 64.5 * H ? -1 : 1;
  const cntL = side0.reduce((s, x) => s + (x < 0 ? 1 : 0), 0);
  const res = { case: name, n0, obstacleDeletedAtStart: S.stats().obstacleDeleted };
  if (seed) Object.assign(res, { seed: seed.sides, count_left: cntL, count_right: n0 - cntL, count_diff_pct: (100 * (n0 - cntL - cntL)) / cntL,
    d_columnDepth_left_m: colDepth(1, 63), d_columnDepth_right_m: colDepth(65, 128) });
  const E0 = E().e;
  let err = null, maxSpeed = 0, minSdf = Infinity, insideMax = 0, subs = 0;
  const crossed = new Uint8Array(n0);
  const force = [], thinSub = { substeps: 0, withBoth: 0, faceBoth: 0, faceOpp: 0, sumA: 0, sumB: 0, faces: 0, wsZero: 0, absA: [], absB: [] };
  S.setSubstepHook((d, dt) => {
    subs++;
    let ins = 0;
    for (let i = 0; i < S.count; i++) {
      const x = S.px[i];
      if (barrier && ((side0[i] < 0 && x > 65 * H) || (side0[i] > 0 && x < 64 * H))) crossed[i] = 1;
      const sp = Math.hypot(S.vx[i], S.vy[i]); if (sp > maxSpeed) maxSpeed = sp;
      if (barrier) { const sd = S.obstacleSdf(x, S.py[i]); if (sd < minSdf) minSdf = sd; if (sd < 0) ins++; }
    }
    if (ins > insideMax) insideMax = ins;
    if (name !== "base_dam") {
      const P = S.IG.pacc, T = S.cellType, WR = S.wallRelease, sc = (RHO * H / dt) * H;
      let FL = 0, FR = 0;
      if (P) for (let j = 1; j < ny - 1; j++) {
        const cl = 63 + j * nx, cr = 65 + j * nx;
        if (T[cl] === 1 && !(WR[cl] & 2) && T[cl + 1] !== 1) FL += P[cl] * sc;
        if (T[cr] === 1 && !(WR[cr] & 1) && T[cr - 1] !== 1) FR += P[cr] * sc;
      }
      force.push([+d.t.toFixed(6), dt, FL, FR, FL - FR]);
      if (barrier) {
        const IG = S.IG; thinSub.substeps++;
        thinSub.faces += IG.lastFaces; thinSub.wsZero += IG.lastWsZero;
        if (IG.lastBoth > 0) { thinSub.withBoth++; thinSub.faceBoth += IG.lastBoth; thinSub.faceOpp += IG.lastOpp; thinSub.absA.push(IG.lastA); thinSub.absB.push(IG.lastB); }
        IG.lastFaces = IG.lastBoth = IG.lastOpp = IG.lastWsZero = 0; IG.lastA = IG.lastB = 0;
      }
    }
  });
  const eSer = [[0, E0]];
  try {
    while (S.simTime < T_END - 1e-9) { S.step(); const e = E(); eSer.push([+S.simTime.toFixed(6), e.e]); if (e.vmax > maxSpeed) maxSpeed = e.vmax; }
  } catch (x) { err = String(x.message || x); }
  let leaks = 0; for (let i = 0; i < n0; i++) leaks += crossed[i];
  let finalL = 0; for (let i = 0; i < S.count; i++) if (S.px[i] < 64.5 * H) finalL++;
  const med = (a) => { if (!a.length) return null; const b = a.slice().sort((x, y) => x - y); return b[b.length >> 1]; };
  Object.assign(res, { nanGuard: err, simTime: S.simTime, particles: [n0, S.count], E0, E30: eSer[eSer.length - 1][1], E30_le_E0: eSer[eSer.length - 1][1] <= E0,
    capHits: S.capLog.length, pcgCapHits: S.pcgCapLog.length, pcgSolves: S.IG.solves, pcgMaxIters: S.IG.maxIt, maxParticleSpeed: maxSpeed, substeps: subs,
    leaks_everCrossed: barrier ? leaks : null, finalLeftCount: finalL, minSdfOverR: barrier ? minSdf / S.radius : null, insideMaxPerSubstep: barrier ? insideMax : null,
    obstacleDeletedEnd: S.stats().obstacleDeleted, pushStuck: S.stats().pushStuck });
  if (barrier) res.thinGhost = { substeps: thinSub.substeps, substepsWithBothSidesMass: thinSub.withBoth, thinFacesPerSubstep_lastFill: thinSub.faces / thinSub.substeps,
    bothMassFacesPerSubstep_mean: thinSub.withBoth ? thinSub.faceBoth / thinSub.withBoth : 0, oppositeSignFraction: thinSub.faceBoth ? thinSub.faceOpp / thinSub.faceBoth : null,
    median_abs_v_sourceA: med(thinSub.absA), median_abs_v_sourceB: med(thinSub.absB), facesWithNoMass_lastFill_perSubstep: thinSub.wsZero / thinSub.substeps,
    allCalls: { fillCallsWithThin: S.IG.thinCalls, thinFaceVisits: S.IG.thinFaces, bothMassVisits: S.IG.bothWet, oppositeSignVisits: S.IG.opp } };
  res.energySeries = eSer; res.forceSeries = force;
  return res;
}

if (CASE) { process.stdout.write(JSON.stringify(runCase(CASE))); process.exit(0); }

function child(name) {
  return new Promise((ok, no) => {
    const f = path.join(os.tmpdir(), "item1co_" + name + ".json"), fd = fs.openSync(f, "w");
    const p = spawn(process.execPath, ["--max-old-space-size=4096", __filename, "--case", name], { cwd: ROOT, stdio: ["ignore", fd, "inherit"] });
    p.on("close", (c) => { fs.closeSync(fd); c === 0 ? ok(JSON.parse(fs.readFileSync(f, "utf8"))) : no(new Error(name + " exit " + c)); });
  });
}
(async () => {
  const names = ["base_dam", "tilt_b1", "tilt_none", "tilt_b1_eqN"];
  const R = {}; (await Promise.all(names.map(child))).forEach((r, k) => (R[names[k]] = r));
  // baseline energy
  fs.writeFileSync(path.join(V, "item1_baseline_energy.csv"), "t,E_none\n" + R.base_dam.energySeries.map((r) => r.join(",")).join("\n") + "\n");
  // tilt series
  for (const n of ["tilt_b1", "tilt_none", "tilt_b1_eqN"]) {
    const suf = n === "tilt_b1" ? "" : n === "tilt_none" ? "_nobarrier" : "_eqN";
    if (n !== "tilt_none") fs.writeFileSync(path.join(V, `item1_api_tilt_force${suf}.csv`), "t,dt,F_left_N_per_m,F_right_N_per_m,F_net_N_per_m\n" + R[n].forceSeries.map((r) => r.join(",")).join("\n") + "\n");
    fs.writeFileSync(path.join(V, `item1_api_tilt_energy${suf}.csv`), "t,E\n" + R[n].energySeries.map((r) => r.join(",")).join("\n") + "\n");
  }
  const strip = (r) => { const c = Object.assign({}, r); delete c.energySeries; delete c.forceSeries; return c; };
  fs.writeFileSync(path.join(os.tmpdir(), "item1co_summary.json"), JSON.stringify(Object.fromEntries(names.map((n) => [n, strip(R[n])])), null, 1));
  console.log(JSON.stringify(Object.fromEntries(names.map((n) => [n, strip(R[n])])), null, 1));
})().catch((e) => { console.error(e); process.exit(1); });
