// Item 1, Phase C: frozen gates for one-cell-thick obstacles (Validation), on the CURRENT sim.js.
// Cases b1 (1-cell barrier, column 64) and b2 (2-cell, columns 64-65), full height, dam break,
// damping 0, speed clamp OFF, NaN guard, reviewed defaults, fresh Node process each, 30 s, CFL 1.
// Usage: node test/item1_gates.js [--sim path] [--tag NAME]  -> validation/item1_gates<tag>.json (+ _energy.csv)
// The solver is loaded from an OS-temp copy with: MAX_SPEED = Infinity, a NaN guard at the start of
// every substep, measurement of the obstacle push-out rows (PE / KE / normal-velocity removal), a
// PCG solve counter and the accumulated pressure of the final active-set pass (for the hydrostatic
// report). "b1_plain" (clamp-off copy without measurement) checks the measurement is read-only
// (identical final state hash).
const fs = require("fs"), os = require("os"), path = require("path"), crypto = require("crypto"), { spawn } = require("child_process");
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf("--" + k); return i >= 0 ? args[i + 1] : d; };
const ROOT = path.join(__dirname, ".."), SIM = path.resolve(opt("sim", path.join(ROOT, "sim.js"))), TAG = opt("tag", "");
const NX = 128, NY = 72, H = 3 / 128, G = 9.81, RHO = 1000, T_END = +opt("tend", "30"), FLOOR = (NY + 1) * H;
const BAR = { b1: [64, 65], b2: [64, 66] };
const CASES = { b1: ["b1", true], b2: ["b2", true], b1_plain: ["b1", false] };

function load(measure) {
  let src = fs.readFileSync(SIM, "utf8");
  const patch = (a, b) => { if (src.split(a).length !== 2) throw new Error("anchor: " + a.slice(0, 70)); src = src.replace(a, b); };
  patch("  const MAX_SPEED = 14;", "  const MAX_SPEED = Infinity; // [item1 gates] clamp OFF");
  patch("  function substep(dt) {\n    pinHeld();\n", "  function substep(dt) {\n    for (let i = 0; i < count; i++) if (!(Number.isFinite(px[i]) && Number.isFinite(py[i]) && Number.isFinite(pvx[i]) && Number.isFinite(pvy[i]))) throw new Error(\"NaN guard: particle \" + i + \" at t=\" + simTime);\n    pinHeld();\n");
  if (measure) {
    patch("  function substep(dt) {\n", `  const IG = { ser: [], vr: 0, solves: 0, maxIt: 0, pacc: null };
  function igE() { let k = 0, q = 0; const fy = (ny - 1) * h, g = params.gravity; for (let i = 0; i < count; i++) { k += 0.5 * (pvx[i] * pvx[i] + pvy[i] * pvy[i]); q += g * (fy - py[i]); } return [k / count, q / count]; }
  function substep(dt) {\n`);
    patch("    if (obsActive) pushOutObstacles();\n", "    if (obsActive) { const a = igE(); IG.vr = 0; pushOutObstacles(); const b = igE(); IG.ser.push([simTime + dt, b[1] - a[1], b[0] - a[0] - IG.vr, IG.vr]); }\n");
    patch("        if (vn < 0) { vx -= vn * o.nx; vy -= vn * o.ny; }\n", "        if (vn < 0) { const e0 = vx * vx + vy * vy; vx -= vn * o.nx; vy -= vn * o.ny; IG.vr += 0.5 * (vx * vx + vy * vy - e0) / count; }\n");
    patch("    lastPcg.iters = it;\n", "    lastPcg.iters = it;\n    IG.solves++; if (it > IG.maxIt) IG.maxIt = it;\n");
    // accumulated pressure of the final active-set pass (passes restart from the snapshot; put-back
    // and re-release solves after the loop are increments on top of the final pass)
    patch("    if (m === 0 || bNorm2 === 0 || maxR * resScale <= PCG_ABS_TOL) return;\n", "    if (m === 0 || bNorm2 === 0 || maxR * resScale <= PCG_ABS_TOL) { IG.solves++; return; }\n");
    patch("  function solvePressurePCG(dt) {\n", "  function solvePressurePCG(dt) {\n    solvePressurePCG_raw(dt);\n    if (!IG.pacc || IG.pacc.length !== pcgP.length) IG.pacc = new Float64Array(pcgP.length);\n    for (let c = 0; c < pcgP.length; c++) IG.pacc[c] += pcgP[c];\n  }\n  function solvePressurePCG_raw(dt) {\n");
    patch("      u.set(uSnap); v.set(vSnap); enforceBoundary(); solvePressurePCG(dt);", "      u.set(uSnap); v.set(vSnap); if (IG.pacc) IG.pacc.fill(0); enforceBoundary(); solvePressurePCG(dt);");
    patch("    setSubstepHook(fn) { substepHook = fn || null; },", "    setSubstepHook(fn) { substepHook = fn || null; },\n    get IG() { return IG; }, get wallRelease() { return wallRelease; },");
  }
  const tmp = path.join(os.tmpdir(), `sim_item1g_${process.pid}_${measure ? "m" : "p"}.js`);
  fs.writeFileSync(tmp, src); const m = require(tmp); fs.unlinkSync(tmp); return m;
}

function runCase(name) {
  const [kind, measure] = CASES[name];
  const S = load(measure).FlipSolver;
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  const [a, b] = BAR[kind];
  S.addObstacle({ type: "rect", cx: 0.5 * (a + b) * H, cy: 0.5 * (NY + 2) * H, w: (b - a) * H, h: (NY + 2) * H, angle: 0 });
  S.reset("dam");
  const fy = (S.ny - 1) * H, nx = S.nx, ny = S.ny;
  const energy = () => { let k = 0, q = 0, vm = 0; for (let i = 0; i < S.count; i++) { const v2 = S.vx[i] ** 2 + S.vy[i] ** 2; k += 0.5 * v2; q += G * (fy - S.py[i]); if (v2 > vm) vm = v2; } return { e: (k + q) / S.count, ke: k / S.count, vmax: Math.sqrt(vm) }; };
  const E0 = energy().e, n0 = S.count, xR = b * H;
  let err = null, maxSpeed = 0, maxDown = 0, everDown = new Uint8Array(n0), minSdf = Infinity, insideMax = 0, subs = 0;
  const forceS = [], depthS = []; // last 5 s
  const surf = () => { // validate.js columnDepth over upstream columns 1..a-1
    const T = S.cellType, topRow = new Int32Array(nx).fill(-1), topY = new Float64Array(nx).fill(Infinity);
    for (let i = 1; i < a; i++) { let j = ny - 2; if (T[i + j * nx] !== 1) continue; while (j > 1 && T[i + (j - 1) * nx] === 1) j--; topRow[i] = j; }
    for (let k = 0; k < S.count; k++) { const i = Math.floor(S.px[k] / H), j = Math.floor(S.py[k] / H); if (i < a && topRow[i] === j && S.py[k] < topY[i]) topY[i] = S.py[k]; }
    let s = 0, n = 0; for (let i = 1; i < a; i++) { s += topY[i] === Infinity ? 0 : FLOOR - topY[i] + 0.25 * H; n++; } return s / n;
  };
  S.setSubstepHook((d, dt) => {
    subs++;
    let down = 0, ins = 0;
    for (let i = 0; i < S.count; i++) {
      if (S.px[i] > xR) { down++; everDown[i] = 1; }
      const sd = S.obstacleSdf(S.px[i], S.py[i]); if (sd < minSdf) minSdf = sd; if (sd < 0) ins++;
      const sp = Math.hypot(S.vx[i], S.vy[i]); if (sp > maxSpeed) maxSpeed = sp;
    }
    if (down > maxDown) maxDown = down; if (ins > insideMax) insideMax = ins;
    if (measure && d.t >= T_END - 5 - 1e-9) {
      // wet face = left face of the barrier (x = a h), cells (a-1, j); constrained faces carry p, released faces p = 0
      const P = S.IG.pacc, T = S.cellType, WR = S.wallRelease;
      let F = 0, wet = 0;
      for (let j = 1; j < ny - 1; j++) { const c = (a - 1) + j * nx; if (T[c] !== 1 || (WR[c] & 2)) continue; F += (P[c] * RHO * H / dt) * H; wet++; }
      forceS.push([d.t, F, wet]);
    }
  });
  const eSer = [[0, E0]];
  let maxE = E0, maxEt = 0;
  try {
    while (S.simTime < T_END - 1e-9) {
      S.step();
      const e = energy(); eSer.push([+S.simTime.toFixed(6), e.e]);
      if (e.e > maxE) { maxE = e.e; maxEt = S.simTime; }
      if (e.vmax > maxSpeed) maxSpeed = e.vmax;
      if (measure && S.simTime >= T_END - 5 - 1e-9) depthS.push(surf());
    }
  } catch (x) { err = String(x.message || x); }
  let finalDown = 0; for (let i = 0; i < S.count; i++) if (S.px[i] > xR) finalDown++;
  let ever = 0; for (let i = 0; i < n0; i++) ever += everDown[i];
  const hsh = crypto.createHash("sha256"); for (const A of [S.px, S.py, S.vx, S.vy]) hsh.update(Buffer.from(A.buffer, A.byteOffset, S.count * 4));
  const res = { case: name, kind, barrier_cells: [a, b], sim: SIM, nanGuard: err, simTime: S.simTime, particles: [n0, S.count], E0, E30: eSer[eSer.length - 1][1], maxE, maxEt,
    capHits: S.capLog.length, pcgCapHits: S.pcgCapLog.length, maxParticleSpeed: maxSpeed, substeps: subs,
    leaks: { maxDownstream: maxDown, finalDownstream: finalDown, everDownstream: ever }, minSdfOverR: minSdf / S.radius, insideMaxPerSubstep: insideMax,
    finalStateHash: hsh.digest("hex").slice(0, 16) };
  if (measure) {
    const IG = S.IG, ser = IG.ser;
    let cpe = 0, cke = 0, cvr = 0; for (const r of ser) { cpe += r[1]; cke += r[2]; cvr += r[3]; }
    // closed 0.5 s windows [t - 0.5, t] over substep end times
    const pref = [0]; for (const r of ser) pref.push(pref[pref.length - 1] + r[1] + r[2] + r[3]);
    let i0 = 0, maxW = -Infinity, maxWt = null, maxVrSub = -Infinity;
    for (let j = 0; j < ser.length; j++) { while (ser[i0][0] < ser[j][0] - 0.5 - 1e-9) i0++; const w = pref[j + 1] - pref[i0]; if (w > maxW) { maxW = w; maxWt = [ser[i0][0], ser[j][0]]; } if (ser[j][3] > maxVrSub) maxVrSub = ser[j][3]; }
    const fMean = forceS.length ? forceS.reduce((s, q) => s + q[1], 0) / forceS.length : null;
    const d = depthS.length ? depthS.reduce((s, q) => s + q, 0) / depthS.length : null;
    const Fref = d != null ? 0.5 * RHO * G * d * d : null;
    Object.assign(res, { pcgSolves: IG.solves, pcgMaxIters: IG.maxIt,
      obstacleRows: { pushPE: cpe, pushKE: cke, velRemoval: cvr, pushPlusRemoval: cpe + cke + cvr, maxWindow05s: maxW, maxWindow05sAt: maxWt, maxRemovalPerSubstep: maxVrSub },
      hydrostatic: { window: [T_END - 5, T_END], samples: forceS.length, F_mean_N_per_m: fMean, F_min: forceS.length ? Math.min(...forceS.map((q) => q[1])) : null, F_max: forceS.length ? Math.max(...forceS.map((q) => q[1])) : null,
        wetCells_mean: forceS.length ? forceS.reduce((s, q) => s + q[2], 0) / forceS.length : null, d_upstream_m: d, F_ref_half_rho_g_d2: Fref, ratio: fMean != null && Fref ? fMean / Fref : null, band_2h_over_d: d ? (2 * H) / d : null },
      energySeries: eSer });
  }
  return res;
}

if (args[0] === "--case") { const r = runCase(args[1]); process.stdout.write(JSON.stringify(r), () => process.exit(0)); }
else {
  const self = path.resolve(__filename), pass = args.slice();
  const run = (c) => new Promise((res, rej) => { const p = spawn(process.execPath, [self, "--case", c, ...pass], { stdio: ["ignore", "pipe", "inherit"] }); let o = ""; p.stdout.on("data", (x) => (o += x)); p.on("close", (code) => (code ? rej(new Error(c + " exited " + code)) : res(JSON.parse(o)))); });
  Promise.all(Object.keys(CASES).map(run)).then((R) => {
    const by = Object.fromEntries(R.map((r) => [r.case, r]));
    const b1 = by.b1, b2 = by.b2;
    const rise = (r) => (r.maxE - r.E0) / r.E0;
    const g = {};
    for (const r of [b1, b2]) {
      const o = r.obstacleRows;
      g[r.case] = {
        noNaN: { value: r.nanGuard, pass: r.nanGuard === null },
        substepCapHits0: { value: r.capHits, pass: r.capHits === 0 },
        maxSpeedBelow14: { value: r.maxParticleSpeed, pass: r.maxParticleSpeed < 14 },
        pcgCapHits0: { value: r.pcgCapHits, solves: r.pcgSolves, maxIters: r.pcgMaxIters, pass: r.pcgCapHits === 0 },
        E30_below_E0: { value: [r.E30, r.E0], pass: r.E30 < r.E0 },
        obstacleCumulative_le_1e4: { value: o.pushPlusRemoval, pass: o.pushPlusRemoval <= 1e-4 },
        obstacleEvery05sWindow_le_1e4: { value: o.maxWindow05s, at: o.maxWindow05sAt, pass: o.maxWindow05s <= 1e-4 },
        removal_le_0: { value: o.velRemoval, maxPerSubstep: o.maxRemovalPerSubstep, pass: o.velRemoval <= 0 && o.maxRemovalPerSubstep <= 0 },
      };
    }
    g.b1.leaks0_minSdf_ge0 = { value: { leaks: b1.leaks, minSdfOverR: b1.minSdfOverR }, pass: b1.leaks.maxDownstream === 0 && b1.leaks.everDownstream === 0 && b1.minSdfOverR >= 0 };
    g.energyRise_b1_vs_b2 = { b1_rise_frac: rise(b1), b2_rise_frac: rise(b2), b1_at: b1.maxEt, b2_at: b2.maxEt, excess_frac: rise(b1) - rise(b2), limit: 0.002, pass: rise(b1) - rise(b2) <= 0.002 };
    g.instrumentationReadOnly = { b1: b1.finalStateHash, b1_plain: by.b1_plain.finalStateHash, identical: b1.finalStateHash === by.b1_plain.finalStateHash };
    const strip = (r) => { const c = Object.assign({}, r); delete c.energySeries; return c; };
    const out = { settings: { grid: "128x72", h: H, damping: 0, cfl: 1, speedClamp: "off", nanGuard: true, tEnd: T_END, sim: SIM }, gates: g, b1: strip(b1), b2: strip(b2), b1_plain: strip(by.b1_plain) };
    fs.writeFileSync(path.join(ROOT, "validation", `item1_gates${TAG}.json`), JSON.stringify(out, null, 1));
    const lines = ["t,E_b1,E_b2"]; for (let i = 0; i < Math.max(b1.energySeries.length, b2.energySeries.length); i++) { const p = b1.energySeries[i], q = b2.energySeries[i]; lines.push([(p || q)[0], p ? p[1] : "", q ? q[1] : ""].join(",")); }
    fs.writeFileSync(path.join(ROOT, "validation", `item1_gates_energy${TAG}.csv`), lines.join("\n") + "\n");
    console.log(JSON.stringify(g, null, 1));
    console.log("hydrostatic b1", JSON.stringify(b1.hydrostatic), "\nhydrostatic b2", JSON.stringify(b2.hydrostatic));
  }).catch((e) => { console.error(e); process.exit(1); });
}
