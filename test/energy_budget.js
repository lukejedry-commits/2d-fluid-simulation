// Stage-by-stage energy budget over the pre-impact rise window (diagnostics only).
// Usage: node test/energy_budget.js [--cells 48|24] [--sep neutral|legacy] [--boundary mirror|legacy] [--sim sim.js]
// Writes an instrumented copy of sim.js to the OS temp dir (the repo file is not
// modified), runs the exact validation column once to find the unsmoothed pre-impact
// rise window (same metric as validate.js), then re-runs it and integrates the change of
// particle energy (KE + PE, J/kg) or mass-weighted grid KE across every substep stage.
// Grid KE = 0.5 * sum(face weight * velocity^2) / N, with the weights from P2G.
const fs = require("fs"), os = require("os"), path = require("path");
const EM = require("./energy_metrics.js");
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf("--" + k); return i >= 0 ? args[i + 1] : d; };
const cellsH = +opt("cells", 48), SEP = opt("sep", "neutral"), BND = opt("boundary", "mirror");
const SIM = path.resolve(opt("sim", path.join(__dirname, "..", "sim.js")));
let src = fs.readFileSync(SIM, "utf8");
function patch(old, neu) { if (src.split(old).length !== 2) throw new Error("energy_budget: pattern not found once: " + old.slice(0, 60)); src = src.replace(old, neu); }
patch(`    if (mirror) enforceBoundary();\n`, `    if (mirror) { const k0 = bOn ? kg() : 0; enforceBoundary(); if (bOn) { bP2G = kg() - k0; add("enforceBoundary (inside P2G, before FLIP reference)", bP2G); } }\n`);
patch(`  function solvePressure(dt) {\n    buildFluidList();\n`, `  let driftOffOnce = false;\n  function solvePressure(dt) {\n    buildFluidList();\n    if (driftOffOnce) cellDrift.fill(0);\n`);
const sub = src.match(/  function substep\(dt\) \{\n[\s\S]*?\n  \}\n/);
if (!sub) throw new Error("energy_budget: substep() not found");
patch(sub[0], `  const B = {}; let bOn = false; let bP2G = 0;
  const add = (k, x) => { B[k] = (B[k] || 0) + x; };
  function kp() { let k = 0, q = 0; const fy = (ny - 1) * h, g = params.gravity; for (let i = 0; i < count; i++) { k += 0.5 * (pvx[i] * pvx[i] + pvy[i] * pvy[i]); q += g * (fy - py[i]); } return [k / count, q / count]; }
  function kg() { let k = 0; for (let i = 0; i < u.length; i++) k += uW[i] * u[i] * u[i]; for (let i = 0; i < v.length; i++) k += vW[i] * v[i] * v[i]; return (0.5 * k) / count; }
` + sub[0].replace("  function substep(dt) {", "  function substepPlain(dt) {") + `  function substep(dt) {
    if (!bOn) return substepPlain(dt);
    pinHeld();
    let a = kp();
    advectParticles(dt); let b = kp(); add("advection (particle PE)", b[1] - a[1]);
    separateParticles(SEPARATION_ITERS); a = kp(); add("separation (particles)", a[0] + a[1] - b[0] - b[1]);
    collideWalls(); b = kp(); add("wall clamp (particles)", b[0] + b[1] - a[0] - a[1]);
    let g0 = 0; { const k0 = b[0]; bP2G = 0; particlesToGrid(); markCells(); g0 = kg(); add("transfer P2G (incl. ghost fold)", g0 - k0 - bP2G); }
    applyGravity(dt); let g = kg(); add("gravity (grid)", g - g0);
    applyViscosity(); let g2 = kg(); add("viscosity", g2 - g);
    enforceBoundary(); g = kg(); add("enforceBoundary (pre-projection)", g - g2);
    const us = Float32Array.from(u), vs = Float32Array.from(v);
    driftOffOnce = true; solvePressure(dt); driftOffOnce = false; const gnd = kg();
    u.set(us); v.set(vs);
    solvePressure(dt); g2 = kg(); add("pressure: projection (no drift)", gnd - g); add("pressure: drift target", g2 - gnd);
    extrapolateVelocities(); g = kg(); add("extrapolation", g - g2);
    enforceBoundary(); g2 = kg(); add("enforceBoundary (post-projection)", g2 - g);
    if (substepHook) substepHook(divergenceStats(dt), dt);
    gridToParticles(); a = kp(); add("transfer G2P (incl. FLIP)", a[0] - g2);
  }
`);
patch(`    setSubstepHook(fn) { substepHook = fn || null; },`, `    setSubstepHook(fn) { substepHook = fn || null; },\n    budget: B, setBudget(on) { bOn = on; if (on) for (const k in B) delete B[k]; },`);
const tmp = path.join(os.tmpdir(), `sim_budget_${process.pid}.js`);
fs.writeFileSync(tmp, src);
function mk() { delete require.cache[tmp]; const { FlipSolver: S } = require(tmp); S.init(128, 72, 3 / 128); S.params.damping = 0; S.params.separation = SEP; S.params.boundary = BND; S.resetColumn(24, cellsH); return S; }
let S = mk();
const h = S.cellSize, fy = (S.ny - 1) * h, xMax = (S.nx - 2) * h - S.radius;
const E = () => { let k = 0, q = 0; for (let i = 0; i < S.count; i++) { k += 0.5 * (S.vx[i] ** 2 + S.vy[i] ** 2); q += 9.81 * (fy - S.py[i]); } return (k + q) / S.count; };
const t = [0], e = [E()]; let reach = null;
S.setSubstepHook((d) => { if (reach === null) { let m = -1; for (let i = 0; i < S.count; i++) if (fy - S.py[i] <= h) m = Math.max(m, S.px[i] - h); if (m >= xMax - 1e-6) reach = d.t; } });
while (reach === null || S.simTime < reach + 0.05) { S.step(); t.push(S.simTime); e.push(E()); if (S.simTime > 2) break; }
// same window rule as validate.js: frame samples, t in [0.2 s, frame time of wall arrival]
const reachFrame = Math.ceil(reach * 60 - 1e-6) / 60;
const pre = EM.riseIn(t, e, 0.2, reachFrame);
S = mk();
while (S.simTime < pre.fromT - 1e-9) S.step();
const Ea = E(); S.setBudget(true);
while (S.simTime < pre.toT - 1e-9) S.step();
const Eb = E(); S.setBudget(false);
fs.unlinkSync(tmp);
let sum = 0; const rows = {};
for (const [k, v] of Object.entries(S.budget)) { rows[k] = +v.toFixed(5); sum += v; }
const out = { column: `H/a=${cellsH / 24}`, separation: SEP, boundary: BND, preImpactRisePctE0: +pre.risePctE0.toFixed(3), window_s: [+pre.fromT.toFixed(3), +pre.toT.toFixed(3)], E0: +e[0].toFixed(5), dE_J_per_kg: +(Eb - Ea).toFixed(5), sumOfRows: +sum.toFixed(5), rows_J_per_kg: rows };
console.log(JSON.stringify(out));
