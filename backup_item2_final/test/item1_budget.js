// Item 1, Phase A/B: per-substep, per-stage energy budget for full-height barriers
// (b1 = 1 cell at column 64, b2 = 2 cells at columns 64-65), dam break, damping 0,
// speed clamp OFF, NaN guard ON, reviewed defaults, fresh Node process per case.
//
// Usage: node test/item1_budget.js [--sim path/to/sim.js] [--tag NAME] [--tend 1]
//   writes validation/item1_stage_budget<tag>.json and .csv
//
// sim.js is never modified. Each case loads an instrumented copy written to the OS temp dir:
//   MAX_SPEED -> Infinity (clamp off), CFL patched for the 0.5 cases, and measurement hooks.
// Energy per unit mass, E = mean over particles of (1/2 |v|^2 + g * height above the floor face).
// Every stage is booked as the change of a "state energy" between two consecutive marks, so the
// rows telescope exactly to the measured particle dE of the substep. Two definitions of the
// energy while velocities live on the grid:
//   A ("probe", primary): KE of the particles if G2P ran right now (read-only G2P sample with the
//     solver's FLIP/APIC blend; before uOld is saved the FLIP reference is the grid itself).
//     This sees every face the particles actually read, incl. massless ghost faces.
//   B ("faceMass"): grid KE = 1/2 sum_f m_f u_f^2 / N with m_f = P2G weights (uW, vW). Blind to
//     massless faces (ghosts), so boundary-fill effects show up in its G2P row instead.
// Particle stages (advection, separation, wall clamp, obstacle push-out) are particle energies
// in both methods. P2G = (grid state energy right after P2G divide) - (particle energy before);
// G2P = (particle energy after G2P) - (grid state energy before G2P).
// The active-set restarts inside solvePressureSeparating restore the snapshot, so bookings made
// by discarded passes are rolled back: only the pass that produced the final field counts.
const fs = require("fs"), os = require("os"), path = require("path"), crypto = require("crypto"), { spawn } = require("child_process");
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf("--" + k); return i >= 0 ? args[i + 1] : d; };
const ROOT = path.join(__dirname, "..");
const SIM = path.resolve(opt("sim", path.join(ROOT, "sim.js")));
const TAG = opt("tag", "");
const T_END = +opt("tend", "1");
const NX = 128, NY = 72, H = 3 / 128, G = 9.81;
const CASES = { b1_cfl1: ["b1", 1], b1_cfl05: ["b1", 0.5], b2_cfl1: ["b2", 1], b2_cfl05: ["b2", 0.5] };
const BAR = { b1: [64, 65], b2: [64, 66] };

function instrumented(cfl) {
  let src = fs.readFileSync(SIM, "utf8");
  const patch = (a, b) => { if (src.split(a).length !== 2) throw new Error("anchor not unique/found: " + a.slice(0, 80)); src = src.replace(a, b); };
  patch("  const MAX_SPEED = 14;", "  const MAX_SPEED = Infinity; // [item1 harness] clamp OFF");
  if (cfl !== 1) patch("  const CFL = 1.0;", "  const CFL = " + cfl + ";");
  // --- instrumentation core (inserted before substep) ---
  const core = `
  // ===== [item1 harness] stage budget instrumentation =====
  let IB_SUB = null, IB_E = [0, 0], IB_PLACE = "", IB_SELF = false, IB_CK = null;
  const IB_PCG = [];
  function ibPart() { let k = 0, q = 0; const fy = (ny - 1) * h, g = params.gravity; for (let i = 0; i < count; i++) { k += 0.5 * (pvx[i] * pvx[i] + pvy[i] * pvy[i]); q += g * (fy - py[i]); } return [k / count, q / count]; }
  function ibProbeC(G, GO, stride, offX, offY, VEL, CA, CB, out) {
    const ih = invH, flip = params.flipRatio, pic = 1 - flip;
    for (let p = 0; p < count; p++) {
      const fx = px[p] * ih - offX, fy = py[p] * ih - offY, bx = Math.floor(fx - 0.5), by = Math.floor(fy - 0.5), dx = fx - bx, dy = fy - by;
      const wx0 = 0.5 * (1.5 - dx) * (1.5 - dx), wx1 = 0.75 - (dx - 1) * (dx - 1), wx2 = 0.5 * (dx - 0.5) * (dx - 0.5);
      const wy0 = 0.5 * (1.5 - dy) * (1.5 - dy), wy1 = 0.75 - (dy - 1) * (dy - 1), wy2 = 0.5 * (dy - 0.5) * (dy - 0.5);
      let vel = 0, dv = 0, k = by * stride + bx;
      for (let b = 0; b < 3; b++, k += stride) {
        const wy = b === 0 ? wy0 : b === 1 ? wy1 : wy2, w0 = wx0 * wy, w1 = wx1 * wy, w2 = wx2 * wy;
        vel += w0 * G[k] + w1 * G[k + 1] + w2 * G[k + 2];
        if (GO) dv += w0 * (G[k] - GO[k]) + w1 * (G[k + 1] - GO[k + 1]) + w2 * (G[k + 2] - GO[k + 2]);
      }
      const vn = pic * vel + flip * (VEL[p] + dv);
      out[p] = vn;
    }
  }
  let IB_TX = null, IB_TY = null;
  function ibProbe() {
    if (!IB_TX || IB_TX.length < count) { IB_TX = new Float64Array(maxN); IB_TY = new Float64Array(maxN); }
    // self-referenced FLIP (G == GOLD) while P2G is not finished; afterwards GOLD = uOld/vOld
    ibProbeC(u, IB_SELF ? null : uOld, nx + 1, 0, 0.5, pvx, c00, c01, IB_TX);
    ibProbeC(v, IB_SELF ? null : vOld, nx, 0.5, 0, pvy, c10, c11, IB_TY);
    let k = 0; for (let p = 0; p < count; p++) k += 0.5 * (IB_TX[p] * IB_TX[p] + IB_TY[p] * IB_TY[p]);
    return k / count;
  }
  function ibGridKE() { let k = 0; for (let i = 0; i < u.length; i++) k += uW[i] * u[i] * u[i]; for (let i = 0; i < v.length; i++) k += vW[i] * v[i] * v[i]; return (0.5 * k) / count; }
  let IB_PE = 0;
  function ibBook(label, a, b) { if (!IB_SUB) return; const r = IB_SUB[label] || (IB_SUB[label] = [0, 0]); r[0] += a; r[1] += b; }
  function IBP(label) { if (!IB_SUB) return; const e = ibPart(), s = e[0] + e[1]; ibBook(label, s - IB_E[0], s - IB_E[1]); IB_E = [s, s]; IB_PE = e[1]; }
  function IBG(label) { if (!IB_SUB) return; const a = ibProbe() + IB_PE, b = ibGridKE() + IB_PE; ibBook(label, a - IB_E[0], b - IB_E[1]); IB_E = [a, b]; }
  function IBGP(label) { if (!IB_SUB) return; IBG((IB_PLACE ? IB_PLACE + "." : "") + label); }
  function ibCkSave() { IBGP("other"); IB_CK = { sub: JSON.stringify(IB_SUB), e: IB_E.slice() }; }
  function ibCkRestore() { if (!IB_SUB || !IB_CK) return; IB_SUB = JSON.parse(IB_CK.sub); IB_E = IB_CK.e.slice(); }
  // ===== end instrumentation core =====
`;
  patch("  // ---- time stepping ------------------------------------------------------\n", core + "  // ---- time stepping ------------------------------------------------------\n");
  // P2G: book the transfer before the in-P2G boundary pass; FLIP reference = grid until uOld saved
  patch("    if (mirror) enforceBoundary();\n    uOld.set(u);\n    vOld.set(v);\n",
    "    IB_SELF = true; IBG(\"P2G\"); IB_PLACE = \"p2g\";\n    if (mirror) enforceBoundary();\n    uOld.set(u);\n    vOld.set(v);\n    IB_SELF = false;\n");
  // enforceBoundary split by code path
  patch("  function enforceBoundary() {\n    enforceWalls();\n    if (obsActive) enforceObstacles(); // obstacle faces follow the same rule as wall faces\n  }",
    "  function enforceBoundary() {\n    IBGP(\"pending\");\n    enforceWalls(); IBGP(\"enforce.walls\");\n    if (obsActive) enforceObstacles(); // obstacle faces follow the same rule as wall faces\n  }");
  { // pre-fix build: fillGhosts(u, ghostU) / fixed build (item 1): fillGhosts(u, ghostU, uW)
    const pre = "    fillGhosts(u, ghostU);\n    fillGhosts(v, ghostV);\n  }", fixed = "    fillGhosts(u, ghostU, uW);\n    fillGhosts(v, ghostV, vW);\n  }";
    const a = src.includes(fixed) ? fixed : pre;
    patch(a, "    IBGP(\"enforce.obstacleInterface\");\n" + a.replace(/\n  \}$/, "") + "\n    IBGP(\"enforce.obstacleGhostFill\");\n  }");
  }
  // active set: checkpoint at entry, roll back on every snapshot restore
  patch("    uSnap.set(u);\n    vSnap.set(v);\n", "    uSnap.set(u);\n    vSnap.set(v);\n    ibCkSave();\n");
  patch("      u.set(uSnap); v.set(vSnap); enforceBoundary(); solvePressurePCG(dt);",
    "      u.set(uSnap); v.set(vSnap); ibCkRestore(); enforceBoundary(); solvePressurePCG(dt);");
  // PCG wrapper: everything since the last mark inside the active set = put-back etc.
  patch("  function solvePressurePCG(dt) {\n", `  function solvePressurePCG(dt) {
    IBGP("other");
    solvePressurePCG_raw(dt);
    IBG("projection (PCG solve + pressure apply, incl. drift rhs)");
    if (IB_SUB) IB_PCG.push([+simTime.toFixed(6), lastPcg.iters, lastPcg.capHit ? 1 : 0, lastPcg.relRes, lastPcg.maxResDt, pcgCount]);
  }
  function solvePressurePCG_raw(dt) {
`);
  // substep with particle-stage marks + NaN guard
  const sub = src.match(/  function substep\(dt\) \{\n[\s\S]*?\n  \}\n/);
  const expected = ["pinHeld();", "advectParticles(dt);", "separateParticles(SEPARATION_ITERS);", "collideWalls();", "if (obsActive) pushOutObstacles();",
    "particlesToGrid();", "markCells();", "applyGravity(dt);", "applyViscosity();", "enforceBoundary();", "solvePressure(dt);",
    "extrapolateVelocities();", "enforceBoundary();", "if (substepHook) substepHook(divergenceStats(dt), dt);", "gridToParticles();"];
  const body = sub[0].split("\n").slice(1, -2).map((l) => l.trim());
  if (JSON.stringify(body) !== JSON.stringify(expected)) throw new Error("substep() changed:\n" + body.join("\n"));
  patch(sub[0], `  const IB_SERIES = []; let IB_ON = true;
  let IB_VR = 0;
  function substep(dt) {
    for (let i = 0; i < count; i++) if (!(Number.isFinite(px[i]) && Number.isFinite(py[i]) && Number.isFinite(pvx[i]) && Number.isFinite(pvy[i]))) throw new Error("NaN guard: particle " + i + " non-finite at t=" + simTime);
    pinHeld();
    IB_SUB = IB_ON ? {} : null; IB_PLACE = "";
    if (IB_SUB) { const e = ibPart(); IB_E = [e[0] + e[1], e[0] + e[1]]; IB_PE = e[1]; }
    const E0 = IB_SUB ? IB_E[0] : 0;
    advectParticles(dt); IBP("gravity+advection");
    separateParticles(SEPARATION_ITERS); IBP("separation (particles)");
    collideWalls(); IBP("wall clamp (particles)");
    if (obsActive) {
      const k0 = IB_SUB ? ibPart() : null; IB_VR = 0;
      pushOutObstacles();
      if (IB_SUB) { const k1 = ibPart(); ibBook("push-out PE", k1[1] - k0[1], k1[1] - k0[1]); ibBook("push-out KE (excl. removal)", k1[0] - k0[0] - IB_VR, k1[0] - k0[0] - IB_VR); ibBook("velocity removal", IB_VR, IB_VR); IB_E = [k1[0] + k1[1], k1[0] + k1[1]]; IB_PE = k1[1]; }
    }
    particlesToGrid();
    markCells();
    IB_PLACE = ""; applyGravity(dt); IBG("gravity+advection");
    applyViscosity(); IBG("viscosity");
    IB_PLACE = "pre"; enforceBoundary();
    IB_PLACE = "activeSet"; solvePressure(dt); IBGP("other");
    IB_PLACE = ""; extrapolateVelocities(); IBG("extrapolation");
    IB_PLACE = "post"; enforceBoundary();
    IB_PLACE = "";
    if (substepHook) substepHook(divergenceStats(dt), dt);
    gridToParticles(); IBP("G2P");
    if (IB_SUB) { const e = ibPart(); IB_SERIES.push({ t0: simTime, dt, dE: e[0] + e[1] - E0, rows: IB_SUB }); }
  }
`);
  patch("        if (vn < 0) { vx -= vn * o.nx; vy -= vn * o.ny; }\n",
    "        if (vn < 0) { const e0 = vx * vx + vy * vy; vx -= vn * o.nx; vy -= vn * o.ny; IB_VR += 0.5 * (vx * vx + vy * vy - e0) / count; }\n");
  patch("    setSubstepHook(fn) { substepHook = fn || null; },", "    setSubstepHook(fn) { substepHook = fn || null; },\n    get ibSeries() { return IB_SERIES; }, get ibPcg() { return IB_PCG; }, ibEnable(on) { IB_ON = !!on; },");
  const tmp = path.join(os.tmpdir(), `sim_item1_${process.pid}_${cfl}.js`);
  fs.writeFileSync(tmp, src);
  const m = require(tmp); fs.unlinkSync(tmp); return m;
}

function runCase(name) {
  const [kind, cfl] = CASES[name];
  const mod = instrumented(cfl), S = mod.FlipSolver;
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  const [a, b] = BAR[kind];
  S.addObstacle({ type: "rect", cx: 0.5 * (a + b) * H, cy: 0.5 * (NY + 2) * H, w: (b - a) * H, h: (NY + 2) * H, angle: 0 });
  S.reset("dam");
  const E0 = (() => { const fy = (S.ny - 1) * H; let k = 0, q = 0; for (let i = 0; i < S.count; i++) { k += 0.5 * (S.vx[i] ** 2 + S.vy[i] ** 2); q += G * (fy - S.py[i]); } return (k + q) / S.count; })();
  const xL = a * H;
  let contact = null, err = null, maxSpeed = 0;
  const sub = []; // per substep extra: [t_end, minDist/h, vmax]
  S.setSubstepHook((d) => {
    let md = Infinity, vm = 0;
    for (let i = 0; i < S.count; i++) { const dd = S.obstacleSdf(S.px[i], S.py[i]); if (dd < md) md = dd; const s = Math.hypot(S.vx[i], S.vy[i]); if (s > vm) vm = s; }
    if (contact === null && md < 2 * H) contact = d.t - 0; // substep end time of the first substep with a particle within 2h
    sub.push([d.t, md / H, vm]);
    if (vm > maxSpeed) maxSpeed = vm;
  });
  try { while (S.simTime < T_END - 1e-9) S.step(); } catch (e) { err = String(e.message || e); }
  const ser = S.ibSeries;
  // contact = first substep (by its start time) where a particle came within 2h (pre-G2P positions)
  const cIdx = sub.findIndex((r) => r[1] < 2);
  const tContact = cIdx >= 0 ? ser[cIdx].t0 : null;
  const labels = new Set(); for (const s of ser) for (const k in s.rows) labels.add(k);
  return { case: name, kind, cfl, barrier_cells: [a, b], err, E0, tContact, contactSubstep: cIdx, simTime: S.simTime, capHits: S.capLog.length, maxParticleSpeed: maxSpeed,
    labels: [...labels], series: ser.map((s, i) => ({ t0: s.t0, dt: s.dt, dE: s.dE, minDist_h: sub[i] ? sub[i][1] : null, vmax: sub[i] ? sub[i][2] : null, rows: s.rows })),
    pcg: S.ibPcg, pcgCapLog: S.pcgCapLog.length };
}

if (args[0] === "--case") {
  const r = runCase(args[1]);
  process.stdout.write(JSON.stringify(r), () => process.exit(0));
} else {
  const self = path.resolve(__filename);
  const passArgs = args.slice();
  const CACHE = path.join(os.tmpdir(), "item1_cache"); fs.mkdirSync(CACHE, { recursive: true });
  const run = (c) => new Promise((res, rej) => {
    const cf = path.join(CACHE, (TAG || "_pre") + "_" + T_END + "_" + c + ".json");
    if (args.includes("--resume") && fs.existsSync(cf)) return res(JSON.parse(fs.readFileSync(cf, "utf8")));
    const p = spawn(process.execPath, [self, "--case", c, ...passArgs], { stdio: ["ignore", "pipe", "inherit"] });
    let o = ""; p.stdout.on("data", (d) => (o += d)); p.on("close", (code) => { if (code) return rej(new Error(c + " exited " + code)); fs.writeFileSync(cf, o); res(JSON.parse(o)); });
  });
  Promise.all(Object.keys(CASES).map(run)).then((R) => {
    const by = Object.fromEntries(R.map((r) => [r.case, r]));
    // aggregate rows to reported stages (method A and B), both windows
    const stageOf = (label) => {
      if (label === "gravity+advection") return "gravity+advection";
      if (label === "P2G") return "P2G";
      if (label === "G2P") return "G2P";
      if (label.startsWith("projection")) return "projection (PCG solve + apply, incl. drift rhs)";
      if (label === "extrapolation") return "extrapolation";
      if (label === "viscosity") return "viscosity";
      if (label.startsWith("push-out PE")) return "push-out PE";
      if (label.startsWith("push-out KE")) return "push-out KE (excl. removal)";
      if (label === "velocity removal") return "velocity removal";
      if (label.startsWith("separation")) return "separation (particles)";
      if (label.startsWith("wall clamp")) return "wall clamp (particles)";
      if (label.startsWith("activeSet.enforce") || label === "activeSet.other" || label === "activeSet.pending") return "separating active set (restore/enforce w. release bits/put-back)";
      if (label.endsWith("enforce.walls")) return "boundary: box-wall faces (enforceWalls)";
      if (label.endsWith("enforce.obstacleInterface")) return "boundary: obstacle interface normal faces";
      if (label.endsWith("enforce.obstacleGhostFill")) return "boundary: obstacle ghost-face fill";
      if (label.endsWith(".pending") || label.endsWith(".other")) return "other (between marks)";
      return "other: " + label;
    };
    const STAGES = ["gravity+advection", "P2G", "boundary: box-wall faces (enforceWalls)", "boundary: obstacle interface normal faces", "boundary: obstacle ghost-face fill",
      "viscosity", "separating active set (restore/enforce w. release bits/put-back)", "projection (PCG solve + apply, incl. drift rhs)", "extrapolation", "G2P",
      "push-out PE", "push-out KE (excl. removal)", "velocity removal", "separation (particles)", "wall clamp (particles)", "other (between marks)"];
    function agg(r, t0, t1) {
      const s = r.series.filter((q) => q.t0 >= t0 - 1e-12 && q.t0 < t1 - 1e-12);
      const out = {}, place = {};
      for (const st of STAGES) out[st] = [0, 0];
      let dE = 0, dur = 0;
      for (const q of s) {
        dE += q.dE; dur += q.dt;
        for (const [k, v] of Object.entries(q.rows)) { const st = stageOf(k); (out[st] || (out[st] = [0, 0])); out[st][0] += v[0]; out[st][1] += v[1]; place[k] = (place[k] || 0) + v[0]; }
      }
      const n = s.length, rows = {};
      let sumA = 0, sumB = 0;
      for (const [st, v] of Object.entries(out)) { sumA += v[0]; sumB += v[1]; rows[st] = { A_total: v[0], A_perSubstep: n ? v[0] / n : 0, A_rate: dur ? v[0] / dur : 0, B_total: v[1], B_perSubstep: n ? v[1] / n : 0, B_rate: dur ? v[1] / dur : 0 }; }
      return { window: [t0, t1], substeps: n, duration: dur, dE_measured: dE, sumA, sumB, closureA: sumA - dE, closureB: sumB - dE, rows, byPlaceA: place };
    }
    const out = { method: "see header of test/item1_budget.js; A = read-only G2P probe (primary), B = face-mass grid KE", sim: SIM, tEnd: T_END, cases: {} };
    for (const r of R) {
      const tc = r.tContact;
      const pc = r.pcg;
      const inW = (t, a, b) => t >= a && t <= b;
      const pcgSum = (a, b) => { const s = pc.filter((q) => inW(q[0], a, b)); return { solves: s.length, capHits: s.filter((q) => q[2]).length, maxIters: Math.max(0, ...s.map((q) => q[1])), meanIters: s.length ? s.reduce((x, q) => x + q[1], 0) / s.length : 0, maxRelRes: Math.max(0, ...s.map((q) => q[3])), maxMaxResDt: Math.max(0, ...s.map((q) => q[4])) }; };
      out.cases[r.case] = { kind: r.kind, cfl: r.cfl, barrier_cells: r.barrier_cells, nanGuard: r.err, simTime: r.simTime, capHits: r.capHits, maxParticleSpeed: r.maxParticleSpeed,
        tContact: tc, contactWindow: tc != null ? agg(r, tc, tc + 0.2) : null, first1s: agg(r, 0, 1),
        // clean windows for the unstable case (requested windows are dominated by the last substeps before NaN):
        //   onset  = [contact, first substep start where E > 1.01 E0)   (b2/fixed: E never exceeds -> contact + 0.2 s)
        //   growth = [t(E > 1.01 E0), t(E > 10 E0))
        ...(() => {
          let E = r.E0, t101 = null, t10 = null;
          for (const q of r.series) { if (q.t0 >= tc && t101 === null && E > 1.01 * r.E0) t101 = q.t0; if (q.t0 >= tc && t10 === null && E > 10 * r.E0) t10 = q.t0; E += q.dE; }
          return { E0: r.E0, t_E_over_1p01E0: t101, t_E_over_10E0: t10, onsetWindow: tc != null ? agg(r, tc, t101 !== null ? t101 : tc + 0.2) : null, growthWindow: t101 !== null && t10 !== null ? agg(r, t101, t10) : null };
        })(), pcgContact: tc != null ? pcgSum(tc, tc + 0.2) : null, pcgFirst1s: pcgSum(0, 1), pcgCapLogTotal: r.pcgCapLog };
    }
    // p = log2(rate_CFL1 / rate_CFL05) per stage
    for (const k of ["b1", "b2"]) for (const w of ["contactWindow", "first1s", "onsetWindow", "growthWindow"]) {
      const a = out.cases[k + "_cfl1"][w], b = out.cases[k + "_cfl05"][w];
      if (!a || !b) continue;
      const p = {};
      for (const st of Object.keys(a.rows)) { const ra = a.rows[st].A_rate, rb = b.rows[st].A_rate; p[st] = ra && rb && ra / rb > 0 ? Math.log2(ra / rb) : null; }
      out.cases[k + "_cfl1"][w].p_A = p;
    }
    const OUT = path.join(ROOT, "validation");
    fs.writeFileSync(path.join(OUT, `item1_stage_budget${TAG}.json`), JSON.stringify(out, null, 1));
    const lines = ["case,window,stage,A_perSubstep,A_rate,B_perSubstep,B_rate,p_A"];
    for (const [cn, c] of Object.entries(out.cases)) for (const w of ["contactWindow", "first1s", "onsetWindow", "growthWindow"]) { const W = c[w]; if (!W) continue; for (const [st, v] of Object.entries(W.rows)) lines.push([cn, w, JSON.stringify(st), v.A_perSubstep, v.A_rate, v.B_perSubstep, v.B_rate, W.p_A ? W.p_A[st] : ""].join(",")); lines.push([cn, w, "\"MEASURED dE\"", W.dE_measured / W.substeps, W.dE_measured / W.duration, "", "", ""].join(",")); }
    fs.writeFileSync(path.join(OUT, `item1_stage_budget${TAG}.csv`), lines.join("\n") + "\n");
    // per-substep series for plotting/forensics (compact)
    const serLines = ["case,t0,dt,dE,minDist_h,vmax," + STAGES.map((s) => JSON.stringify("A:" + s)).join(",")];
    for (const r of R) for (const q of r.series) { const agg1 = {}; for (const [k, v] of Object.entries(q.rows)) { const st = stageOf(k); agg1[st] = (agg1[st] || 0) + v[0]; } const fx = (x, d) => (x == null || !Number.isFinite(x) ? String(x) : x.toExponential(d)); serLines.push([r.case, q.t0.toFixed(6), fx(q.dt, 4), fx(q.dE, 6), q.minDist_h == null ? "" : q.minDist_h.toFixed(3), q.vmax == null ? "" : q.vmax.toFixed(3), ...STAGES.map((s) => fx(agg1[s] || 0, 6))].join(",")); }
    fs.writeFileSync(path.join(OUT, `item1_stage_series${TAG}.csv`), serLines.join("\n") + "\n");
    for (const [cn, c] of Object.entries(out.cases)) {
      console.log(`\n== ${cn}  nan:${c.nanGuard}  tContact:${c.tContact}  caps:${c.capHits}  vmax:${c.maxParticleSpeed}  pcg1s:${JSON.stringify(c.pcgFirst1s)}`);
      for (const w of ["contactWindow", "first1s", "onsetWindow", "growthWindow"]) { const W = c[w]; if (!W) continue; console.log(`  [${w}] n=${W.substeps} dE=${W.dE_measured} closureA=${W.closureA} closureB=${W.closureB}`);
        for (const [st, v] of Object.entries(W.rows)) if (v.A_total || v.B_total) console.log(`    ${st.padEnd(70)} A/sub ${v.A_perSubstep.toExponential(3).padStart(11)}  A rate ${v.A_rate.toExponential(3).padStart(11)}  B rate ${v.B_rate.toExponential(3).padStart(11)}${W.p_A && W.p_A[st] != null ? "  p=" + W.p_A[st].toFixed(2) : ""}`); }
    }
  }).catch((e) => { console.error(e); process.exit(1); });
}
