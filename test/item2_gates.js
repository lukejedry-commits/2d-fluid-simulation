// Item 2 (push-out PE lift) - Run 2 frozen gates: dam break past the floor box / the ramp (one
// obstacle per run), CFL 1, substep cap 8, damping 0, 30 s, reviewed defaults.
// Usage: node test/item2_gates.js [--sim path] [--tag NAME] [--cases a,b,...]
//          parent: runs each case in a fresh child process (children never spawn anything)
//        node test/item2_gates.js --case NAME [--sim path]   (child; JSON on stdout)
// Cases: box, ramp (instrumented, CFL 1); box_plain, ramp_plain (plain sim.js, final-state hash
//        must equal the instrumented run: measurement is read-only); box_cfl05, ramp_cfl05.
// Budget (method A of test/item1_budget.js, read-only G2P probe): every stage books the change of
// the particle energy (particle stages) or of the probe energy (grid stages), so the rows add up
// exactly to the measured particle dE of each substep. Bookkeeping adopted for item 2: P2G and the
// obstacle ghost fill are ONE combined row ("P2G + obstacle ghost fill").
// Gate rows: "push-out PE" + "push-out KE (excl. removal)" (= push-out row) and the removal rows
// "velocity removal (pre-advection strip)" (new code) + "velocity removal (push-out backstop)".
// Gate-4 counters: every velocity change made by limitObstacleApproach() is checked against the
// obstacle normal from outside the function (tangential part of dv, outward normal part).
const fs = require("fs"), os = require("os"), path = require("path"), crypto = require("crypto"), { spawn } = require("child_process");
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf("--" + k); return i >= 0 ? args[i + 1] : d; };
const ROOT = path.join(__dirname, ".."), V = path.join(ROOT, "validation");
const SIM = path.resolve(opt("sim", path.join(ROOT, "sim.js"))), TAG = opt("tag", "");
const NX = 128, NY = 72, H = 3 / 128, G = 9.81, T_END = +opt("tend", "30"), FLOOR = (NY + 1) * H;
const CASE = opt("case", null);
const CASES = { box: ["box", 1, true], ramp: ["ramp", 1, true], box_plain: ["box", 1, false], ramp_plain: ["ramp", 1, false], box_cfl05: ["box", 0.5, true], ramp_cfl05: ["ramp", 0.5, true] };

function load(cfl, instrument) {
  if (cfl === 1 && !instrument) return require(SIM);
  let src = fs.readFileSync(SIM, "utf8");
  const patch = (a, b) => { if (src.split(a).length !== 2) throw new Error("anchor not unique/found: " + a.slice(0, 80)); src = src.replace(a, b); };
  if (cfl !== 1) patch("  const CFL = 1.0;", "  const CFL = " + cfl + ";");
  if (!instrument) return req(src, cfl, "p");
  const hasClip = src.includes("clipSepPush(X[i]");
  const hasLimiter = src.includes("    if (obsActive) limitObstacleApproach(dt);\n");
  patch("  // ---- time stepping ------------------------------------------------------\n", `
  // ===== [item2 harness] stage budget instrumentation (method A) =====
  let IB_SUB = null, IB_E = [0, 0], IB_PLACE = "", IB_SELF = false, IB_CK = null, IB_PE = 0;
  function ibPart() { let k = 0, q = 0; const fy = (ny - 1) * h, g = params.gravity; for (let i = 0; i < count; i++) { k += 0.5 * (pvx[i] * pvx[i] + pvy[i] * pvy[i]); q += g * (fy - py[i]); } return [k / count, q / count]; }
  let IB_TX = null, IB_TY = null;
  function ibProbeC(G, GO, stride, offX, offY, VEL, out) {
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
      out[p] = pic * vel + flip * (VEL[p] + dv);
    }
  }
  function ibProbe() {
    if (!IB_TX || IB_TX.length < count) { IB_TX = new Float64Array(maxN); IB_TY = new Float64Array(maxN); }
    ibProbeC(u, IB_SELF ? null : uOld, nx + 1, 0, 0.5, pvx, IB_TX);
    ibProbeC(v, IB_SELF ? null : vOld, nx, 0.5, 0, pvy, IB_TY);
    let k = 0; for (let p = 0; p < count; p++) k += 0.5 * (IB_TX[p] * IB_TX[p] + IB_TY[p] * IB_TY[p]);
    return k / count;
  }
  function ibBook(label, a) { if (!IB_SUB) return; IB_SUB[label] = (IB_SUB[label] || 0) + a; }
  function IBP(label) { if (!IB_SUB) return; const e = ibPart(), s = e[0] + e[1]; ibBook(label, s - IB_E[0]); IB_E = [s, s]; IB_PE = e[1]; }
  function IBG(label) { if (!IB_SUB) return; const a = ibProbe() + IB_PE; ibBook(label, a - IB_E[0]); IB_E = [a, a]; }
  function IBGP(label) { if (!IB_SUB) return; IBG((IB_PLACE ? IB_PLACE + "." : "") + label); }
  function ibCkSave() { IBGP("other"); IB_CK = { sub: Object.assign({}, IB_SUB), e: IB_E.slice() }; }
  function ibCkRestore() { if (!IB_SUB || !IB_CK) return; IB_SUB = Object.assign({}, IB_CK.sub); IB_E = IB_CK.e.slice(); }
  const IB_TOT = {}; const IB_SER = [];
  const IP = { pushed: 0, depthSum: 0, depthMax: 0, vr: 0, stuck0: 0 };
  const IL = { events: 0, full: 0, partial: 0, tangIn: 0, outIn: 0, tangExt: 0, outExt: 0, changedParticles: 0, maxRelTang: 0 };
  const IW = { writes: 0, dKEt: 0, sig2: 0, sig2exact: 0, sumVt2: 0 }; // gate 4b (stored values, per-particle sums)
  const IC = { pairs: 0, clipped: 0, maxMagDiffRel: 0, maxComOverH: 0, maxGrow: 0, storedPairViol: 0, storedPairMaxRel: 0, a: null, velChanged: 0, sumX: 0, sumY: 0, sig2X: 0, sig2Y: 0 };
  const IBand = { advN: 0, advDepthMax: 0, sepN: 0, sepDepthMax: 0 };
  let IB_S0X = null, IB_S0Y = null, IB_V0X = null, IB_V0Y = null;
  const ibO = { nx: 0, ny: 0 };
  function ibBand() { let n = 0, dm = 0; for (let i = 0; i < count; i++) { const d = obstacleSdfAt(px[i], py[i], pr, ibO); if (d < pr) { n++; if (pr - d > dm) dm = pr - d; } } return [n, dm]; }
  // ===== end =====
  // ---- time stepping ------------------------------------------------------\n`);
  patch("    if (mirror) enforceBoundary();\n    uOld.set(u);\n    vOld.set(v);\n",
    "    IB_SELF = true; IBG(\"P2G\"); IB_PLACE = \"p2g\";\n    if (mirror) enforceBoundary();\n    uOld.set(u);\n    vOld.set(v);\n    IB_SELF = false;\n");
  patch("  function enforceBoundary() {\n    enforceWalls();\n    if (obsActive) enforceObstacles(); // obstacle faces follow the same rule as wall faces\n  }",
    "  function enforceBoundary() {\n    IBGP(\"pending\");\n    enforceWalls(); IBGP(\"enforce.walls\");\n    if (obsActive) enforceObstacles(); // obstacle faces follow the same rule as wall faces\n  }");
  { const a = "    fillGhosts(u, ghostU, uW);\n    fillGhosts(v, ghostV, vW);\n  }";
    patch(a, "    IBGP(\"enforce.obstacleInterface\");\n" + a.replace(/\n  \}$/, "") + "\n    IBGP(\"enforce.obstacleGhostFill\");\n  }"); }
  patch("    uSnap.set(u);\n    vSnap.set(v);\n", "    uSnap.set(u);\n    vSnap.set(v);\n    ibCkSave();\n");
  patch("      u.set(uSnap); v.set(vSnap); enforceBoundary(); solvePressurePCG(dt);", "      u.set(uSnap); v.set(vSnap); ibCkRestore(); enforceBoundary(); solvePressurePCG(dt);");
  patch("  function solvePressurePCG(dt) {\n", "  function solvePressurePCG(dt) {\n    IBGP(\"other\");\n    solvePressurePCG_raw(dt);\n    IBG(\"projection\");\n  }\n  function solvePressurePCG_raw(dt) {\n");
  // push-out backstop: count + depth + normal-velocity removal (read-only)
  patch("      if (d >= pr) continue;\n      let vx = VX[i], vy = VY[i];\n", "      if (d >= pr) continue;\n      IP.pushed++; IP.depthSum += pr - d; if (pr - d > IP.depthMax) IP.depthMax = pr - d;\n      let vx = VX[i], vy = VY[i];\n");
  patch("        if (vn < 0) { vx -= vn * o.nx; vy -= vn * o.ny; }\n", "        if (vn < 0) { const e0 = vx * vx + vy * vy; vx -= vn * o.nx; vy -= vn * o.ny; IP.vr += 0.5 * (vx * vx + vy * vy - e0) / count; }\n");
  if (hasLimiter) // gate-4 check inside the new code (normal o used by the code itself)
    patch("        if (vn < lim) { const dv = lim - vn; vx += dv * o.nx; vy += dv * o.ny; changed = true; }\n",
      "        if (vn < lim) { const v0x = vx, v0y = vy, sp = Math.hypot(vx, vy); const dv = lim - vn; vx += dv * o.nx; vy += dv * o.ny; changed = true;\n" +
      "          IL.events++; if (d <= rc) IL.full++; else IL.partial++; const ddx = vx - v0x, ddy = vy - v0y, dn = ddx * o.nx + ddy * o.ny, tx = ddx - dn * o.nx, ty = ddy - dn * o.ny;\n" +
      "          if (Math.hypot(tx, ty) > 1e-12 * sp) IL.tangIn++; if (Math.abs(Math.max(vx * o.nx + vy * o.ny, 0) - Math.max(vn, 0)) > 1e-12 * sp) IL.outIn++; }\n");
  if (hasClip) { patch('              if (obsActive) { clipSepPush(X[i], Y[i], X[j], Y[j], dx, dy); dx = sepClipOut.x; dy = sepClipOut.y; }\n', '              if (obsActive) { const a0 = [X[i], Y[i], X[j], Y[j]], ox = dx, oy = dy; clipSepPush(X[i], Y[i], X[j], Y[j], dx, dy); dx = sepClipOut.x; dy = sepClipOut.y;\n                IC.pairs++; if (dx !== ox || dy !== oy) IC.clipped++; IC.a = a0;\n                { const mi = Math.hypot(-dx, -dy), mj = Math.hypot(dx, dy), m0 = Math.hypot(ox, oy), md = Math.abs(mi - mj) / (Math.max(mi, mj) || 1), com = Math.hypot(-dx + dx, -dy + dy) / h;\n                  if (md > IC.maxMagDiffRel) IC.maxMagDiffRel = md; if (com > IC.maxComOverH) IC.maxComOverH = com; if (m0 > 0 && mi / m0 > IC.maxGrow) IC.maxGrow = mi / m0; } }\n'); patch('              X[i] -= dx;\n              Y[i] -= dy;\n              X[j] += dx;\n              Y[j] += dy;\n            } else if (hi) {', '              X[i] -= dx;\n              Y[i] -= dy;\n              X[j] += dx;\n              Y[j] += dy;\n              if (obsActive && IC.a) { const a = IC.a, dAx = X[i] - a[0], dAy = Y[i] - a[1], dBx = X[j] - a[2], dBy = Y[j] - a[3];\n                const q = Math.pow(2, -24), tol = q * (Math.abs(a[0]) + Math.abs(a[2])), tolY = q * (Math.abs(a[1]) + Math.abs(a[3])), sx = dAx + dBx, sy = dAy + dBy, ex = Math.abs(sx), ey = Math.abs(sy);\n                if (ex > tol || ey > tolY) IC.storedPairViol++; const rr = Math.max(tol ? ex / tol : 0, tolY ? ey / tolY : 0); if (rr > IC.storedPairMaxRel) IC.storedPairMaxRel = rr;\n                IC.sumX += sx; IC.sumY += sy; IC.sig2X += q * q / 3 * (X[i] * X[i] + X[j] * X[j]); IC.sig2Y += q * q / 3 * (Y[i] * Y[i] + Y[j] * Y[j]); IC.a = null; }\n            } else if (hi) {'); }
  const sub = src.match(/  function substep\(dt\) \{\n[\s\S]*?\n  \}\n/);
  const body = sub[0].split("\n").slice(1, -2).map((l) => l.trim());
  const expected = ["pinHeld();", ...(hasLimiter ? ["if (obsActive) limitObstacleApproach(dt);"] : []), "advectParticles(dt);", "separateParticles(SEPARATION_ITERS);", "collideWalls();", "if (obsActive) pushOutObstacles();",
    "particlesToGrid();", "markCells();", "applyGravity(dt);", "applyViscosity();", "enforceBoundary();", "solvePressure(dt);",
    "extrapolateVelocities();", "enforceBoundary();", "if (substepHook) substepHook(divergenceStats(dt), dt);", "gridToParticles();"];
  if (JSON.stringify(body) !== JSON.stringify(expected)) throw new Error("substep() changed:\n" + body.join("\n"));
  patch(sub[0], `  function substep(dt) {
    for (let i = 0; i < count; i++) if (!(Number.isFinite(px[i]) && Number.isFinite(py[i]) && Number.isFinite(pvx[i]) && Number.isFinite(pvy[i]))) throw new Error("NaN guard: particle " + i + " non-finite at t=" + simTime);
    pinHeld();
    IB_SUB = {}; IB_PLACE = "";
    { const e = ibPart(); IB_E = [e[0] + e[1], e[0] + e[1]]; IB_PE = e[1]; }
    const E0 = IB_E[0]; let nStrip = 0;
    ${hasLimiter ? `if (obsActive) {
      if (!IB_V0X || IB_V0X.length < count) { IB_V0X = new Float64Array(maxN); IB_V0Y = new Float64Array(maxN); }
      IB_V0X.set(pvx.subarray(0, count)); IB_V0Y.set(pvy.subarray(0, count));
      limitObstacleApproach(dt);
      for (let i = 0; i < count; i++) {
        const dx = pvx[i] - IB_V0X[i], dy = pvy[i] - IB_V0Y[i]; if (dx === 0 && dy === 0) continue;
        nStrip++; IL.changedParticles++;
        obstacleSdfAt(px[i], py[i], 1e9, ibO); // independent normal (single-obstacle runs)
        const sp = Math.hypot(IB_V0X[i], IB_V0Y[i]), dn = dx * ibO.nx + dy * ibO.ny, tx = dx - dn * ibO.nx, ty = dy - dn * ibO.ny, rt = Math.hypot(tx, ty) / (sp || 1);
        { const tX = -ibO.ny, tY = ibO.nx, vt0 = IB_V0X[i] * tX + IB_V0Y[i] * tY, vt1 = pvx[i] * tX + pvy[i] * tY, q2 = Math.pow(2, -48);
          IW.writes++; IW.dKEt += 0.5 * (vt1 * vt1 - vt0 * vt0); IW.sumVt2 += vt0 * vt0 * (IB_V0X[i] ** 2 + IB_V0Y[i] ** 2); IW.sig2exact += q2 / 3 * vt0 * vt0 * (tX * tX * IB_V0X[i] ** 2 + tY * tY * IB_V0Y[i] ** 2); }
        if (rt > IL.maxRelTang) IL.maxRelTang = rt; if (rt > 1e-12) IL.tangExt++;
        const o0 = Math.max(IB_V0X[i] * ibO.nx + IB_V0Y[i] * ibO.ny, 0), o1 = Math.max(pvx[i] * ibO.nx + pvy[i] * ibO.ny, 0); if (Math.abs(o1 - o0) > 1e-12 * sp) IL.outExt++;
      }
      IBP("velocity removal (pre-advection strip)");
    }` : ""}
    advectParticles(dt); IBP("gravity+advection");
    if (obsActive) { if (!IB_S0X || IB_S0X.length < count) { IB_S0X = new Float32Array(maxN); IB_S0Y = new Float32Array(maxN); } IB_S0X.set(pvx.subarray(0, count)); IB_S0Y.set(pvy.subarray(0, count)); }
    let bA = [0, 0]; if (obsActive) { bA = ibBand(); IBand.advN += bA[0]; if (bA[1] > IBand.advDepthMax) IBand.advDepthMax = bA[1]; }
    separateParticles(SEPARATION_ITERS); IBP("separation (particles)");
    if (obsActive) for (let i = 0; i < count; i++) if (pvx[i] !== IB_S0X[i] || pvy[i] !== IB_S0Y[i]) IC.velChanged++;
    collideWalls(); IBP("wall clamp (particles)");
    let bS = [0, 0]; if (obsActive) { bS = ibBand(); IBand.sepN += bS[0]; if (bS[1] > IBand.sepDepthMax) IBand.sepDepthMax = bS[1]; }
    let pPE = 0, pKE = 0, pVR = 0; const np0 = IP.pushed;
    if (obsActive) {
      const k0 = ibPart(); IP.vr = 0;
      pushOutObstacles();
      const k1 = ibPart(); pPE = k1[1] - k0[1]; pVR = IP.vr; pKE = k1[0] - k0[0] - pVR;
      ibBook("push-out PE", pPE); ibBook("push-out KE (excl. removal)", pKE); ibBook("velocity removal (push-out backstop)", pVR);
      IB_E = [k1[0] + k1[1], k1[0] + k1[1]]; IB_PE = k1[1];
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
    const e = ibPart(), dE = e[0] + e[1] - E0;
    let sumRows = 0; for (const k in IB_SUB) { IB_TOT[k] = (IB_TOT[k] || 0) + IB_SUB[k]; sumRows += IB_SUB[k]; }
    IB_SER.push([simTime + dt, dt, dE, sumRows, pPE, pKE, pVR, IB_SUB["velocity removal (pre-advection strip)"] || 0, IP.pushed - np0, nStrip, bA[0], bS[0], IB_SUB["separation (particles)"] || 0]);
    IB_SUB = null;
  }
`);
  patch("    setSubstepHook(fn) { substepHook = fn || null; },", "    setSubstepHook(fn) { substepHook = fn || null; },\n    get ibTot() { return IB_TOT; }, get ibSer() { return IB_SER; }, get ibPush() { return IP; }, get ibLim() { return IL; }, get ibW() { return IW; }, get ibClip() { return IC; }, get ibBandStats() { return IBand; }, hasLimiter: " + hasLimiter + ",");
  return req(src, cfl, "i");
}
function req(src, cfl, tag) { const tmp = path.join(os.tmpdir(), `sim_item2_${process.pid}_${cfl}_${tag}.js`); fs.writeFileSync(tmp, src); const m = require(tmp); fs.unlinkSync(tmp); return m; }

function runCase(name) {
  const [kind, cfl, inst] = CASES[name];
  const mod = load(cfl, inst), S = mod.FlipSolver, OS = mod.ObstacleShapes;
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  if (kind === "box") S.addObstacle({ type: "rect", cx: (70 + 21) * H, cy: FLOOR - 4 * H, w: 42 * H, h: 8 * H, angle: 0 }); // cells 70..111, 8 high (Run 1/2 box)
  else S.addObstacle(OS.ramp(112 * H, FLOOR - 8 * H, 70 * H, FLOOR)); // Run 2 ramp, same footprint, rises to the right
  S.reset("dam");
  const n0 = S.count, fy = (S.ny - 1) * H;
  const E = () => { let k = 0, q = 0; for (let i = 0; i < S.count; i++) { k += 0.5 * (S.vx[i] ** 2 + S.vy[i] ** 2); q += G * (fy - S.py[i]); } return (k + q) / S.count; };
  const E0 = E();
  let insideMax = 0, minSdf = Infinity, countChanged = false, subs = 0, err = null;
  S.setSubstepHook(() => { // substep end (positions final; G2P does not move particles)
    subs++; let ins = 0;
    for (let i = 0; i < S.count; i++) { const d = S.obstacleSdf(S.px[i], S.py[i]); if (d < minSdf) minSdf = d; if (d < 0) ins++; }
    if (ins > insideMax) insideMax = ins; if (S.count !== n0) countChanged = true;
  });
  try { while (S.simTime < T_END - 1e-9) S.step(); } catch (x) { err = String(x.message || x); }
  const st = S.stats();
  const hsh = crypto.createHash("sha256"); for (const A of [S.px, S.py, S.vx, S.vy]) hsh.update(Buffer.from(A.buffer, A.byteOffset, S.count * 4));
  const res = { case: name, kind, cfl, sim: SIM, instrumented: inst, nanGuard: err, nan: st.nan, simTime: S.simTime, particles: [n0, S.count], countChanged, obstacleDeleted: st.obstacleDeleted, pushStuck: st.pushStuck,
    insideMaxAtSubstepEnd: insideMax, minSdfOverR: minSdf / S.radius, capHits: S.capLog.length, pcgCapHits: S.pcgCapLog.length, substeps: subs, E0, E30: E(), finalStateHash: hsh.digest("hex").slice(0, 16) };
  if (!inst) return res;
  res.hasLimiter = S.hasLimiter; res.rows = S.ibTot; res.push = S.ibPush; res.limiter = S.ibLim; res.gate4b = S.ibW; res.clip = S.ibClip; res.count = S.count; res.band = S.ibBandStats;
  res.series = S.ibSer;
  return res;
}

if (CASE) { if (!CASES[CASE]) { console.error("unknown case " + CASE); process.exit(2); } const r = runCase(CASE); process.stdout.write(JSON.stringify(r), () => process.exit(0)); }
else {
  const names = (opt("cases", Object.keys(CASES).join(","))).split(",");
  const child = (c) => new Promise((ok, no) => {
    const f = path.join(os.tmpdir(), `item2_${TAG || "cur"}_${c}.json`), fd = fs.openSync(f, "w");
    const p = spawn(process.execPath, ["--max-old-space-size=4096", __filename, "--case", c, "--sim", SIM, "--tend", String(T_END)], { stdio: ["ignore", fd, "inherit"] });
    p.on("close", (code) => { fs.closeSync(fd); code === 0 ? ok(JSON.parse(fs.readFileSync(f, "utf8"))) : no(new Error(c + " exited " + code)); });
  });
  Promise.all(names.map(child)).then((R) => {
    const out = Object.fromEntries(R.map((r) => [r.case, r]));
    fs.writeFileSync(path.join(os.tmpdir(), `item2_${TAG || "cur"}_all.json`), JSON.stringify(out));
    console.error("wrote " + path.join(os.tmpdir(), `item2_${TAG || "cur"}_all.json`));
  }).catch((e) => { console.error(e); process.exit(1); });
}
