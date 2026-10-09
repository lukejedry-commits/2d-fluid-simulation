// Validation runs for the FLIP/APIC solver (Martin & Moyce column collapse).
// Usage: node test/validate.js [--solver pcg|gs] [--sep neutral|legacy] [--boundary mirror|legacy|separating] [--suffix _pcg] [--sim path/to/sim.js] [--only h2|h1]
//   writes validation/<file><suffix>.csv and validation/summary<suffix>.json
// Each column runs in a fresh Node process; resetColumn also reseeds. Not referenced by index.html.
const fs = require("fs");
const path = require("path");
const { spawnSync } = require("child_process");
const args = process.argv.slice(2);
const opt = (k, d) => { const i = args.indexOf("--" + k); return i >= 0 ? args[i + 1] : d; };
const SOLVER = opt("solver", "pcg");
const SEP = opt("sep", "neutral");
const BOUNDARY = opt("boundary", "separating");
const SUFFIX = opt("suffix", "_" + SOLVER);
const ONLY = opt("only", "");
const EM = require("./energy_metrics.js");
const SIM = opt("sim", path.join(__dirname, "..", "sim.js"));
const OUT = path.join(__dirname, "..", "validation");
fs.mkdirSync(OUT, { recursive: true });

if (!ONLY) {
  // Parent: spawn one fresh process per column, then merge summaries / comparison CSVs.
  const self = path.resolve(__filename);
  const baseArgs = [];
  for (let i = 0; i < args.length; i++) {
    if (args[i] === "--only") { i++; continue; }
    baseArgs.push(args[i]);
  }
  function child(tag) {
    const r = spawnSync(process.execPath, [self, ...baseArgs, "--only", tag], { encoding: "utf8", maxBuffer: 80e6 });
    if (r.status !== 0) {
      process.stderr.write(r.stderr || "");
      process.stderr.write(r.stdout || "");
      process.stderr.write("validate.js child " + tag + " exited " + r.status + "\n");
      process.exit(r.status || 1);
    }
    return JSON.parse(r.stdout);
  }
  const a = child("h2"), b = child("h1");
  const summary = {
    h2: a.h2, h1: b.h1,
    comparison: { h2: a.comparison.h2, h1: b.comparison.h1 },
    harness: { reseed: true, freshProcessPerCase: true },
  };
  const rows = [];
  for (const tag of ["h2", "h1"]) {
    const pth = path.join(OUT, "comparison_mm_" + tag + SUFFIX + ".csv");
    if (fs.existsSync(pth)) {
      rows.push(...fs.readFileSync(pth, "utf8").trim().split("\n").slice(1));
      fs.unlinkSync(pth);
    }
  }
  fs.writeFileSync(path.join(OUT, "comparison_mm" + SUFFIX + ".csv"),
    "case,clock,shift_T,ref_a_in,T_ref,T_sim_release,Z_ref,Z_sim,ratio,result,note\n" + rows.join("\n") + "\n");
  fs.writeFileSync(path.join(OUT, "summary" + SUFFIX + ".json"), JSON.stringify(summary, null, 2));
  const brief = (c) => Object.fromEntries(["release", "shifted"].map((k) => [k, Object.fromEntries(Object.entries(c[k]).map(([aa, ss]) => [aa, { points: ss.points, pass: ss.pass, minRatio: +ss.minRatio.toFixed(4), minAtT: ss.minAtT, minOffWall: ss.minOffWall }]))]));
  console.log(JSON.stringify({
    h2: summary.h2, h1: summary.h1,
    comparison: {
      h2: { T144: summary.comparison.h2.T144, shift: summary.comparison.h2.shift, ...brief(summary.comparison.h2) },
      h1: { T144: summary.comparison.h1.T144, shift: summary.comparison.h1.shift, ...brief(summary.comparison.h1) },
    },
    harness: summary.harness,
  }, null, 1));
  process.exit(0);
}

const { FlipSolver: S } = require(path.resolve(SIM));

const NX = 128, NY = 72, H_CELL = 3 / 128; // tank 3 m x 1.6875 m
const G = 9.81;
const SIM_SECONDS = 30;
const FRONT_DT = 0.02;

const sorted = (arr) => Float64Array.from(arr).sort();
function percentile(arr, p) {
  const a = sorted(arr);
  if (!a.length) return NaN;
  return a[Math.min(a.length - 1, Math.max(0, Math.ceil(p * a.length) - 1))];
}
const avg = (x) => (x.length ? x.reduce((s, v) => s + v, 0) / x.length : NaN);
function csv(file, header, rows) {
  fs.writeFileSync(path.join(OUT, file.replace(".csv", SUFFIX + ".csv")), header + "\n" + rows.map((r) => r.join(",")).join("\n") + "\n");
}

function run(cellsA, cellsH, tag) {
  S.init(NX, NY, H_CELL);
  S.params.damping = 0; // validation spec
  S.params.pressureSolver = SOLVER;
  S.params.separation = SEP;
  S.params.boundary = BOUNDARY;
  const p = S.params;
  const config = { solver: p.pressureSolver, flipRatio: p.flipRatio, gsIters: p.pressureIters, pcgMaxIters: p.pcgMaxIters, pcgRelTol: S.PCG_REL_TOL, pcgAbsTol: S.PCG_ABS_TOL, viscosity: p.viscosity, gravity: p.gravity, omegaGS: S.OVER_RELAX, cfl: S.CFL, maxSub: S.MAX_SUBSTEPS, damping: p.damping, separation: p.separation, boundary: p.boundary };
  S.resetColumn(cellsA, cellsH);
  S.capLog.length = 0;
  S.pcgCapLog.length = 0;

  const h = S.cellSize;
  const a = cellsA * h;
  const Hc = cellsH * h;
  const floorY = (S.ny - 1) * h;
  const wallX = h;
  const tankW = (S.nx - 2) * h;
  const xFrontMax = tankW - S.radius;
  const Tscale = Math.sqrt(G * Hc) / a;
  const N0count = S.count;

  function frontX() {
    let m = -1;
    for (let i = 0; i < S.count; i++) {
      if (floorY - S.py[i] <= h) {
        const x = S.px[i] - wallX;
        if (x > m) m = x;
      }
    }
    return m;
  }
  const frontRaw = [{ t: 0, x: frontX() }];
  let frontDone = false;
  let frontReachT = null;
  const divRows = [];
  const subIters = [];
  let sub = 0;
  let hookMs = 0;
  let maxCflSub = 0, maxCflFullSub = 0;
  // separation-pass centre-of-mass shift per substep (y down; lift = -dy)
  const sepC = { maxAbsRawDy: 0, maxAbsRawDx: 0, maxAbsEffDy: 0, maxAbsEffDx: 0, maxEffLift: 0, maxEffLiftT: null, sumEffLift: 0, sumEffLiftPreImpact: 0 };
  const wallA = { maxPasses: 0, sumPasses: 0, n: 0, capHits: 0, capNegHits: 0, worstNegP: 0, worstUInto: 0, compFailSteps: 0, totalCompFail: 0 };
  S.setSubstepHook(function (d, dt) {
    const th = performance.now();
    sub++;
    subIters.push(d.pcgIters);
    if (d.cfl > maxCflSub) maxCflSub = d.cfl;
    if (d.cflFull > maxCflFullSub) maxCflFullSub = d.cflFull;
    sepC.maxAbsRawDy = Math.max(sepC.maxAbsRawDy, Math.abs(d.sepRawDy));
    sepC.maxAbsRawDx = Math.max(sepC.maxAbsRawDx, Math.abs(d.sepRawDx));
    sepC.maxAbsEffDy = Math.max(sepC.maxAbsEffDy, Math.abs(d.sepEffDy));
    sepC.maxAbsEffDx = Math.max(sepC.maxAbsEffDx, Math.abs(d.sepEffDx));
    if (-d.sepEffDy > sepC.maxEffLift) { sepC.maxEffLift = -d.sepEffDy; sepC.maxEffLiftT = d.t; }
    sepC.sumEffLift += -d.sepEffDy;
    if (!frontDone) sepC.sumEffLiftPreImpact += -d.sepEffDy;
    wallA.n++; wallA.sumPasses += d.sepPasses; if (d.sepPasses > wallA.maxPasses) wallA.maxPasses = d.sepPasses;
    if (d.sepCapHit) wallA.capHits++;
    if (d.sepCapHit && d.sepNegPLeft) wallA.capNegHits++;
    if (d.sepWorstNegP < wallA.worstNegP) wallA.worstNegP = d.sepWorstNegP;
    if (d.sepWorstUInto > wallA.worstUInto) wallA.worstUInto = d.sepWorstUInto;
    if (d.sepCompFail) { wallA.compFailSteps++; wallA.totalCompFail += d.sepCompFail; }
    divRows.push([sub, d.t.toFixed(6), dt.toExponential(4), d.fluidCells, d.maxDiv.toExponential(4), d.meanDiv.toExponential(4), d.maxRes.toExponential(4), d.meanRes.toExponential(4), d.pcgIters, d.pcgRelRes.toExponential(3), d.pcgCapHit ? 1 : 0, d.cfl.toFixed(4), d.sepRawDy.toExponential(3), d.sepEffDy.toExponential(3)]);
    if (!frontDone) {
      const x = frontX();
      frontRaw.push({ t: d.t, x });
      if (x >= xFrontMax - 1e-6) { frontDone = true; frontReachT = d.t; }
    }
    hookMs += performance.now() - th;
  });

  const energyRows = [], volumeRows = [];
  const tArr = [], Etot = [], KEa = [], PEa = [];
  const N0 = cellsA * cellsH;
  let maxRelN = 0, maxRelNt = 0, N10 = null, N30 = null;
  const rhoSamples = [];
  let rhoMin = Infinity, rhoMinT = 0, rhoMax = -Infinity, rhoMaxT = 0, rho30 = null, rhoOut = 0;
  // [2b] Validation grades density on the legacy sampler; also track renorm (physics).
  let rhoMinLeg = Infinity, rhoMinLegT = 0, rhoMaxLeg = -Infinity, rhoMaxLegT = 0, rho30Leg = null, rhoOutLeg = 0;
  let rhoMinRen = Infinity, rhoMaxRen = -Infinity, rho30Ren = null;
  let nanMax = 0, escMax = 0;
  const depthSamples = [];
  let depth30 = null, meanHeight30 = null;
  const countEvents = [];
  let lastCount = S.count;
  let frameIdx = 0;
  const frameMs = [];

  // particle pairs closer than 0.8 r and 1.2 r (floor-row overpacking metric)
  function closePairs() {
    const n = S.count, r = S.radius, cs = 1.6 * r, gx = Math.ceil((S.nx * h) / cs) + 2, map = new Map();
    for (let i = 0; i < n; i++) { const k = Math.floor(S.px[i] / cs) + gx * Math.floor(S.py[i] / cs); let L = map.get(k); if (!L) map.set(k, (L = [])); L.push(i); }
    let p08 = 0, p12 = 0, p12floor = 0;
    for (let i = 0; i < n; i++) {
      const cx = Math.floor(S.px[i] / cs), cy = Math.floor(S.py[i] / cs);
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
        const L = map.get(cx + dx + gx * (cy + dy)); if (!L) continue;
        for (const j of L) { if (j <= i) continue; const d = Math.hypot(S.px[j] - S.px[i], S.py[j] - S.py[i]);
          if (d < 0.8 * r) p08++; if (d < 1.2 * r) { p12++; if (floorY - S.py[i] < h && floorY - S.py[j] < h) p12floor++; } }
      }
    }
    let nb = 0; for (let i = 0; i < n; i++) if (floorY - S.py[i] < h) nb++;
    return { lt0p8r: p08, lt1p2r: p12, lt1p2rBothInFloorRow: p12floor, floorRowParticlesPerCellOver4: nb / (4 * (S.nx - 2)) };
  }
  const pairs = {};
  function columnDepth() {
    // free-surface height of the floor-connected water per column (FLUID cells walked
    // up from the floor; highest particle centre in the top cell + 0.25 h, i.e. half
    // the 0.5 h lattice spacing one particle represents), averaged over all columns
    const nx = S.nx, ny = S.ny, T = S.cellType;
    const topY = new Float64Array(nx).fill(Infinity);
    const topRow = new Int32Array(nx).fill(-1);
    for (let i = 1; i < nx - 1; i++) {
      let j = ny - 2;
      if (T[i + j * nx] !== 1) continue;
      while (j > 1 && T[i + (j - 1) * nx] === 1) j--;
      topRow[i] = j;
    }
    for (let k = 0; k < S.count; k++) {
      const i = Math.floor(S.px[k] / h);
      const j = Math.floor(S.py[k] / h);
      if (topRow[i] === j && S.py[k] < topY[i]) topY[i] = S.py[k];
    }
    let sum = 0;
    for (let i = 1; i < nx - 1; i++) sum += topY[i] === Infinity ? 0 : floorY - topY[i] + 0.25 * h;
    return sum / (nx - 2);
  }

  {
    let pe0 = 0;
    for (let i = 0; i < S.count; i++) pe0 += G * (floorY - S.py[i]);
    pe0 /= S.count;
    tArr.push(0); KEa.push(0); PEa.push(pe0); Etot.push(pe0);
    energyRows.push(["0.00000", 0, "0.000000", pe0.toFixed(6), pe0.toFixed(6), S.count, 0, "0.0000"]);
  }
  while (S.simTime < SIM_SECONDS - 1e-9) {
    const h0 = hookMs;
    const t0 = performance.now();
    S.step();
    frameMs.push(performance.now() - t0 - (hookMs - h0)); // solver cost without diagnostics
    frameIdx++;
    const st = S.stats();
    const t = S.simTime;
    let ke = 0, pe = 0;
    for (let i = 0; i < S.count; i++) {
      ke += 0.5 * (S.vx[i] * S.vx[i] + S.vy[i] * S.vy[i]);
      pe += G * (floorY - S.py[i]);
    }
    ke /= S.count;
    pe /= S.count;
    tArr.push(t); KEa.push(ke); PEa.push(pe); Etot.push(ke + pe);
    energyRows.push([t.toFixed(5), frameIdx, ke.toFixed(6), pe.toFixed(6), (ke + pe).toFixed(6), S.count, st.substeps, st.umax.toFixed(4)]);
    const N = st.fluidCells;
    const rel = (N - N0) / N0;
    if (Math.abs(rel) > maxRelN) { maxRelN = Math.abs(rel); maxRelNt = t; }
    if (N10 === null && t >= 10 - 1e-9) N10 = N;
    if (t >= SIM_SECONDS - 1e-9) N30 = N;
    // Grade on legacy sampler (2b); physics may use renorm.
    const rhoLeg = (st.interiorDensityLegacy != null ? st.interiorDensityLegacy : st.interiorDensity) / st.restDensity;
    const rhoRen = st.interiorDensity / st.restDensity;
    const rho = rhoLeg; // Validation grades on legacy
    rhoSamples.push(rho);
    if (rho < rhoMin) { rhoMin = rho; rhoMinT = t; }
    if (rho > rhoMax) { rhoMax = rho; rhoMaxT = t; }
    if (rho < 0.98 || rho > 1.02) rhoOut++;
    rho30 = rho;
    if (rhoLeg < rhoMinLeg) { rhoMinLeg = rhoLeg; rhoMinLegT = t; }
    if (rhoLeg > rhoMaxLeg) { rhoMaxLeg = rhoLeg; rhoMaxLegT = t; }
    if (rhoLeg < 0.98 || rhoLeg > 1.02) rhoOutLeg++;
    rho30Leg = rhoLeg;
    if (rhoRen < rhoMinRen) rhoMinRen = rhoRen;
    if (rhoRen > rhoMaxRen) rhoMaxRen = rhoRen;
    rho30Ren = rhoRen;
    nanMax = Math.max(nanMax, st.nan);
    escMax = Math.max(escMax, st.outside);
    if (S.count !== lastCount) countEvents.push({ t, from: lastCount, to: S.count });
    lastCount = S.count;
    volumeRows.push([t.toFixed(5), frameIdx, N, rel.toFixed(5), rho.toFixed(5), S.count, st.nan, st.outside]);
    if (t >= SIM_SECONDS - 5 - 1e-9) depthSamples.push(columnDepth());
    for (const tp of [1.2, 5, 30]) if (Math.abs(t - tp) < 0.5 / 60 && !pairs["t" + tp]) pairs["t" + tp] = closePairs();
    if (t >= SIM_SECONDS - 1e-9) {
      depth30 = columnDepth();
      let sh = 0;
      for (let i = 0; i < S.count; i++) sh += floorY - S.py[i];
      meanHeight30 = sh / S.count;
    }
  }
  S.setSubstepHook(null);

  // front on the 0.02 s grid (linear interpolation between substep ends, <= 1 cell apart)
  const front = [];
  const tEnd = frontReachT !== null ? frontReachT : frontRaw[frontRaw.length - 1].t;
  let k = 0;
  for (let n = 0; ; n++) {
    const t = n * FRONT_DT;
    if (t > tEnd + 1e-12) break;
    while (k < frontRaw.length - 2 && frontRaw[k + 1].t < t) k++;
    const A = frontRaw[k], B = frontRaw[k + 1] || A;
    const w = B.t > A.t ? Math.min(1, Math.max(0, (t - A.t) / (B.t - A.t))) : 0;
    front.push({ t, x: A.x + w * (B.x - A.x) });
  }
  if (frontReachT !== null && front[front.length - 1].t < frontReachT - 1e-9) front.push({ t: frontReachT, x: frontRaw[frontRaw.length - 1].x, arrival: true });
  const frontRows = [];
  let maxdZdT = 0, maxdZdTat = null;
  const dzOver = [];
  for (let i = 0; i < front.length; i++) {
    const f = front[i];
    f.Z = f.x / a; f.T = f.t * Tscale;
    let d = "";
    if (i > 0) {
      const v = (f.Z - front[i - 1].Z) / (f.T - front[i - 1].T);
      d = v.toFixed(4);
      if (v > 2) dzOver.push({ T: +f.T.toFixed(3), dZdT: +v.toFixed(4) });
      if (v > maxdZdT) { maxdZdT = v; maxdZdTat = f.T; }
    }
    frontRows.push([f.t.toFixed(5), f.T.toFixed(5), f.x.toFixed(5), f.Z.toFixed(5), d, f.arrival ? "arrival" : ""]);
  }
  // T_1.44: first sim time with Z >= 1.44 (linear interpolation between samples)
  let T144 = null;
  for (let i = 1; i < front.length; i++) {
    if (front[i].Z >= 1.44) {
      const A = front[i - 1], B = front[i];
      T144 = A.Z >= 1.44 ? A.T : A.T + ((1.44 - A.Z) * (B.T - A.T)) / (B.Z - A.Z);
      break;
    }
  }
  csv(`front_${tag}.csv`, "t_s,T,x_front_m,Z,dZdT,note", frontRows);
  csv(`energy_${tag}.csv`, "t_s,frame,KE_J_per_kg,PE_J_per_kg,total_J_per_kg,particles,substeps,u_max_m_per_s", energyRows);
  csv(`volume_${tag}.csv`, "t_s,frame,fluid_cells,rel_change_vs_N0,interior_density_over_rest,particles,nan,escaped", volumeRows);
  csv(`div_${tag}.csv`, "substep,t_s,dt_s,fluid_cells,raw_max_absdiv_dt,raw_mean_absdiv_dt,solver_max_residual_dt,solver_mean_residual_dt,pcg_iters,pcg_rel_residual,pcg_cap_hit,cfl,sep_com_dy_raw_m,sep_com_dy_after_wall_clamp_m", divRows);

  // energy: early peak + sustained-rise metric (0.5 s centred moving average vs running minimum, t >= 0.5 s)
  let pk = 0;
  for (let i = 0; i < Etot.length && tArr[i] <= 1; i++) if (Etot[i] > Etot[pk]) pk = i;
  const win = 30;
  const sm = Etot.map((_, i) => { let s = 0, n = 0; for (let j = Math.max(0, i - win / 2); j <= Math.min(Etot.length - 1, i + win / 2); j++) { s += Etot[j]; n++; } return s / n; });
  let runMin = Infinity, maxRise = 0, maxRiseT = 0;
  for (let i = 0; i < sm.length; i++) {
    if (tArr[i] < 0.5) continue;
    if (sm[i] < runMin) runMin = sm[i];
    if (sm[i] - runMin > maxRise) { maxRise = sm[i] - runMin; maxRiseT = tArr[i]; }
  }
  // unsmoothed rises above the running minimum (frame-level, 60 Hz), in three windows
  const riseIn = (t0, t1) => EM.riseIn(tArr, Etot, t0, t1);
  const reach = frontReachT !== null ? frontReachT : 1;
  const riseUnsmoothed = { preImpact_0p2_to_wall: riseIn(0.2, reach), impactPlus0p5_to_5s: riseIn(reach + 0.5, 5), after5s: riseIn(5, 1e9) };
  // individual post-impact spikes, wall impact .. 5 s (zig-zag, 0.02 % E0 reversal, >= 0.05 % E0)
  const postImpactSpikes = EM.spikes(tArr, Etot, reach, 5);
  // every total-energy upswing > 0.02 % E0 between 8 and 30 s (0.01 % E0 reversal)
  const lateUpswings = EM.spikes(tArr, Etot, 8, 30, 0.01, 0.02);
  const col = (i) => divRows.map((r) => +r[i]);
  const rawMax = col(4), rawMean = col(5), resMax = col(6), resMean = col(7);
  const at = (tq) => { let i = tArr.findIndex((t) => t >= tq - 1e-9); if (i < 0) i = tArr.length - 1; return { t: tArr[i], KE: KEa[i], PE: PEa[i], E: Etot[i] }; };
  const fs30 = sorted(frameMs.slice(60));

  return {
    tag, cellsA, cellsH, a, H: Hc, Tscale, h, config,
    particles: N0count, particlesEnd: S.count, countEvents, nanRecoveries: S.nanRecoveries, nanMax, escMax,
    front: { reachT: frontReachT, reachTstar: frontReachT !== null ? frontReachT * Tscale : null, samples: front.length, maxdZdT, maxdZdTat, over2: dzOver, T144, series: front.map((f) => ({ t: f.t, T: f.T, Z: f.Z })) },
    energy: { E0: Etot[0], earlyPeak: { t: tArr[pk], E: Etot[pk], risePct: (100 * (Etot[pk] - Etot[0])) / Etot[0], backBelowE0At: tArr.find((t, i) => t > tArr[pk] && Etot[i] <= Etot[0]) }, maxRiseAfter05_smoothed: maxRise, maxRiseAfter05T: maxRiseT, maxRiseAfter05Pct: (100 * maxRise) / Etot[0], smoothedNote: "0.5 s centred moving average vs running min from t=0.5 s; the window straddles the wall impact, so it hides the pre-impact rise. Use riseUnsmoothed.", riseUnsmoothed, postImpactSpikes, lateUpswings_8_30s: lateUpswings, at1: at(1), at5: at(5), at10: at(10), at20: at(20), at30: at(30) },
    volume: { N0, maxRelN, maxRelNt, N10, N30, drift10_30: (N30 - N10) / N0, absDrift10_30: Math.abs(N30 - N10) / N0, rhoMin, rhoMinT, rhoMax, rhoMaxT, rho30, rhoSamples: rhoSamples.length, rhoOutOfBand: rhoOut, densitySamplers: { gradeOn: "legacy", legacy: { rhoMin: rhoMinLeg, rhoMinT: rhoMinLegT, rhoMax: rhoMaxLeg, rhoMaxT: rhoMaxLegT, rho30: rho30Leg, rhoOutOfBand: rhoOutLeg }, renorm: { rhoMin: rhoMinRen, rhoMax: rhoMaxRen, rho30: rho30Ren } } },
    depth: { last5Mean: avg(depthSamples), last5Min: Math.min(...depthSamples), last5Max: Math.max(...depthSamples), at30: depth30, fromParticleVolume: (S.count * h * h * 0.25) / tankW, target: (a * Hc) / tankW },
    height: { meanParticleHeight30: meanHeight30, halfLast5Depth: avg(depthSamples) / 2, halfDepth30: depth30 / 2, E30: at(30).E, gTimesMeanHeight: G * meanHeight30 },
    div: { substeps: divRows.length, rawMax: Math.max(...rawMax), rawP95OfMax: percentile(rawMax, 0.95), rawMean: avg(rawMean), solverMax: Math.max(...resMax), solverP95OfMax: percentile(resMax, 0.95), solverMean: avg(resMean) },
    pcg: { mean: avg(subIters), p95: percentile(subIters, 0.95), max: Math.max(...subIters), capHits: S.pcgCapLog.length, capLog: S.pcgCapLog.slice(0, 50) },
    cost: { meanMsPerFrame: avg(frameMs.slice(60)), p95MsPerFrame: percentile(frameMs.slice(60), 0.95), maxMsPerFrame: fs30[fs30.length - 1], note: "excludes first 60 frames (JIT warm-up) and diagnostics-hook time" },
    cfl: { maxAdvective: maxCflSub, maxWithGravityTerm: maxCflFullSub },
    closePairs: pairs,
    wallActiveSet: {
      ...wallA,
      meanPasses: wallA.n ? wallA.sumPasses / wallA.n : 0,
      p95Passes: (() => { const a = Float64Array.from(S.sepPassLog).sort(); return a.length ? a[Math.min(a.length - 1, Math.ceil(0.95 * a.length) - 1)] : 0; })(),
      pass: BOUNDARY !== "separating" ? null : (wallA.capNegHits === 0 && wallA.worstUInto <= 1e-6 && wallA.totalCompFail === 0 && wallA.worstNegP >= -1e-12),
      note: "separating-wall active set (2g only). Contact p = cell p when constrained, 0 when released. Null pass when boundary != separating.",
    },
    separationCom: { ...sepC, g_x_sumEffLift_J_per_kg: G * sepC.sumEffLift, g_x_sumEffLiftPreImpact_J_per_kg: G * sepC.sumEffLiftPreImpact, note: "per-substep centre-of-mass shift of the free particles caused by the push-apart pass (m, y down). raw = right after the pass; eff = wall clamp applied to the positions before and after the pass. lift = -dy." },
    cap: { hits: S.capLog.length, log: S.capLog.slice(0, 200) },
    frames: frameIdx, simTime: S.simTime,
  };
}

const tag = ONLY;
const spec = tag === "h2" ? [24, 48, "h2", 2, 1.19] : [24, 24, "h1", 1, 0.8];
const runResult = run(spec[0], spec[1], spec[2]);

const ref = fs.readFileSync(path.join(OUT, "mm_reference.csv"), "utf8").trim().split("\n").slice(1).map((l) => {
  const c = l.split(",");
  return { n2: +c[1], a_in: c[2], Z: +c[3], T: +c[4] };
});
function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + ((ys[i] - ys[i - 1]) * (x - xs[i - 1])) / (xs[i] - xs[i - 1]);
  return null;
}
function compare(run, n2, Tref) {
  const shift = Tref - run.front.T144;
  const sT = run.front.series.map((f) => f.T), sZ = run.front.series.map((f) => f.Z);
  const series = {};
  for (const p of ref.filter((q) => q.n2 === n2)) (series[p.a_in] = series[p.a_in] || []).push(p);
  const out = { T144: run.front.T144, Tref, shift, release: {}, shifted: {} };
  for (const [clock, off] of [["release", 0], ["shifted", shift]]) {
    for (const [aIn, pts] of Object.entries(series)) {
      pts.sort((x, y) => x.T - y.T);
      const rows = [];
      for (const p of pts) {
        if (p.T < 0.5 || p.Z > 5.33) continue;
        const tSim = p.T - off;
        let zs = interp(sT, sZ, tSim);
        let note = "";
        if (zs === null) { zs = sZ[sZ.length - 1]; note = "sim front already at far wall"; }
        rows.push({ T: p.T, Tsim: tSim, Zref: p.Z, Zsim: zs, ratio: zs / p.Z, pass: zs >= 0.9 * p.Z, note });
      }
      const min = rows.reduce((m, r) => (r.ratio < m.ratio ? r : m), rows[0]);
      const offWall = rows.filter((r) => !r.note);
      const minOff = offWall.length ? offWall.reduce((m, r) => (r.ratio < m.ratio ? r : m), offWall[0]) : null;
      out[clock][aIn] = { points: rows.length, pass: rows.every((r) => r.pass), minRatio: min.ratio, minAtT: min.T, minAtZref: min.Zref, minOffWall: minOff && { ratio: minOff.ratio, T: minOff.T, Zref: minOff.Zref }, rows };
    }
  }
  return out;
}
const cmpOne = compare(runResult, spec[3], spec[4]);
const comparisonRows = [];
for (const clock of ["release", "shifted"])
  for (const [aIn, ss] of Object.entries(cmpOne[clock]))
    for (const r of ss.rows) comparisonRows.push([tag, clock, clock === "shifted" ? cmpOne.shift.toFixed(4) : 0, aIn, r.T, r.Tsim.toFixed(4), r.Zref, r.Zsim.toFixed(4), r.ratio.toFixed(4), r.pass ? "pass" : "FAIL", r.note]);
// per-case comparison file; parent merges into comparison_mm<suffix>.csv
fs.writeFileSync(path.join(OUT, `comparison_mm_${tag}${SUFFIX}.csv`), "case,clock,shift_T,ref_a_in,T_ref,T_sim_release,Z_ref,Z_sim,ratio,result,note\n" + comparisonRows.map((r) => r.join(",")).join("\n") + "\n");
const strip = (r) => ({ ...r, front: { ...r.front, series: undefined }, cap: { hits: r.cap.hits }, pcg: { ...r.pcg, capLog: undefined } });
const payload = { comparison: {} };
payload[tag] = strip(runResult);
payload.comparison[tag] = cmpOne;
console.log(JSON.stringify(payload));
