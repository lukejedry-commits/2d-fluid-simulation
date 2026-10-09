// Item 1, Phase D: regression of the fixed build (plain sim.js, no patches, speed clamp ON).
// Usage: node test/item1_regression.js            -> validation/item1_regression.json
//   1. test/validate.js outputs with suffix _item1 (run `node test/validate.js --suffix _item1` first)
//      byte-compared with the citation dump (validation/*_citation.csv, summary_citation.json; the
//      summary is compared key by key, wall-clock "cost" fields reported separately).
//   2. Run 1 (displaced depth) via test/obstacle_gates.js child cases r1_none / r1_box:
//      delta d = d_box - d_none must lie in [0.0381, 0.0849] (frozen band).
//   3. b2 (2-cell barrier) and b1 (1-cell) on the plain sim.js: dam break, damping 0, CFL 1, 30 s:
//      leaks, min SDF, E(30) < E0, substep-cap hits, max speed.
const fs = require("fs"), path = require("path"), { spawn } = require("child_process");
const ROOT = path.join(__dirname, ".."), V = path.join(ROOT, "validation");
const NX = 128, NY = 72, H = 3 / 128, G = 9.81, FLOOR = (NY + 1) * H, T_END = 30;
const CASE = (() => { const i = process.argv.indexOf("--case"); return i >= 0 ? process.argv[i + 1] : null; })();

function barrierCase(kind) {
  const { FlipSolver: S } = require(path.join(ROOT, "sim.js"));
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  const [a, b] = kind === "b1" ? [64, 65] : [64, 66];
  S.addObstacle({ type: "rect", cx: 0.5 * (a + b) * H, cy: 0.5 * (NY + 2) * H, w: (b - a) * H, h: (NY + 2) * H, angle: 0 });
  S.reset("dam");
  const fy = (S.ny - 1) * H, xR = b * H, n0 = S.count;
  const E = () => { let e = 0; for (let i = 0; i < S.count; i++) e += 0.5 * (S.vx[i] ** 2 + S.vy[i] ** 2) + G * (fy - S.py[i]); return e / S.count; };
  const E0 = E(); let maxE = E0, maxEt = 0, maxDown = 0, minSdf = Infinity, nan = 0, vmax = 0;
  const ever = new Uint8Array(n0);
  S.setSubstepHook(() => {
    let down = 0;
    for (let i = 0; i < S.count; i++) {
      if (S.px[i] > xR) { down++; ever[i] = 1; }
      const sd = S.obstacleSdf(S.px[i], S.py[i]); if (sd < minSdf) minSdf = sd;
    }
    if (down > maxDown) maxDown = down;
  });
  while (S.simTime < T_END - 1e-9) {
    S.step();
    const e = E(); if (e > maxE) { maxE = e; maxEt = S.simTime; }
    const st = S.stats(); nan = Math.max(nan, st.nan); vmax = Math.max(vmax, st.vmax);
  }
  let everDown = 0; for (let i = 0; i < n0; i++) everDown += ever[i];
  return { case: kind, build: "plain sim.js (speed clamp ON)", damping: 0, cfl: S.params.cfl, simTime: S.simTime, E0, E30: E(), maxE, maxEt, nan,
    capHits: S.capLog.length, pcgCapHits: S.pcgCapLog.length, vmaxFrames: vmax, leaks: { maxDownstream: maxDown, everDownstream: everDown }, minSdfOverR: minSdf / S.radius };
}

if (CASE) { process.stdout.write(JSON.stringify(barrierCase(CASE))); process.exit(0); }

function child(args) {
  return new Promise((res, rej) => {
    const p = spawn(process.execPath, args, { cwd: ROOT, stdio: ["ignore", "pipe", "pipe"] });
    let out = "", err = ""; p.stdout.on("data", (d) => (out += d)); p.stderr.on("data", (d) => (err += d));
    p.on("close", (c) => (c === 0 ? res(JSON.parse(out)) : rej(new Error(args.join(" ") + " exit " + c + "\n" + err))));
  });
}

(async () => {
  const out = { generated: new Date().toString() };
  // 1. validate.js vs citation
  const files = ["comparison_mm", "div_h1", "div_h2", "energy_h1", "energy_h2", "front_h1", "front_h2", "volume_h1", "volume_h2"];
  const cmp = {};
  for (const f of files) {
    const A = path.join(V, f + "_citation.csv"), B = path.join(V, f + "_item1.csv");
    cmp[f + ".csv"] = fs.existsSync(B) ? (fs.readFileSync(A).equals(fs.readFileSync(B)) ? "identical" : "DIFFERENT") : "missing (run validate.js --suffix _item1)";
  }
  const sA = JSON.parse(fs.readFileSync(path.join(V, "summary_citation.json"))), sBp = path.join(V, "summary_item1.json");
  const diffs = [], cost = [];
  if (fs.existsSync(sBp)) {
    const sB = JSON.parse(fs.readFileSync(sBp));
    const walk = (x, y, p) => {
      if (x && y && typeof x === "object" && typeof y === "object") { for (const k of new Set([...Object.keys(x), ...Object.keys(y)])) walk(x[k], y[k], p + "." + k); }
      else if (JSON.stringify(x) !== JSON.stringify(y)) (/\.cost\./.test(p) ? cost : diffs).push({ key: p, citation: x, item1: y });
    };
    walk(sA, sB, "");
  }
  out.validate = { files: cmp, summaryDiffsExcludingCost: diffs, summaryCostDiffs: cost.length + " wall-clock fields differ (timing only)",
    identical: Object.values(cmp).every((v) => v === "identical") && diffs.length === 0 && fs.existsSync(sBp) };
  // 2 + 3 in parallel child processes
  const og = path.join(ROOT, "test", "obstacle_gates.js"), me = __filename;
  const [r1n, r1b, b2, b1] = await Promise.all([child([og, "--case", "r1_none"]), child([og, "--case", "r1_box"]), child([me, "--case", "b2"]), child([me, "--case", "b1"])]);
  const dd = r1b.d_last5s - r1n.d_last5s;
  out.run1 = { d_none: r1n.d_last5s, d_box: r1b.d_last5s, delta_d: dd, band: [0.0381, 0.0849], pass: dd >= 0.0381 && dd <= 0.0849, box_minSdfOverR: r1b.minSdfOverR, box_insideMax: r1b.insideMax, nan: Math.max(r1n.nanMax, r1b.nanMax) };
  out.b2 = Object.assign(b2, { pass: b2.leaks.everDownstream === 0 && b2.minSdfOverR >= 0 && b2.E30 < b2.E0 && b2.nan === 0 });
  out.b1 = Object.assign(b1, { pass: b1.leaks.everDownstream === 0 && b1.minSdfOverR >= 0 && b1.E30 < b1.E0 && b1.nan === 0 });
  fs.writeFileSync(path.join(V, "item1_regression.json"), JSON.stringify(out, null, 2));
  console.log(JSON.stringify({ validateIdentical: out.validate.identical, files: cmp, run1: out.run1, b2: out.b2, b1: out.b1 }, null, 1));
})().catch((e) => { console.error(e); process.exit(1); });
