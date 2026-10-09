// Diagnostic for Run 3 (report-only): total mechanical energy per unit mass vs time for the
// full-height barriers of test/obstacle_gates.js (and a no-barrier baseline). Fresh Node
// process per case, damping 0, reviewed defaults, 30 s. Plain sim.js (CFL 1), no changes.
// Writes validation/thin_barrier_energy_obst.json and .csv.
const path = require("path"), fs = require("fs"), { spawn } = require("child_process");
const NX = 128, NY = 72, H = 3 / 128, T_END = 30;
const B = { none: null, b05c: [64.25, 64.75], b05e: [64.6, 65.1], b1: [64, 65], b2: [64, 66] };
const k = process.argv[2] === "--case" ? process.argv[3] : null;
if (k) {
  const { FlipSolver: S } = require(path.join(__dirname, "..", "sim.js"));
  S.init(NX, NY, H); S.params.damping = 0; console.warn = () => {};
  const b = B[k];
  if (b) S.addObstacle({ type: "rect", cx: 0.5 * (b[0] + b[1]) * H, cy: 0.5 * (NY + 2) * H, w: (b[1] - b[0]) * H, h: (NY + 2) * H, angle: 0 });
  S.reset("dam");
  const fy = (S.ny - 1) * H;
  const E = () => { let ke = 0, pe = 0, vmax = 0, fast = 0; for (let i = 0; i < S.count; i++) { const v2 = S.vx[i] ** 2 + S.vy[i] ** 2; ke += 0.5 * v2; pe += 9.81 * (fy - S.py[i]); vmax = Math.max(vmax, Math.sqrt(v2)); if (v2 > 25) fast++; } return { t: S.simTime, ke: ke / S.count, e: (ke + pe) / S.count, vmax, fast }; };
  const rows = [E()]; let next = 1;
  while (S.simTime < T_END - 1e-9) { S.step(); if (S.simTime >= next - 1e-9) { rows.push(E()); next++; } }
  process.stdout.write(JSON.stringify({ case: k, barrier_cells: b, capHits: S.capLog.length, nan: S.stats().nan, rows }), () => process.exit(0));
} else {
  const run = (c) => new Promise((res, rej) => { const p = spawn(process.execPath, [__filename, "--case", c], { stdio: ["ignore", "pipe", "inherit"] }); let o = ""; p.stdout.on("data", (d) => (o += d)); p.on("close", (code) => (code ? rej(new Error(c + " " + code)) : res(JSON.parse(o)))); });
  Promise.all(Object.keys(B).map(run)).then((R) => {
    const out = path.join(__dirname, "..", "validation");
    const summary = R.map((r) => ({ case: r.case, barrier_cells: r.barrier_cells, capHits: r.capHits, nan: r.nan, E0: r.rows[0].e, E_max: Math.max(...r.rows.map((q) => q.e)), E30: r.rows[r.rows.length - 1].e, KE30: r.rows[r.rows.length - 1].ke, vmax_max: Math.max(...r.rows.map((q) => q.vmax)), fastOver5mps_max: Math.max(...r.rows.map((q) => q.fast)) }));
    fs.writeFileSync(path.join(out, "thin_barrier_energy_obst.json"), JSON.stringify({ note: "E = KE + g*(height above floor), J/kg, sampled every 1 s sim time; speed clamp is 14 m/s", summary, series: R }, null, 1));
    const lines = ["t," + R.map((r) => `E_${r.case},KE_${r.case}`).join(",")];
    for (let i = 0; i < R[0].rows.length; i++) lines.push(R[0].rows[i].t.toFixed(3) + "," + R.map((r) => (r.rows[i] ? r.rows[i].e.toFixed(5) + "," + r.rows[i].ke.toFixed(5) : ",")).join(","));
    fs.writeFileSync(path.join(out, "thin_barrier_energy_obst.csv"), lines.join("\n") + "\n");
    console.log(JSON.stringify(summary, null, 1));
  });
}
