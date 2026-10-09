// Fills the E(t) series missing from the obstacle_gates outputs (needed for plots only).
// Same setups as test/obstacle_gates.js (Run 2 cases at CFL 1 / 0.5, and the 1-cell barrier at
// CFL 0.5) with the plain sim.js (CFL patched in an OS-temp copy, nothing else). Each case is a
// fresh Node process; E = KE + g*(height above floor) per unit mass, sampled every frame.
// The final-state hash / cap-hit count is compared with the original gate run to prove the
// re-run reproduces it. Writes validation/energy_series_obst.json.
const fs = require("fs"), os = require("os"), path = require("path"), crypto = require("crypto"), { spawn } = require("child_process");
const ROOT = path.join(__dirname, ".."), SIM = path.join(ROOT, "sim.js");
const NX = 128, NY = 72, H = 3 / 128, T_END = 30, G = 9.81, FLOOR = (NY + 1) * H;
function load(cfl) {
  if (cfl === 1) return require(SIM);
  let s = fs.readFileSync(SIM, "utf8");
  if (s.split("  const CFL = 1.0;").length !== 2) throw new Error("CFL anchor");
  s = s.replace("  const CFL = 1.0;", "  const CFL = " + cfl + ";");
  const tmp = path.join(os.tmpdir(), `sim_obst_series_${process.pid}.js`); fs.writeFileSync(tmp, s);
  const m = require(tmp); fs.unlinkSync(tmp); return m;
}
const CASES = {
  r2_none: [1, "none"], r2_box: [1, "box"], r2_ramp: [1, "ramp"],
  r2_none_cfl05: [0.5, "none"], r2_box_cfl05: [0.5, "box"], r2_ramp_cfl05: [0.5, "ramp"],
  r3_b1_cfl05: [0.5, "b1"],
};
const name = process.argv[2] === "--case" ? process.argv[3] : null;
if (name) {
  const [cfl, kind] = CASES[name], mod = load(cfl), S = mod.FlipSolver, O = mod.ObstacleShapes;
  console.warn = () => {};
  S.init(NX, NY, H); S.params.damping = 0;
  if (kind === "box") S.addObstacle({ type: "rect", cx: 91 * H, cy: FLOOR - 4 * H, w: 42 * H, h: 8 * H, angle: 0 });
  if (kind === "ramp") S.addObstacle(O.ramp(112 * H, FLOOR - 8 * H, 70 * H, FLOOR));
  if (kind === "b1") S.addObstacle({ type: "rect", cx: 64.5 * H, cy: 0.5 * (NY + 2) * H, w: H, h: (NY + 2) * H, angle: 0 });
  S.reset("dam");
  const fy = (S.ny - 1) * S.cellSize;
  const E = () => { let k = 0, q = 0, vm = 0; for (let i = 0; i < S.count; i++) { const v2 = S.vx[i] ** 2 + S.vy[i] ** 2; k += 0.5 * v2; q += G * (fy - S.py[i]); if (v2 > vm) vm = v2; } return [+S.simTime.toFixed(6), (k + q) / S.count, k / S.count, Math.sqrt(vm), S.capLog.length]; };
  const rows = [E()];
  while (S.simTime < T_END - 1e-9) { S.step(); rows.push(E()); }
  const hsh = crypto.createHash("sha256");
  for (const A of [S.px, S.py, S.vx, S.vy]) hsh.update(Buffer.from(A.buffer, A.byteOffset, S.count * 4));
  process.stdout.write(JSON.stringify({ case: name, cfl, kind, finalStateHash: hsh.digest("hex").slice(0, 16), capHits: S.capLog.length, columns: ["t", "E", "KE", "vmax", "capHitsCum"], rows }), () => process.exit(0));
} else {
  const run = (c) => new Promise((res, rej) => { const p = spawn(process.execPath, [__filename, "--case", c], { stdio: ["ignore", "pipe", "inherit"] }); let o = ""; p.stdout.on("data", (d) => (o += d)); p.on("close", (code) => (code ? rej(new Error(c)) : res(JSON.parse(o)))); });
  Promise.all(Object.keys(CASES).map(run)).then((R) => {
    const cache = path.join(os.tmpdir(), "obst_gates_cache"), out = {};
    for (const r of R) {
      let ref = null; try { ref = JSON.parse(fs.readFileSync(path.join(cache, r.case + ".json"), "utf8")); } catch (e) {}
      r.reference = ref ? { finalStateHash: ref.finalStateHash || null, capHits: ref.capHits, E30: ref.E30 ?? null } : null;
      r.reproduces = ref ? (ref.finalStateHash ? ref.finalStateHash === r.finalStateHash : ref.capHits === r.capHits) : null;
      out[r.case] = r;
      console.log(r.case, "hash", r.finalStateHash, "ref", r.reference && r.reference.finalStateHash, "caps", r.capHits, "ref", r.reference && r.reference.capHits, "E30", r.rows[r.rows.length - 1][1].toFixed(6), "refE30", r.reference && r.reference.E30, "reproduces", r.reproduces);
    }
    fs.writeFileSync(path.join(ROOT, "validation", "energy_series_obst.json"), JSON.stringify(out));
  });
}
