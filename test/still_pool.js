// Hard assert for step 2f: flat pool at rest for 1 s.
// Pass: max |u| <= 1e-4 m/s AND |ΔE| <= 1e-6 J/kg. Fail → stop, do not dam-break.
const fs = require("fs");
const path = require("path");
const { FlipSolver: S } = require(path.join(__dirname, "..", "sim.js"));

function run(mass) {
  S.init(128, 72, 3 / 128);
  S.params.damping = 0;
  S.params.pressureSolver = "pcg";
  S.params.separation = "neutral";
  S.params.boundary = "separating";
  S.params.pressureMass = mass;
  // Lattice flat pool (full width, 12 cells deep) — no spawn jitter.
  S.resetColumn(128, 12);

  const h = S.cellSize;
  const fy = (S.ny - 1) * h;
  const g = S.params.gravity;
  const energy = () => {
    let k = 0, q = 0;
    for (let i = 0; i < S.count; i++) {
      k += 0.5 * (S.vx[i] * S.vx[i] + S.vy[i] * S.vy[i]);
      q += g * (fy - S.py[i]);
    }
    return (k + q) / S.count;
  };
  const maxSpeed = () => {
    let m = 0;
    for (let i = 0; i < S.count; i++) m = Math.max(m, Math.hypot(S.vx[i], S.vy[i]));
    return m;
  };

  const E0 = energy();
  let maxU = 0;
  while (S.simTime < 1 - 1e-12) {
    S.step();
    maxU = Math.max(maxU, maxSpeed());
  }
  const E1 = energy();
  const dE = E1 - E0;
  const passU = maxU <= 1e-4;
  const passE = Math.abs(dE) <= 1e-6;
  return {
    pressureMass: mass,
    particles: S.count,
    t: S.simTime,
    E0, E1, dE, maxAbsU: maxU,
    limits: { maxAbsU: 1e-4, absDE: 1e-6 },
    passU, passE, pass: passU && passE,
    clampAccum: {
      hits: S.clampAccum.hits,
      faces: S.clampAccum.faces,
      fraction: S.clampAccum.faces ? S.clampAccum.hits / S.clampAccum.faces : 0,
    },
  };
}

const weighted = run("weighted");
const equal = run("equal");
const out = {
  test: "still-pool",
  note: "Lattice pool resetColumn(128,12), damping 0, separating+neutral+PCG, 1 s.",
  weighted,
  equal,
  pass: weighted.pass,
  stopDamBreak: !weighted.pass,
};
const outPath = path.join(__dirname, "..", "validation", "still_pool_2f.json");
fs.mkdirSync(path.dirname(outPath), { recursive: true });
fs.writeFileSync(outPath, JSON.stringify(out, null, 2));
console.log(JSON.stringify(out, null, 2));
process.exit(weighted.pass ? 0 : 1);
