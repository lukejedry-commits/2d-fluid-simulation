/* 2D SPH water — vanilla JS. Air is empty space; we only integrate water. */

const canvas = document.getElementById("view");
const ctx = canvas.getContext("2d", { alpha: false });

const hudN = document.getElementById("hud-n");
const hudFps = document.getElementById("hud-fps");
const btnPause = document.getElementById("btn-pause");

// Beat: 100 px = 1 m so gravity is 9.81 * 100 px/s², never 9.81 px/s².
const PIXELS_PER_METER = 100;
// Larger discs + spacing → visible water-column depth when stacked.
const DRAW_RADIUS = 6.5; // px, equal for every particle
const SPACING = 13; // dam-break packing, px (~2*R)
const H_PX = 2.2 * SPACING; // smoothing length
const H_M = H_PX / PIXELS_PER_METER;
const H_PX2 = H_PX * H_PX;
const H_M2 = H_M * H_M;
// Mass is fixed for the default pack spacing at DEFAULT_REST_DENSITY.
// Do NOT scale mass with the rest-density dial — that cancels the dial's effect
// (rho and restDensity move together). Keep mass at 420 packing so startup is
// calm at the default, and sliding rest density actually changes pressure.
const DEFAULT_REST_DENSITY = 420;
const SPACING_M2 = (SPACING / PIXELS_PER_METER) * (SPACING / PIXELS_PER_METER);
const particleMass = DEFAULT_REST_DENSITY * SPACING_M2;

const DEFAULT_N = 700;
const MAX_N = 2500;
const FIXED_DT = 1 / 120;
const MAX_SUBSTEPS = 4;
const MAX_NB = 56;
const MAX_SPEED = 2200; // px/s clamp so a hitch cannot explode the state
const RESTITUTION = 0.28;
const WALL_FRICTION = 0.35;

// Soft particle–particle separation (multi-iter Jacobi/GS).
const COLLIDE_DIST = 2.1 * DRAW_RADIUS;
const COLLIDE_DIST2 = COLLIDE_DIST * COLLIDE_DIST;
const COLLIDE_PUSH = 0.95; // near-full separation over iterations
const COLLIDE_DAMP = 0.14;
const COLLIDE_ITERS_WCSPH = 4; // legacy soft-collide count
const COLLIDE_ITERS_PBF = 2; // light overlap cleanup under constraint solvers
function softCollideIters() {
  return SOLVER === "dfsph" || SOLVER === "pbf" ? COLLIDE_ITERS_PBF : COLLIDE_ITERS_WCSPH;
}
const FLOOR_SPRING_K = 1800; // ay spring (px/s² per px penetration)

// Held water-drop: hex disk with axial radius 3 → 37 particles.
const DROP_HEX_R = 3;
const MAX_THROW = 1600; // px/s clamp on release velocity
// Expected neighbor count in bulk hex pack within h (for free-surface cohesion).
const BULK_NB = 18;
// XSPH velocity blending (Schechter/Bridson): coherent wave crests, less particle noise.
const XSPH_EPS = 0.1;

// ---------------------------------------------------------------------------
// Solver switch — flip locally without a rewrite. Default ships as PBF.
//   "pbf"   = Position Based Fluids density constraint (Macklin/Müller)
//   "dfsph" = lite 2D DFSPH (constant-density + divergence-free loops)
// ---------------------------------------------------------------------------
const SOLVER = "pbf"; // "pbf" | "dfsph"

// PBF: iterative λ density constraint (replaces hard WCSPH pressure push).
const PBF_SOLVER_ITERS = 4; // 3–5 typical
const PBF_LAMBDA_EPS = 1.0e-6; // relaxation in λ denominator
const PBF_TENSILE_K = 0.1; // artificial pressure / tensile instability
const PBF_TENSILE_N = 4;
const PBF_TENSILE_DQ = 0.3; // Δq as fraction of h (W(dq)/W(0) reference)

// Lite DFSPH: simplified constant-density + divergence-free pressure solves.
const DFSPH_DENSITY_ITERS = 3;
const DFSPH_DIVERGENCE_ITERS = 2;
const DFSPH_EPS = 1.0e-6;

// 2D Müller-style kernels (distances in meters). Recomputed from H_M.
const POLY6_COEFF = 4 / (Math.PI * Math.pow(H_M, 8));
const POLY6_SELF = POLY6_COEFF * Math.pow(H_M2, 3);
const SPIKY_DW = -30 / (Math.PI * Math.pow(H_M, 5)); // dW_spiky / dr
const VISC_LAP = 40 / (Math.PI * Math.pow(H_M, 5));

// Typed SoA — not thousands of particle objects.
const x = new Float32Array(MAX_N);
const y = new Float32Array(MAX_N);
const vx = new Float32Array(MAX_N);
const vy = new Float32Array(MAX_N);
const ax = new Float32Array(MAX_N);
const ay = new Float32Array(MAX_N);
const density = new Float32Array(MAX_N);
const pressure = new Float32Array(MAX_N);
const nbIndex = new Uint32Array(MAX_N * MAX_NB);
const nbCount = new Uint16Array(MAX_N);
const nbDist = new Float32Array(MAX_N * MAX_NB);
// Accumulated Jacobi corrections for collide / PBF
const dxCorr = new Float32Array(MAX_N);
const dyCorr = new Float32Array(MAX_N);
// Held-drop flags + offsets from grab point (mouse).
const held = new Uint8Array(MAX_N);
const holdOx = new Float32Array(MAX_N);
const holdOy = new Float32Array(MAX_N);
// Solver scratch: PBF λ / DFSPH κ, DFSPH α, positions at step start for v rebuild.
const lambdaArr = new Float32Array(MAX_N);
const alphaArr = new Float32Array(MAX_N);
const x0 = new Float32Array(MAX_N);
const y0 = new Float32Array(MAX_N);

let count = 0;
let width = 1;
let height = 1;
let dpr = 1;
let paused = false;
let gravity = 9.81;
let restDensity = DEFAULT_REST_DENSITY;
let stiffness = 6000;
let viscosity = 0.03;
// Surface tension / cohesion (m/s² scale before × PIXELS_PER_METER).
// Depth-aware: strong on free surface / drops, ~0 in bulk. Default 2.5 beads drops without pinning the lake.
let surfaceTension = 2.5;
let resetTarget = DEFAULT_N;

let grid = [];
let gridCols = 0;
let gridRows = 0;

let mouseDown = false;
let mx = 0;
let my = 0;
let lastMx = 0;
let lastMy = 0;
let lastMouseStamp = 0;
let mvx = 0;
let mvy = 0;
let holding = false;

let accTime = 0;
let lastStamp = 0;
let fpsEma = 60;
let fpsFrames = 0;
let fpsStamp = 0;

function setupCanvas() {
  dpr = Math.max(1, window.devicePixelRatio || 1);
  width = window.innerWidth;
  height = window.innerHeight;
  canvas.width = Math.floor(width * dpr);
  canvas.height = Math.floor(height * dpr);
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function spawnOne(px, py, pvx, pvy) {
  if (count >= MAX_N) return false;
  const i = count++;
  x[i] = px;
  y[i] = py;
  vx[i] = pvx;
  vy[i] = pvy;
  density[i] = restDensity;
  pressure[i] = 0;
  held[i] = 0;
  return true;
}

function spawnParticles(n, mode) {
  // Dam-break: hexagonal lattice — taller column (more freeboard) for wave response.
  clearHold();
  count = 0;
  const want = Math.max(1, Math.min(MAX_N, n | 0));
  const rad = DRAW_RADIUS;
  const left = rad + 12;
  const floor = height - rad - 4;
  const rowPitch = SPACING * 0.8660254; // √3/2 for hex rows
  // Prefer vertical stack over width: ~0.58 canvas for dam, shallow pool mode unchanged.
  const targetH = height * (mode === "pool" ? 0.22 : 0.58);
  const maxRows = Math.max(4, Math.floor(targetH / rowPitch));

  let row = 0;
  let col = 0;
  while (count < want) {
    const ox = (row & 1) * (SPACING * 0.5); // hex offset odd rows
    const px = left + col * SPACING + ox;
    const py = floor - row * rowPitch;
    if (py < rad + 4) {
      // ran out of vertical room — widen horizontally
      row = 0;
      col++;
      continue;
    }
    spawnOne(px, py, 0, 0);
    row++;
    if (row >= maxRows) {
      row = 0;
      col++;
    }
  }
}

function clearHold() {
  if (holding) {
    for (let i = 0; i < count; i++) held[i] = 0;
  }
  holding = false;
}

/** Spawn a compact hex-disk drop centered at (cx,cy) and mark it held. */
function spawnHeldDrop(cx, cy) {
  if (holding) return;
  const rowPitch = SPACING * 0.8660254;
  let any = false;
  // Axial hex coords: |q|+|r|+|s| <= 2*R with s=-q-r → disk of radius DROP_HEX_R.
  for (let q = -DROP_HEX_R; q <= DROP_HEX_R; q++) {
    for (let r = -DROP_HEX_R; r <= DROP_HEX_R; r++) {
      const s = -q - r;
      if (Math.abs(s) > DROP_HEX_R) continue;
      if (count >= MAX_N) break;
      const ox = SPACING * (q + r * 0.5);
      const oy = rowPitch * r;
      if (!spawnOne(cx + ox, cy + oy, 0, 0)) break;
      const i = count - 1;
      held[i] = 1;
      holdOx[i] = ox;
      holdOy[i] = oy;
      any = true;
    }
    if (count >= MAX_N) break;
  }
  holding = any;
}

function releaseHeldDrop() {
  if (!holding) {
    mouseDown = false;
    return;
  }
  // Clamp throw speed so a flick cannot explode the sim.
  let tvx = mvx;
  let tvy = mvy;
  const sp2 = tvx * tvx + tvy * tvy;
  if (sp2 > MAX_THROW * MAX_THROW) {
    const s = MAX_THROW / Math.sqrt(sp2);
    tvx *= s;
    tvy *= s;
  }
  for (let i = 0; i < count; i++) {
    if (!held[i]) continue;
    vx[i] = tvx;
    vy[i] = tvy;
    held[i] = 0;
  }
  holding = false;
  mouseDown = false;
}

/** Pin held particles to mouse + grab offsets; give them mouse velocity. */
function pinHeldParticles() {
  if (!holding) return;
  for (let i = 0; i < count; i++) {
    if (!held[i]) continue;
    x[i] = mx + holdOx[i];
    y[i] = my + holdOy[i];
    vx[i] = mvx;
    vy[i] = mvy;
    ax[i] = 0;
    ay[i] = 0;
  }
}

function buildHash() {
  // Beat: uniform grid / spatial hash. Cell size = h.
  const cols = Math.max(1, Math.ceil(width / H_PX));
  const rows = Math.max(1, Math.ceil(height / H_PX));
  const nCells = cols * rows;
  if (cols !== gridCols || rows !== gridRows || grid.length !== nCells) {
    gridCols = cols;
    gridRows = rows;
    grid = new Array(nCells);
    for (let c = 0; c < nCells; c++) grid[c] = [];
  } else {
    for (let c = 0; c < nCells; c++) grid[c].length = 0;
  }
  for (let i = 0; i < count; i++) {
    let cx = (x[i] / H_PX) | 0;
    let cy = (y[i] / H_PX) | 0;
    if (cx < 0) cx = 0;
    else if (cx >= cols) cx = cols - 1;
    if (cy < 0) cy = 0;
    else if (cy >= rows) cy = rows - 1;
    grid[cy * cols + cx].push(i);
  }
}

function findNeighbors() {
  const cols = gridCols;
  const rows = gridRows;
  for (let i = 0; i < count; i++) {
    let nnb = 0;
    const off = i * MAX_NB;
    let cxi = (x[i] / H_PX) | 0;
    let cyi = (y[i] / H_PX) | 0;
    if (cxi < 0) cxi = 0;
    else if (cxi >= cols) cxi = cols - 1;
    if (cyi < 0) cyi = 0;
    else if (cyi >= rows) cyi = rows - 1;
    const x0 = cxi > 0 ? cxi - 1 : 0;
    const x1 = cxi + 1 < cols ? cxi + 1 : cols - 1;
    const y0 = cyi > 0 ? cyi - 1 : 0;
    const y1 = cyi + 1 < rows ? cyi + 1 : rows - 1;
    const xi = x[i];
    const yi = y[i];
    for (let cy = y0; cy <= y1; cy++) {
      const rowBase = cy * cols;
      for (let cx = x0; cx <= x1; cx++) {
        const cell = grid[rowBase + cx];
        for (let n = 0, len = cell.length; n < len; n++) {
          const j = cell[n];
          if (j === i) continue;
          const dx = xi - x[j];
          const dy = yi - y[j];
          const r2 = dx * dx + dy * dy;
          if (r2 >= H_PX2 || nnb >= MAX_NB) continue;
          nbIndex[off + nnb] = j;
          nbDist[off + nnb] = Math.sqrt(r2) / PIXELS_PER_METER;
          nnb++;
        }
      }
    }
    nbCount[i] = nnb;
  }
}

function computeDensityPressure() {
  // Density always; WCSPH gas-EOS pressure only if somehow needed.
  // Under PBF / DFSPH the constraint solvers own incompressibility — keep p=0
  // so applyForces does not double-push against λ / κ corrections.
  const useEos = false; // pressure EOS path retired; solvers handle density
  for (let i = 0; i < count; i++) {
    let rho = particleMass * POLY6_SELF;
    const off = i * MAX_NB;
    const nnb = nbCount[i];
    for (let k = 0; k < nnb; k++) {
      const r = nbDist[off + k];
      const d = H_M2 - r * r;
      rho += particleMass * POLY6_COEFF * d * d * d;
    }
    density[i] = rho;
    if (useEos) {
      const p = stiffness * (rho / restDensity - 1);
      pressure[i] = p > 0 ? p : 0;
    } else {
      pressure[i] = 0;
    }
  }
}

function applyForces() {
  const gPx = gravity * PIXELS_PER_METER;
  const mu = viscosity;
  const sigma = surfaceTension;
  const floorY = height - DRAW_RADIUS;
  const kFloor = FLOOR_SPRING_K;
  for (let i = 0; i < count; i++) {
    // Held particles are pinned — skip gravity/pressure/viscosity/cohesion.
    if (held[i]) {
      ax[i] = 0;
      ay[i] = 0;
      continue;
    }
    let amx = 0;
    let amy = 0;
    const rhoi = density[i] > 1e-4 ? density[i] : 1e-4;
    const pi = pressure[i];
    const invRhoi2 = 1 / (rhoi * rhoi);
    const off = i * MAX_NB;
    const nnb = nbCount[i];
    const vxi = vx[i];
    const vyi = vy[i];
    const xi = x[i];
    const yi = y[i];
    // Depth-aware cohesion: nnb low (free surface / drop) → strong; nnb~BULK_NB+ → ~0.
    const surfaceW = sigma > 0 ? Math.max(0, 1 - nnb / BULK_NB) : 0;
    for (let k = 0; k < nnb; k++) {
      const j = nbIndex[off + k];
      const r = nbDist[off + k];
      if (r < 1e-8) continue;
      const dxm = (xi - x[j]) / PIXELS_PER_METER;
      const dym = (yi - y[j]) / PIXELS_PER_METER;
      const invr = 1 / r;
      const rx = dxm * invr;
      const ry = dym * invr;
      const rhoj = density[j] > 1e-4 ? density[j] : 1e-4;
      const q = H_M - r;
      const dW = SPIKY_DW * q * q;
      // Pressure EOS term: inert under PBF/DFSPH (pressure[] = 0). Kept for clarity.
      const pTerm = particleMass * (pi * invRhoi2 + pressure[j] / (rhoj * rhoj));
      amx += -pTerm * dW * rx;
      amy += -pTerm * dW * ry;
      const lap = VISC_LAP * q;
      const viscScale = (mu * particleMass) / rhoj * lap;
      amx += (viscScale * (vx[j] - vxi)) / PIXELS_PER_METER;
      amy += (viscScale * (vy[j] - vyi)) / PIXELS_PER_METER;
      // Pairwise cohesion within h; gated by surfaceW so bulk pool is not pinned.
      if (surfaceW > 0) {
        const C = q / H_M; // = 1 - r/H_M
        const coh = sigma * C * surfaceW;
        amx -= coh * rx;
        amy -= coh * ry;
      }
    }
    ax[i] = amx * PIXELS_PER_METER;
    ay[i] = amy * PIXELS_PER_METER + gPx;
    // Soft floor spring in forces (not just a velocity kick after integrate).
    if (yi > floorY) {
      ay[i] -= kFloor * (yi - floorY);
    }
  }
}

/** XSPH: v_i += eps Σ (m_j/ρ_j)(v_j−v_i) W_ij — smooths velocity for coherent crests. */
function applyXsph() {
  const eps = XSPH_EPS;
  if (eps <= 0) return;
  // Scratch into dxCorr/dyCorr (collide/PBF overwrite later).
  for (let i = 0; i < count; i++) {
    if (held[i]) {
      dxCorr[i] = 0;
      dyCorr[i] = 0;
      continue;
    }
    let dvx = 0;
    let dvy = 0;
    const off = i * MAX_NB;
    const nnb = nbCount[i];
    const vxi = vx[i];
    const vyi = vy[i];
    for (let k = 0; k < nnb; k++) {
      const j = nbIndex[off + k];
      const r = nbDist[off + k];
      const d = H_M2 - r * r;
      if (d <= 0) continue;
      const W = POLY6_COEFF * d * d * d;
      const rhoj = density[j] > 1e-4 ? density[j] : 1e-4;
      const s = (particleMass / rhoj) * W;
      dvx += s * (vx[j] - vxi);
      dvy += s * (vy[j] - vyi);
    }
    dxCorr[i] = eps * dvx;
    dyCorr[i] = eps * dvy;
  }
  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    vx[i] += dxCorr[i];
    vy[i] += dyCorr[i];
  }
}

function integrate(dt) {
  // Symplectic Euler: v from forces, XSPH blend, then x.
  for (let i = 0; i < count; i++) {
    if (held[i]) continue; // pinned separately
    let vxi = vx[i] + ax[i] * dt;
    let vyi = vy[i] + ay[i] * dt;
    const speed2 = vxi * vxi + vyi * vyi;
    if (speed2 > MAX_SPEED * MAX_SPEED) {
      const s = MAX_SPEED / Math.sqrt(speed2);
      vxi *= s;
      vyi *= s;
    }
    vx[i] = vxi;
    vy[i] = vyi;
  }
  applyXsph();
  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
  }
}

function collideParticlesOnce() {
  // Pairwise soft separation. CRITICAL: never skip r2≈0 — those must push apart.
  // Scale personal space with rest density so the dial can loosen/tighten packing
  // (mass is fixed; without this, soft collide floors spacing near SPACING).
  const cols = gridCols;
  const rows = gridRows;
  const rhoScale = Math.sqrt(DEFAULT_REST_DENSITY / Math.max(restDensity, 1));
  const minD = COLLIDE_DIST * Math.min(1.85, Math.max(0.55, rhoScale));
  const minD2 = minD * minD;
  const push = COLLIDE_PUSH;
  const damp = COLLIDE_DAMP;
  const halfMin = minD * 0.5;

  for (let i = 0; i < count; i++) {
    dxCorr[i] = 0;
    dyCorr[i] = 0;
  }

  for (let i = 0; i < count; i++) {
    let cxi = (x[i] / H_PX) | 0;
    let cyi = (y[i] / H_PX) | 0;
    if (cxi < 0) cxi = 0;
    else if (cxi >= cols) cxi = cols - 1;
    if (cyi < 0) cyi = 0;
    else if (cyi >= rows) cyi = rows - 1;
    const x0 = cxi > 0 ? cxi - 1 : 0;
    const x1 = cxi + 1 < cols ? cxi + 1 : cols - 1;
    const y0 = cyi > 0 ? cyi - 1 : 0;
    const y1 = cyi + 1 < rows ? cyi + 1 : rows - 1;
    for (let cy = y0; cy <= y1; cy++) {
      const rowBase = cy * cols;
      for (let cx = x0; cx <= x1; cx++) {
        const cell = grid[rowBase + cx];
        for (let n = 0, len = cell.length; n < len; n++) {
          const j = cell[n];
          if (j <= i) continue;
          let dx = x[i] - x[j];
          let dy = y[i] - y[j];
          let r2 = dx * dx + dy * dy;
          if (r2 >= minD2) continue;

          let nx, ny, dist, overlap;
          if (r2 < 1e-8) {
            // Nearly coincident: fallback / random normal, still separate by minD/2 each.
            const ang = Math.random() * Math.PI * 2;
            nx = Math.cos(ang);
            ny = Math.sin(ang);
            dist = 0;
            overlap = minD;
            // Direct half-minD push (full separation target).
            const corr0 = halfMin * push;
            // Don't displace held particles — free partner takes the full correction.
            if (held[i] && !held[j]) {
              dxCorr[j] -= nx * corr0 * 2;
              dyCorr[j] -= ny * corr0 * 2;
            } else if (held[j] && !held[i]) {
              dxCorr[i] += nx * corr0 * 2;
              dyCorr[i] += ny * corr0 * 2;
            } else if (!held[i] && !held[j]) {
              dxCorr[i] += nx * corr0;
              dyCorr[i] += ny * corr0;
              dxCorr[j] -= nx * corr0;
              dyCorr[j] -= ny * corr0;
            }
          } else {
            dist = Math.sqrt(r2);
            nx = dx / dist;
            ny = dy / dist;
            overlap = minD - dist;
            const corr = overlap * 0.5 * push;
            if (held[i] && !held[j]) {
              dxCorr[j] -= nx * corr * 2;
              dyCorr[j] -= ny * corr * 2;
            } else if (held[j] && !held[i]) {
              dxCorr[i] += nx * corr * 2;
              dyCorr[i] += ny * corr * 2;
            } else if (!held[i] && !held[j]) {
              dxCorr[i] += nx * corr;
              dyCorr[i] += ny * corr;
              dxCorr[j] -= nx * corr;
              dyCorr[j] -= ny * corr;
            }
          }

          if (!held[i] && !held[j]) {
            const rvn = (vx[i] - vx[j]) * nx + (vy[i] - vy[j]) * ny;
            if (rvn < 0) {
              const impulse = rvn * damp * 0.5;
              vx[i] -= impulse * nx;
              vy[i] -= impulse * ny;
              vx[j] += impulse * nx;
              vy[j] += impulse * ny;
            }
          }
        }
      }
    }
  }

  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    x[i] += dxCorr[i];
    y[i] += dyCorr[i];
  }
}


function poly6W(r) {
  const d = H_M2 - r * r;
  if (d <= 0) return 0;
  return POLY6_COEFF * d * d * d;
}

/** Spiky gradient magnitude scalar: ∇W = (SPIKY_DW * q²) * r_hat, q=h-r. */
function spikyGradMag(r) {
  if (r < 1e-8 || r >= H_M) return 0;
  const q = H_M - r;
  return SPIKY_DW * q * q; // negative
}

/**
 * Position Based Fluids (Macklin/Müller): iterative λ density constraints.
 * Positions are corrected so ρ ≈ ρ₀; velocities rebuilt from Δx/dt.
 */
function solvePBF(dt) {
  const invDt = dt > 1e-8 ? 1 / dt : 0;
  const rho0 = restDensity;
  const mass = particleMass;
  const wQ = poly6W(PBF_TENSILE_DQ * H_M);
  const w0 = POLY6_SELF > 1e-12 ? POLY6_SELF : 1;
  const invW0 = 1 / w0;

  // Snapshot for velocity rebuild
  for (let i = 0; i < count; i++) {
    x0[i] = x[i];
    y0[i] = y[i];
  }

  for (let iter = 0; iter < PBF_SOLVER_ITERS; iter++) {
    buildHash();
    findNeighbors();

    // Density + λ
    for (let i = 0; i < count; i++) {
      if (held[i]) {
        lambdaArr[i] = 0;
        density[i] = rho0;
        continue;
      }
      let rho = mass * POLY6_SELF;
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      const xi = x[i];
      const yi = y[i];
      for (let k = 0; k < nnb; k++) {
        const r = nbDist[off + k];
        rho += mass * poly6W(r);
      }
      density[i] = rho;
      const C = rho / rho0 - 1;
      if (C < 0) {
        // Allow mild under-density at free surface (no tensile pull from λ alone)
        lambdaArr[i] = 0;
        continue;
      }
      // Σ |∇C|² ≈ Σ |(m/ρ₀) ∇W|² + self term
      let sumGradSq = 0;
      let gradSelfX = 0;
      let gradSelfY = 0;
      for (let k = 0; k < nnb; k++) {
        const j = nbIndex[off + k];
        const r = nbDist[off + k];
        if (r < 1e-8) continue;
        const dxm = (xi - x[j]) / PIXELS_PER_METER;
        const dym = (yi - y[j]) / PIXELS_PER_METER;
        const invr = 1 / r;
        const gMag = spikyGradMag(r); // scalar along r_hat from i→? (xi-xj)
        const gx = gMag * dxm * invr;
        const gy = gMag * dym * invr;
        const scale = mass / rho0;
        const gxj = scale * gx;
        const gyj = scale * gy;
        sumGradSq += gxj * gxj + gyj * gyj;
        gradSelfX += gxj;
        gradSelfY += gyj;
      }
      sumGradSq += gradSelfX * gradSelfX + gradSelfY * gradSelfY;
      lambdaArr[i] = -C / (sumGradSq + PBF_LAMBDA_EPS);
    }

    // Position corrections Δp
    for (let i = 0; i < count; i++) {
      dxCorr[i] = 0;
      dyCorr[i] = 0;
    }
    for (let i = 0; i < count; i++) {
      if (held[i]) continue;
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      const xi = x[i];
      const yi = y[i];
      const li = lambdaArr[i];
      let dpx = 0;
      let dpy = 0;
      for (let k = 0; k < nnb; k++) {
        const j = nbIndex[off + k];
        const r = nbDist[off + k];
        if (r < 1e-8) continue;
        const dxm = (xi - x[j]) / PIXELS_PER_METER;
        const dym = (yi - y[j]) / PIXELS_PER_METER;
        const invr = 1 / r;
        const gMag = spikyGradMag(r);
        const gx = gMag * dxm * invr;
        const gy = gMag * dym * invr;
        // Tensile instability correction (artificial pressure)
        let scorr = 0;
        if (PBF_TENSILE_K > 0 && wQ > 1e-12) {
          const ratio = poly6W(r) * invW0;
          const ratioQ = wQ * invW0;
          if (ratioQ > 1e-12) {
            const s = ratio / ratioQ;
            scorr = -PBF_TENSILE_K * Math.pow(s, PBF_TENSILE_N);
          }
        }
        const coeff = (li + lambdaArr[j] + scorr) * (mass / rho0);
        dpx += coeff * gx;
        dpy += coeff * gy;
      }
      // ∇W is in 1/m; convert Δp to pixels
      dxCorr[i] += dpx * PIXELS_PER_METER;
      dyCorr[i] += dpy * PIXELS_PER_METER;
    }
    for (let i = 0; i < count; i++) {
      if (held[i]) continue;
      x[i] += dxCorr[i];
      y[i] += dyCorr[i];
    }
    collideBounds();
  }

  // Rebuild velocities from position change
  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    vx[i] = (x[i] - x0[i]) * invDt;
    vy[i] = (y[i] - y0[i]) * invDt;
  }
}

/**
 * Lite DFSPH: constant-density κ solve + divergence-free velocity correction.
 * Flip SOLVER to "dfsph" to activate. Simplified 2D form of Bender/Koschier DFSPH.
 */
function computeDFSPHAlpha() {
  const mass = particleMass;
  for (let i = 0; i < count; i++) {
    if (held[i]) {
      alphaArr[i] = 0;
      continue;
    }
    const off = i * MAX_NB;
    const nnb = nbCount[i];
    const xi = x[i];
    const yi = y[i];
    let gradSumX = 0;
    let gradSumY = 0;
    let sumNorm = 0;
    for (let k = 0; k < nnb; k++) {
      const j = nbIndex[off + k];
      const r = nbDist[off + k];
      if (r < 1e-8) continue;
      const dxm = (xi - x[j]) / PIXELS_PER_METER;
      const dym = (yi - y[j]) / PIXELS_PER_METER;
      const invr = 1 / r;
      const gMag = spikyGradMag(r);
      const gx = mass * gMag * dxm * invr;
      const gy = mass * gMag * dym * invr;
      gradSumX += gx;
      gradSumY += gy;
      sumNorm += gx * gx + gy * gy;
    }
    const denom = sumNorm + gradSumX * gradSumX + gradSumY * gradSumY;
    alphaArr[i] = denom > DFSPH_EPS ? 1 / denom : 0;
  }
}

function solveDFSPHDensity(dt) {
  const rho0 = restDensity;
  const mass = particleMass;
  const invDt2 = dt > 1e-8 ? 1 / (dt * dt) : 0;
  for (let iter = 0; iter < DFSPH_DENSITY_ITERS; iter++) {
    buildHash();
    findNeighbors();
    // density
    for (let i = 0; i < count; i++) {
      let rho = mass * POLY6_SELF;
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      for (let k = 0; k < nnb; k++) rho += mass * poly6W(nbDist[off + k]);
      density[i] = rho;
    }
    computeDFSPHAlpha();
    // κ from density error (pressure-like)
    for (let i = 0; i < count; i++) {
      if (held[i]) {
        lambdaArr[i] = 0;
        continue;
      }
      const err = density[i] - rho0;
      lambdaArr[i] = err > 0 ? err * alphaArr[i] * invDt2 : 0; // κ_i
    }
    // velocity correction from pressure forces
    for (let i = 0; i < count; i++) {
      if (held[i]) continue;
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      const xi = x[i];
      const yi = y[i];
      const ki = lambdaArr[i];
      let dvx = 0;
      let dvy = 0;
      for (let k = 0; k < nnb; k++) {
        const j = nbIndex[off + k];
        const r = nbDist[off + k];
        if (r < 1e-8) continue;
        const dxm = (xi - x[j]) / PIXELS_PER_METER;
        const dym = (yi - y[j]) / PIXELS_PER_METER;
        const invr = 1 / r;
        const gMag = spikyGradMag(r);
        const gx = gMag * dxm * invr;
        const gy = gMag * dym * invr;
        const coeff = -dt * mass * (ki / Math.max(density[i], 1e-4) + lambdaArr[j] / Math.max(density[j], 1e-4));
        // convert m/s²·s style into px/s: coeff * ∇W(m) * PIXELS
        dvx += coeff * gx * PIXELS_PER_METER;
        dvy += coeff * gy * PIXELS_PER_METER;
      }
      vx[i] += dvx;
      vy[i] += dvy;
    }
  }
}

function solveDFSPHDivergence(dt) {
  const mass = particleMass;
  const invDt = dt > 1e-8 ? 1 / dt : 0;
  for (let iter = 0; iter < DFSPH_DIVERGENCE_ITERS; iter++) {
    buildHash();
    findNeighbors();
    computeDFSPHAlpha();
    // velocity divergence * (−ρ) ~ density change rate
    for (let i = 0; i < count; i++) {
      if (held[i]) {
        lambdaArr[i] = 0;
        continue;
      }
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      const xi = x[i];
      const yi = y[i];
      let div = 0;
      for (let k = 0; k < nnb; k++) {
        const j = nbIndex[off + k];
        const r = nbDist[off + k];
        if (r < 1e-8) continue;
        const dxm = (xi - x[j]) / PIXELS_PER_METER;
        const dym = (yi - y[j]) / PIXELS_PER_METER;
        const invr = 1 / r;
        const gMag = spikyGradMag(r);
        const gx = gMag * dxm * invr;
        const gy = gMag * dym * invr;
        // ∇·v ≈ −Σ m_j (v_i − v_j)·∇W  (SPH form); use (vj-vi) for density rate
        const dvx = (vx[j] - vx[i]) / PIXELS_PER_METER;
        const dvy = (vy[j] - vy[i]) / PIXELS_PER_METER;
        div += mass * (dvx * gx + dvy * gy);
      }
      // κ from divergence (want div≈0)
      lambdaArr[i] = div * alphaArr[i] * invDt;
    }
    for (let i = 0; i < count; i++) {
      if (held[i]) continue;
      const off = i * MAX_NB;
      const nnb = nbCount[i];
      const xi = x[i];
      const yi = y[i];
      const ki = lambdaArr[i];
      let dvx = 0;
      let dvy = 0;
      for (let k = 0; k < nnb; k++) {
        const j = nbIndex[off + k];
        const r = nbDist[off + k];
        if (r < 1e-8) continue;
        const dxm = (xi - x[j]) / PIXELS_PER_METER;
        const dym = (yi - y[j]) / PIXELS_PER_METER;
        const invr = 1 / r;
        const gMag = spikyGradMag(r);
        const gx = gMag * dxm * invr;
        const gy = gMag * dym * invr;
        const coeff = -dt * mass * (ki + lambdaArr[j]);
        dvx += coeff * gx * PIXELS_PER_METER;
        dvy += coeff * gy * PIXELS_PER_METER;
      }
      vx[i] += dvx;
      vy[i] += dvy;
    }
  }
}

function solveDFSPH(dt) {
  // After non-pressure forces + integrate prediction:
  // 1) divergence-free velocity, 2) constant-density, positions already advanced.
  solveDFSPHDivergence(dt);
  solveDFSPHDensity(dt);
  // Advect again with corrected v for this substep's remainder — positions already
  // moved once; re-sync x from last non-pressure step using corrected v delta is
  // approximate. Simpler: re-integrate position from stored x0 with final v.
  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    x[i] = x0[i] + vx[i] * dt;
    y[i] = y0[i] + vy[i] * dt;
  }
  collideBounds();
}

function collideParticles() {
  const iters = softCollideIters();
  for (let iter = 0; iter < iters; iter++) {
    if (iter === 0 || (iter & 1) === 0) buildHash();
    collideParticlesOnce();
  }
}


function collideBounds() {
  const r = DRAW_RADIUS;
  const minX = r;
  const maxX = width - r;
  const minY = r;
  const maxY = height - r;
  const bounce = RESTITUTION;
  const slide = 1 - WALL_FRICTION;
  for (let i = 0; i < count; i++) {
    if (held[i]) continue;
    if (x[i] < minX) {
      x[i] = minX;
      if (vx[i] < 0) vx[i] *= -bounce;
    } else if (x[i] > maxX) {
      x[i] = maxX;
      if (vx[i] > 0) vx[i] *= -bounce;
    }
    if (y[i] < minY) {
      y[i] = minY;
      if (vy[i] < 0) vy[i] *= -bounce;
    } else if (y[i] > maxY) {
      y[i] = maxY;
      if (vy[i] > 0) vy[i] *= -bounce;
      vx[i] *= slide;
    }
  }
}

function step(dt) {
  // Re-pin before forces so density/neighbors see the drop at the cursor.
  pinHeldParticles();
  buildHash();
  findNeighbors();
  computeDensityPressure(); // density for viscosity/cohesion; pressure EOS off
  applyForces(); // gravity + viscosity + depth-aware cohesion (no EOS pressure)

  // Snapshot for DFSPH position rebuild / PBF velocity rebuild
  for (let i = 0; i < count; i++) {
    x0[i] = x[i];
    y0[i] = y[i];
  }

  integrate(dt); // v += a dt, XSPH, x += v dt

  if (SOLVER === "dfsph") {
    solveDFSPH(dt);
  } else {
    // Default: Position Based Fluids density constraints
    solvePBF(dt);
  }

  collideParticles(); // light overlap cleanup
  collideBounds();
  // Final pin after collide/solver so the held blob stays glued to the mouse.
  pinHeldParticles();
}

function render() {
  ctx.fillStyle = "#07090d";
  ctx.fillRect(0, 0, width, height);

  // v1: discs only. Later, when neighbors are dense, we will render a liquid surface
  // (metaballs / marching squares). Empty pixels stay air — we do not mesh air.
  ctx.fillStyle = "rgba(64, 176, 230, 0.55)";
  ctx.beginPath();
  for (let i = 0; i < count; i++) {
    ctx.moveTo(x[i] + DRAW_RADIUS, y[i]);
    ctx.arc(x[i], y[i], DRAW_RADIUS, 0, Math.PI * 2);
  }
  ctx.fill();

  if (mouseDown || holding) {
    ctx.strokeStyle = "rgba(180, 210, 230, 0.35)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.arc(mx, my, H_PX, 0, Math.PI * 2);
    ctx.stroke();
  }

  hudN.textContent = String(count);
}

function tick(stamp) {
  requestAnimationFrame(tick);
  if (!lastStamp) lastStamp = stamp;
  const frameDt = Math.min(0.05, (stamp - lastStamp) / 1000);
  lastStamp = stamp;

  // Decay mouse velocity if the pointer hasn't moved recently (avoids stale throws).
  if (lastMouseStamp && stamp - lastMouseStamp > 50) {
    mvx = 0;
    mvy = 0;
  }

  fpsFrames++;
  if (stamp - fpsStamp >= 400) {
    const inst = (fpsFrames * 1000) / (stamp - fpsStamp);
    fpsEma = fpsEma * 0.7 + inst * 0.3;
    hudFps.textContent = String(Math.round(fpsEma));
    fpsFrames = 0;
    fpsStamp = stamp;
  }

  if (!paused) {
    accTime += frameDt;
    let steps = 0;
    while (accTime >= FIXED_DT && steps < MAX_SUBSTEPS) {
      step(FIXED_DT);
      accTime -= FIXED_DT;
      steps++;
    }
    if (steps === MAX_SUBSTEPS) accTime = 0;
  } else if (holding) {
    // Still follow the cursor while paused.
    pinHeldParticles();
  }
  render();
}

function pointerFromEvent(e) {
  const rect = canvas.getBoundingClientRect();
  const nx = e.clientX - rect.left;
  const ny = e.clientY - rect.top;
  const now = performance.now();
  if (lastMouseStamp > 0) {
    const dt = (now - lastMouseStamp) / 1000;
    if (dt > 1e-4 && dt < 0.1) {
      mvx = (nx - lastMx) / dt;
      mvy = (ny - lastMy) / dt;
    }
  }
  lastMx = nx;
  lastMy = ny;
  lastMouseStamp = now;
  mx = nx;
  my = ny;
}

function setPaused(next) {
  paused = next;
  btnPause.textContent = paused ? "Resume" : "Pause";
}

function resetSim(mode) {
  accTime = 0;
  spawnParticles(resetTarget, mode);
}

function wireUi() {
  const slN = document.getElementById("sl-n");
  const slG = document.getElementById("sl-g");
  const slRho = document.getElementById("sl-rho");
  const slK = document.getElementById("sl-k");
  const slMu = document.getElementById("sl-mu");
  const slSigma = document.getElementById("sl-sigma");
  const valN = document.getElementById("val-n");
  const valG = document.getElementById("val-g");
  const valRho = document.getElementById("val-rho");
  const valK = document.getElementById("val-k");
  const valMu = document.getElementById("val-mu");
  const valSigma = document.getElementById("val-sigma");

  slN.addEventListener("input", function () {
    resetTarget = parseInt(slN.value, 10);
    valN.textContent = String(resetTarget);
    resetSim("dam");
  });
  slG.addEventListener("input", function () {
    gravity = parseFloat(slG.value);
    valG.textContent = gravity.toFixed(2);
  });
  slRho.addEventListener("input", function () {
    restDensity = parseFloat(slRho.value);
    valRho.textContent = String(restDensity | 0);
  });
  slK.addEventListener("input", function () {
    stiffness = parseFloat(slK.value);
    valK.textContent = String(stiffness | 0);
  });
  slMu.addEventListener("input", function () {
    viscosity = parseFloat(slMu.value);
    valMu.textContent = viscosity.toFixed(2);
  });
  slSigma.addEventListener("input", function () {
    surfaceTension = parseFloat(slSigma.value);
    valSigma.textContent = surfaceTension.toFixed(1);
  });

  function syncParamsFromSliders() {
    resetTarget = parseInt(slN.value, 10);
    gravity = parseFloat(slG.value);
    restDensity = parseFloat(slRho.value);
    stiffness = parseFloat(slK.value);
    viscosity = parseFloat(slMu.value);
    surfaceTension = parseFloat(slSigma.value);
    valN.textContent = String(resetTarget);
    valG.textContent = gravity.toFixed(2);
    valRho.textContent = String(restDensity | 0);
    valK.textContent = String(stiffness | 0);
    valMu.textContent = viscosity.toFixed(2);
    valSigma.textContent = surfaceTension.toFixed(1);
  }
  syncParamsFromSliders();

  document.getElementById("btn-pause").addEventListener("click", function () {
    setPaused(!paused);
  });
  document.getElementById("btn-reset").addEventListener("click", function () {
    resetSim("dam");
  });
  document.getElementById("btn-dam").addEventListener("click", function () {
    resetSim("dam");
  });

  canvas.addEventListener("pointerdown", function (e) {
    if (e.button !== 0) return;
    pointerFromEvent(e);
    mouseDown = true;
    canvas.setPointerCapture(e.pointerId);
    // Spawn a predefined drop and hold it until release.
    if (!holding) spawnHeldDrop(mx, my);
  });
  canvas.addEventListener("pointermove", function (e) {
    pointerFromEvent(e);
  });
  canvas.addEventListener("pointerup", function () {
    releaseHeldDrop();
  });
  canvas.addEventListener("pointercancel", function () {
    releaseHeldDrop();
  });

  window.addEventListener("keydown", function (e) {
    if (e.code === "Space") {
      e.preventDefault();
      setPaused(!paused);
    } else if (e.code === "KeyR") {
      resetSim("dam");
    }
  });
  window.addEventListener("resize", setupCanvas);
}

setupCanvas();
wireUi();
resetSim("dam");
requestAnimationFrame(tick);
