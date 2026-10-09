/* 2D FLIP/APIC water — vanilla JS, no libraries, no build step.
 *
 * File layout (independent sections):
 *   0. ObstacleShapes    — pure geometry (signed distance, builders) for fixed obstacles.
 *   1. FlipSolver        — the physics core. Pure JS + typed arrays, no DOM/canvas.
 *                          Works in meters, y points DOWN, gravity is +y.
 *   2. Canvas2DRenderer  — small renderer that only READS the solver's arrays
 *                          (px, py, vx, vy, count). A WebGL renderer can replace it
 *                          by implementing the same { resize, draw } interface.
 *   3. App               — canvas sizing, HUD, sliders, mouse, Object Edit Mode, main loop
 *                          (browser only).
 *
 * Solver substep:
 *   0. dt = CFL*h/(u_max + sqrt(g h)), CFL 1, max 8 substeps per 1/60 s frame
 *   1. advect particles (x += v dt), push overlapping particles apart, clamp to walls,
 *      push particles out of fixed obstacles (analytic SDF)
 *   2. P2G   — APIC transfer with quadratic B-spline weights (mass-weighted)
 *   3. mark cells (solid = walls + rasterized obstacles / fluid / air) + density per cell
 *   4. gravity + light viscosity + damping on the grid, no-through walls
 *   5. pressure projection — PCG + MIC(0) (equal-mass faces by default; "weighted" = step 2f),
 *      plus a drift term that pushes over-dense cells apart so volume is kept
 *   6. extrapolate projected velocities 2 faces into the air
 *   7. G2P   — APIC velocity + affine matrix C, optionally blended with a FLIP delta
 */

// ===========================================================================
// 0. ObstacleShapes — pure geometry for fixed solid obstacles (no DOM)
//    Shapes are plain objects in solver units (metres, y points DOWN):
//      { type: "circle",  cx, cy, r, angle }
//      { type: "rect",    cx, cy, w, h, angle }              oriented box, full size
//      { type: "polygon", cx, cy, angle, points: [[x, y], ...], kind }
//         points are local coordinates around (cx, cy) (the area centroid), rotated
//         by angle; kind = "ramp" | "triangle" | "freehand" | "polygon" (UI label).
//    angle is in radians and turns +x toward +y (clockwise on screen, since y is down).
//    Signed distance: negative inside, positive outside; the polygon inside test is
//    even-odd, so even a self-crossing polygon has a well-defined sign.
// ===========================================================================
const ObstacleShapes = (function () {
  const isNum = (v) => typeof v === "number" && Number.isFinite(v);

  /** Signed shoelace area (> 0 = counter-clockwise in a y-up frame). */
  function polygonArea(pts) {
    let a = 0;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) a += pts[j][0] * pts[i][1] - pts[i][0] * pts[j][1];
    return 0.5 * a;
  }
  function polygonPerimeter(pts) {
    let p = 0;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) p += Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]);
    return p;
  }
  function polygonCentroid(pts) {
    let a = 0, cx = 0, cy = 0;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const cr = pts[j][0] * pts[i][1] - pts[i][0] * pts[j][1];
      a += cr;
      cx += (pts[j][0] + pts[i][0]) * cr;
      cy += (pts[j][1] + pts[i][1]) * cr;
    }
    if (Math.abs(a) < 1e-14) {
      let mx = 0, my = 0;
      for (const p of pts) { mx += p[0]; my += p[1]; }
      return [mx / pts.length, my / pts.length];
    }
    return [cx / (3 * a), cy / (3 * a)];
  }
  function orient(ax, ay, bx, by, cx, cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  }
  /** Proper or touching intersection of segments p1p2 and p3p4. */
  function segmentsIntersect(p1, p2, p3, p4) {
    const d1 = orient(p3[0], p3[1], p4[0], p4[1], p1[0], p1[1]);
    const d2 = orient(p3[0], p3[1], p4[0], p4[1], p2[0], p2[1]);
    const d3 = orient(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
    const d4 = orient(p1[0], p1[1], p2[0], p2[1], p4[0], p4[1]);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return true;
    const on = (a, b, c, d) => d === 0 && Math.min(a[0], b[0]) <= c[0] && c[0] <= Math.max(a[0], b[0]) && Math.min(a[1], b[1]) <= c[1] && c[1] <= Math.max(a[1], b[1]);
    return on(p3, p4, p1, d1) || on(p3, p4, p2, d2) || on(p1, p2, p3, d3) || on(p1, p2, p4, d4);
  }
  /** True if any two non-adjacent edges of the closed polygon intersect. */
  function selfIntersects(pts) {
    const n = pts.length;
    if (n < 4) return false;
    for (let i = 0; i < n; i++) {
      const a = pts[i], b = pts[(i + 1) % n];
      for (let j = i + 2; j < n; j++) {
        if (i === 0 && j === n - 1) continue; // these two share vertex 0
        if (segmentsIntersect(a, b, pts[j], pts[(j + 1) % n])) return true;
      }
    }
    return false;
  }
  /** Ramer-Douglas-Peucker polyline simplification (keeps the end points). */
  function rdp(pts, eps) {
    const n = pts.length;
    if (n < 3) return pts.slice();
    const keep = new Uint8Array(n);
    keep[0] = keep[n - 1] = 1;
    const stack = [[0, n - 1]];
    while (stack.length) {
      const [a, b] = stack.pop();
      const ax = pts[a][0], ay = pts[a][1], bx = pts[b][0], by = pts[b][1];
      const ex = bx - ax, ey = by - ay, ee = ex * ex + ey * ey;
      let best = -1, bi = -1;
      for (let i = a + 1; i < b; i++) {
        const wx = pts[i][0] - ax, wy = pts[i][1] - ay;
        let t = ee > 0 ? (wx * ex + wy * ey) / ee : 0;
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const d = Math.hypot(wx - ex * t, wy - ey * t);
        if (d > best) { best = d; bi = i; }
      }
      if (best > eps && bi > 0) {
        keep[bi] = 1;
        stack.push([a, bi], [bi, b]);
      }
    }
    const out = [];
    for (let i = 0; i < n; i++) if (keep[i]) out.push(pts[i]);
    return out;
  }
  /** Andrew's monotone chain convex hull. */
  function convexHull(pts) {
    const p = pts.map((q) => [q[0], q[1]]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    if (p.length < 3) return p;
    const lower = [], upper = [];
    for (const q of p) {
      while (lower.length >= 2 && orient(...lower[lower.length - 2], ...lower[lower.length - 1], q[0], q[1]) <= 0) lower.pop();
      lower.push(q);
    }
    for (let i = p.length - 1; i >= 0; i--) {
      const q = p[i];
      while (upper.length >= 2 && orient(...upper[upper.length - 2], ...upper[upper.length - 1], q[0], q[1]) <= 0) upper.pop();
      upper.push(q);
    }
    upper.pop();
    lower.pop();
    return lower.concat(upper);
  }
  function dedupe(pts, eps) {
    const out = [];
    for (const p of pts) {
      const q = out[out.length - 1];
      if (!q || Math.hypot(p[0] - q[0], p[1] - q[1]) > eps) out.push([p[0], p[1]]);
    }
    while (out.length > 1 && Math.hypot(out[0][0] - out[out.length - 1][0], out[0][1] - out[out.length - 1][1]) <= eps) out.pop();
    return out;
  }

  /** Validated clean copy of a shape. Throws TypeError on anything unusable. */
  function normalize(s) {
    if (!s || typeof s !== "object") throw new TypeError("obstacle shape must be an object");
    const angle = s.angle == null ? 0 : s.angle;
    if (!isNum(s.cx) || !isNum(s.cy) || !isNum(angle)) throw new TypeError("obstacle needs finite cx, cy, angle");
    if (s.type === "circle") {
      if (!isNum(s.r) || !(s.r > 0)) throw new TypeError("circle needs r > 0");
      return { type: "circle", cx: s.cx, cy: s.cy, r: s.r, angle };
    }
    if (s.type === "rect") {
      if (!isNum(s.w) || !isNum(s.h) || !(s.w > 0) || !(s.h > 0)) throw new TypeError("rect needs w > 0 and h > 0");
      return { type: "rect", cx: s.cx, cy: s.cy, w: s.w, h: s.h, angle };
    }
    if (s.type === "polygon") {
      if (!Array.isArray(s.points)) throw new TypeError("polygon needs points [[x, y], ...]");
      for (const p of s.points) if (!p || !isNum(p[0]) || !isNum(p[1])) throw new TypeError("polygon points must be finite [x, y]");
      const pts = dedupe(s.points, 1e-9);
      if (pts.length < 3 || Math.abs(polygonArea(pts)) < 1e-10) throw new TypeError("polygon is degenerate (needs 3+ points and non-zero area)");
      return { type: "polygon", kind: typeof s.kind === "string" ? s.kind : "polygon", cx: s.cx, cy: s.cy, angle, points: pts };
    }
    throw new TypeError("unknown obstacle type " + s.type + " (circle | rect | polygon)");
  }
  function clone(s) {
    const c = Object.assign({}, s);
    if (s.points) c.points = s.points.map((p) => [p[0], p[1]]);
    return c;
  }

  /** Precompute trig / world vertices / AABB for fast SDF queries. */
  function compile(s) {
    const c = Math.cos(s.angle || 0), sn = Math.sin(s.angle || 0);
    if (s.type === "circle") {
      return { kind: 0, cx: s.cx, cy: s.cy, r: s.r, minX: s.cx - s.r, maxX: s.cx + s.r, minY: s.cy - s.r, maxY: s.cy + s.r };
    }
    if (s.type === "rect") {
      const hw = 0.5 * s.w, hh = 0.5 * s.h;
      const ex = Math.abs(c) * hw + Math.abs(sn) * hh, ey = Math.abs(sn) * hw + Math.abs(c) * hh;
      return { kind: 1, cx: s.cx, cy: s.cy, c, s: sn, hw, hh, minX: s.cx - ex, maxX: s.cx + ex, minY: s.cy - ey, maxY: s.cy + ey };
    }
    const n = s.points.length;
    const vx = new Float64Array(n), vy = new Float64Array(n);
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (let i = 0; i < n; i++) {
      const lx = s.points[i][0], ly = s.points[i][1];
      vx[i] = s.cx + c * lx - sn * ly;
      vy[i] = s.cy + sn * lx + c * ly;
      if (vx[i] < minX) minX = vx[i];
      if (vx[i] > maxX) maxX = vx[i];
      if (vy[i] < minY) minY = vy[i];
      if (vy[i] > maxY) maxY = vy[i];
    }
    return { kind: 2, n, vx, vy, orient: polygonArea(s.points) > 0 ? 1 : -1, minX, maxX, minY, maxY };
  }

  /**
   * Signed distance of compiled shape o at (x, y); writes the outward unit normal
   * (the SDF gradient) to out.nx, out.ny when out is given.
   */
  function sdfCompiled(o, x, y, out) {
    let d, nx, ny;
    if (o.kind === 0) {
      const dx = x - o.cx, dy = y - o.cy, L = Math.hypot(dx, dy);
      d = L - o.r;
      if (L > 1e-12) { nx = dx / L; ny = dy / L; } else { nx = 0; ny = -1; }
    } else if (o.kind === 1) {
      const dx = x - o.cx, dy = y - o.cy;
      const lx = o.c * dx + o.s * dy, ly = -o.s * dx + o.c * dy;
      const qx = Math.abs(lx) - o.hw, qy = Math.abs(ly) - o.hh;
      const sx = lx < 0 ? -1 : 1, sy = ly < 0 ? -1 : 1;
      let mx, my;
      if (qx > 0 || qy > 0) {
        const ox = qx > 0 ? qx : 0, oy = qy > 0 ? qy : 0;
        d = Math.hypot(ox, oy);
        mx = (sx * ox) / d; my = (sy * oy) / d;
      } else if (qx > qy) { d = qx; mx = sx; my = 0; } else { d = qy; mx = 0; my = sy; }
      nx = o.c * mx - o.s * my;
      ny = o.s * mx + o.c * my;
    } else {
      const VX = o.vx, VY = o.vy, n = o.n;
      let best = Infinity, bqx = 0, bqy = 0, bex = 1, bey = 0, inside = false;
      for (let i = 0, j = n - 1; i < n; j = i++) {
        const ax = VX[j], ay = VY[j], ex = VX[i] - ax, ey = VY[i] - ay;
        const wx = x - ax, wy = y - ay, ee = ex * ex + ey * ey;
        let t = ee > 0 ? (wx * ex + wy * ey) / ee : 0;
        t = t < 0 ? 0 : t > 1 ? 1 : t;
        const qx = ax + ex * t, qy = ay + ey * t;
        const d2 = (x - qx) * (x - qx) + (y - qy) * (y - qy);
        if (d2 < best) { best = d2; bqx = qx; bqy = qy; bex = ex; bey = ey; }
        if ((ay > y) !== (VY[i] > y) && x < ax + ((y - ay) * ex) / ey) inside = !inside; // even-odd
      }
      const dist = Math.sqrt(best);
      d = inside ? -dist : dist;
      if (dist > 1e-12) {
        nx = (x - bqx) / dist; ny = (y - bqy) / dist;
        if (inside) { nx = -nx; ny = -ny; }
      } else {
        const L = Math.hypot(bex, bey) || 1;
        nx = (o.orient * bey) / L; ny = (-o.orient * bex) / L;
      }
    }
    if (out) { out.nx = nx; out.ny = ny; }
    return d;
  }
  function sdf(shape, x, y, out) {
    return sdfCompiled(compile(shape), x, y, out);
  }

  /** World-space outline: polygon/rect vertices (circle → null). */
  function worldVertices(s) {
    if (s.type === "circle") return null;
    const c = Math.cos(s.angle || 0), sn = Math.sin(s.angle || 0);
    const local = s.type === "rect"
      ? [[-s.w / 2, -s.h / 2], [s.w / 2, -s.h / 2], [s.w / 2, s.h / 2], [-s.w / 2, s.h / 2]]
      : s.points;
    return local.map((p) => [s.cx + c * p[0] - sn * p[1], s.cy + sn * p[0] + c * p[1]]);
  }
  /** Radius of the smallest centre-anchored circle containing the shape. */
  function boundingRadius(s) {
    if (s.type === "circle") return s.r;
    if (s.type === "rect") return 0.5 * Math.hypot(s.w, s.h);
    let r = 0;
    for (const p of s.points) r = Math.max(r, Math.hypot(p[0], p[1]));
    return r;
  }
  function translate(s, dx, dy) {
    const c = clone(s);
    c.cx += dx; c.cy += dy;
    return c;
  }
  function rotate(s, da) {
    const c = clone(s);
    let a = (c.angle || 0) + da;
    a = Math.atan2(Math.sin(a), Math.cos(a)); // keep in (-pi, pi]
    c.angle = a;
    return c;
  }

  // ---- builders (world coordinates, metres, y down) -------------------------
  function circle(cx, cy, r) { return { type: "circle", cx, cy, r, angle: 0 }; }
  function box(cx, cy, w, h, angle) { return { type: "rect", cx, cy, w, h, angle: angle || 0 }; }
  /** Polygon from world-space points; local points are stored around the area centroid. */
  function polygon(worldPts, kind, angle) {
    const [cx, cy] = polygonCentroid(worldPts);
    return { type: "polygon", kind: kind || "polygon", cx, cy, angle: angle || 0, points: worldPts.map((p) => [p[0] - cx, p[1] - cy]) };
  }
  /** Right-triangle ramp in the drag box: vertical side under the HIGHER drag point. */
  function ramp(x0, y0, x1, y1) {
    const top = Math.min(y0, y1), bot = Math.max(y0, y1);
    const hx = y0 <= y1 ? x0 : x1, lx = y0 <= y1 ? x1 : x0;
    return polygon([[hx, top], [hx, bot], [lx, bot]], "ramp");
  }
  /** Isosceles triangle in the drag box, apex up. */
  function triangle(x0, y0, x1, y1) {
    const l = Math.min(x0, x1), r = Math.max(x0, x1), top = Math.min(y0, y1), bot = Math.max(y0, y1);
    return polygon([[0.5 * (l + r), top], [r, bot], [l, bot]], "triangle");
  }
  /**
   * Freehand stroke (world points) → closed, simplified polygon.
   * RDP with eps = 0.5 cell (then coarser until <= 64 vertices). Rejects strokes that
   * enclose almost no area or are thinner than ~1.5 cells on average (they would
   * not even show up on the grid). A stroke that crosses itself is replaced by its
   * convex hull (note set) so the solid is always a simple polygon.
   * Returns { shape } or { error }.
   */
  function freehand(worldPts, cell) {
    let pts = dedupe(worldPts, 1e-6);
    if (pts.length < 3) return { error: "stroke too short: hold and draw a closed loop" };
    let eps = 0.5 * cell;
    let simp = rdp(pts.concat([pts[0]]), eps);
    while (simp.length > 65 && eps < 20 * cell) { eps *= 1.5; simp = rdp(pts.concat([pts[0]]), eps); }
    simp = dedupe(simp, 1e-6);
    if (simp.length < 3) return { error: "stroke is a line: draw a closed loop" };
    let note = "";
    if (selfIntersects(simp)) {
      simp = convexHull(simp);
      note = "stroke crossed itself: using its outline (convex hull)";
    }
    const area = Math.abs(polygonArea(simp));
    const perim = polygonPerimeter(simp);
    if (simp.length < 3 || area < 4 * cell * cell) return { error: "shape too small (needs a few grid cells of area)" };
    if ((2 * area) / perim < 0.75 * cell) return { error: "shape too thin: draw a fatter closed loop" };
    const shape = polygon(simp, "freehand");
    return note ? { shape, note } : { shape };
  }

  /**
   * Editor rule (item 1): minimum feature thickness of 2 grid cells (2h).
   * Two tests on the solver's own raster (an interior cell is solid iff its centre has
   * sdf < 0); tank walls and everything outside the tank count as solid (a part lying on a
   * wall merges with it), other obstacles do not (they can be removed later):
   *   1. medial thickness (SDF): a solid cell centre passes if it lies within h/2 of a disc
   *      of radius h that fits inside shape ∪ walls (morphological opening by a 2h-wide
   *      disc, searched on an h/4 lattice). A strip thinner than 2h has no such disc at all;
   *      the h/2 tolerance lets 90° corners pass and fails corners sharper than ~84°.
   *   2. grid check: no two adjacent solid cells form a one-cell-thick row or column (both
   *      side neighbours non-solid) - the configuration whose ghost face sees fluid on both
   *      sides (the energy blow-up fixed in fillGhosts).
   * Returns { ok, cells, thin: [[i, j], ...] (failing cells), minCells: 2 }.
   */
  const MIN_FEATURE_CELLS = 2;
  function thicknessCheck(shape, nx, ny, h) {
    const o = compile(normalize(shape)), ih = 1 / h;
    const R = 0.5 * MIN_FEATURE_CELLS * h, tol = 0.5 * h, step = 0.25 * h;
    const x0 = h, x1 = (nx - 1) * h, y0 = h, y1 = (ny - 1) * h; // tank interior
    const wallSdf = (x, y) => { // > 0 inside the tank, < 0 inside walls / beyond
      const dx = Math.min(x - x0, x1 - x), dy = Math.min(y - y0, y1 - y);
      if (dx >= 0 && dy >= 0) return Math.min(dx, dy);
      const ox = Math.max(0, -dx), oy = Math.max(0, -dy);
      return -Math.hypot(ox, oy) || Math.max(dx, dy);
    };
    const i0 = Math.max(1, Math.floor(o.minX * ih - 0.5)), i1 = Math.min(nx - 2, Math.ceil(o.maxX * ih - 0.5));
    const j0 = Math.max(1, Math.floor(o.minY * ih - 0.5)), j1 = Math.min(ny - 2, Math.ceil(o.maxY * ih - 0.5));
    const res = { ok: false, cells: 0, thin: [], minCells: MIN_FEATURE_CELLS };
    if (i0 > i1 || j0 > j1) return res;
    const W = i1 - i0 + 3, Hh = j1 - j0 + 3, solid = new Uint8Array(W * Hh), mine = [];
    for (let j = j0 - 1; j <= j1 + 1; j++) for (let i = i0 - 1; i <= i1 + 1; i++) {
      const k = (i - i0 + 1) + (j - j0 + 1) * W;
      if (i <= 0 || j <= 0 || i >= nx - 1 || j >= ny - 1) solid[k] = 1; // tank wall
      else if (i >= i0 && i <= i1 && j >= j0 && j <= j1 && sdfCompiled(o, (i + 0.5) * h, (j + 0.5) * h, null) < 0) { solid[k] = 1; mine.push(k); }
    }
    res.cells = mine.length;
    // signed distance of shape ∪ walls (min of the two; conservative at the junction)
    const uSdf = (x, y) => Math.min(sdfCompiled(o, x, y, null), wallSdf(x, y));
    const n = Math.ceil((R + tol) / step), rr = ((R + tol) / step) ** 2;
    const bad = new Uint8Array(W * Hh);
    for (const k of mine) {
      const cx = ((k % W) + i0 - 1 + 0.5) * h, cy = (((k / W) | 0) + j0 - 1 + 0.5) * h;
      let ok = uSdf(cx, cy) <= -R;
      for (let b = -n; b <= n && !ok; b++) for (let a = -n; a <= n && !ok; a++) {
        if (a * a + b * b > rr) continue;
        ok = uSdf(cx + a * step, cy + b * step) <= -R;
      }
      if (!ok) bad[k] = 1;
      // one-cell-thick column (vertical pair with both side columns open) / row
      if (solid[k + W] && !solid[k - 1] && !solid[k - 1 + W] && !solid[k + 1] && !solid[k + 1 + W]) { bad[k] = 1; bad[k + W] = 1; }
      if (solid[k + 1] && !solid[k - W] && !solid[k - W + 1] && !solid[k + W] && !solid[k + W + 1]) { bad[k] = 1; bad[k + 1] = 1; }
    }
    for (const k of mine) if (bad[k]) res.thin.push([(k % W) + i0 - 1, ((k / W) | 0) + j0 - 1]);
    res.ok = res.cells > 0 && res.thin.length === 0;
    return res;
  }
  /**
   * Clip the convex polygon corners sharper than maxDeg (default 88 deg) so the cut across
   * each tip is `cut` long (at most 45 % of either edge). The clipped polygon replaces the
   * shape itself, so the exact particle collision shape and the rasterized grid solid agree
   * (clipping only the raster would leave an unrasterized sliver that water passes through).
   * Returns { shape, clipped } (clipped = number of corners cut; non-polygons unchanged).
   */
  function clipAcuteTips(shape, cut, maxDeg) {
    const s = normalize(shape);
    if (s.type !== "polygon") return { shape: s, clipped: 0 };
    const lim = ((maxDeg || 88) * Math.PI) / 180;
    const V = worldVertices(s), n = V.length, sg = polygonArea(V) > 0 ? 1 : -1;
    const out = [];
    let clipped = 0;
    for (let k = 0; k < n; k++) {
      const p = V[(k + n - 1) % n], v = V[k], q = V[(k + 1) % n];
      const e1x = p[0] - v[0], e1y = p[1] - v[1], e2x = q[0] - v[0], e2y = q[1] - v[1];
      const L1 = Math.hypot(e1x, e1y), L2 = Math.hypot(e2x, e2y);
      const cr = (v[0] - p[0]) * (q[1] - v[1]) - (v[1] - p[1]) * (q[0] - v[0]);
      const ang = Math.acos(Math.max(-1, Math.min(1, (e1x * e2x + e1y * e2y) / (L1 * L2 || 1))));
      if (L1 > 0 && L2 > 0 && cr * sg > 0 && ang < lim) {
        const d = (0.5 * cut) / Math.sin(0.5 * ang), s1 = Math.min(d, 0.45 * L1), s2 = Math.min(d, 0.45 * L2);
        out.push([v[0] + (e1x / L1) * s1, v[1] + (e1y / L1) * s1], [v[0] + (e2x / L2) * s2, v[1] + (e2y / L2) * s2]);
        clipped++;
      } else out.push([v[0], v[1]]);
    }
    return clipped ? { shape: normalize(polygon(out, s.kind)), clipped } : { shape: s, clipped: 0 };
  }

  return {
    normalize, clone, compile, sdfCompiled, sdf, worldVertices, boundingRadius, translate, rotate,
    circle, box, polygon, ramp, triangle, freehand,
    rdp, convexHull, selfIntersects, polygonArea, polygonCentroid,
    MIN_FEATURE_CELLS, thicknessCheck, clipAcuteTips,
  };
})();

/**
 * Demo obstacle set used by index.html?demoObstacles=1 and test/obstacles.js:
 * a circle in the dam-break path and a ramp on the floor that launches the wave.
 * Positions scale with the tank (solver must be init()ed).
 */
function demoObstacleShapes(solver) {
  const h = solver.cellSize;
  const W = (solver.nx - 2) * h, H = (solver.ny - 2) * h;
  const floor = (solver.ny - 1) * h;
  return [
    ObstacleShapes.circle(h + 0.5 * W, floor - 0.2 * H, 0.04 * W),
    ObstacleShapes.ramp(h + 0.88 * W, floor - 0.13 * H, h + 0.68 * W, floor),
  ];
}

// ===========================================================================
// 1. FlipSolver — physics core (no DOM access)
// ===========================================================================
const FlipSolver = (function () {
  // ---- tunables -----------------------------------------------------------
  const PARTICLE_RADIUS = 0.3; // in cells -> ~3.2 particles per cell (hex packing)
  const SEPARATION_DIST = 1.6; // min particle spacing, in particle radii (2.0 = touching)
  const SEPARATION_ITERS = 1; // push-apart passes per substep
  const FIXED_DT = 1 / 60; // one simulation frame
  const MAX_SUBSTEPS = 8; // cap per frame; if hit, the frame is logged and the sim falls behind
  const CFL = 1.0; // dt = CFL * h / (u_max + sqrt(g h)); never exceeded
  const OVER_RELAX = 1.9; // SOR factor for Gauss-Seidel (gs solver only)
  const PCG_REL_TOL = 1e-6; // stop when ||r||2 / ||b||2 <= 1e-6 ...
  const PCG_ABS_TOL = 1e-4; // ... or max|r| (as divergence) * dt <= 1e-4
  const MIC_TAU = 0.97; // MIC(0) modification parameter (Bridson)
  const MIC_SIGMA = 0.25; // MIC(0) safety: fall back to the diagonal if the pivot gets small
  const DRIFT_RATE = 1.0; // 1/s [2a]: softened from 4.0; how fast over-dense cells relax back to rest density
  const WALL_FRICTION = 0.0; // legacy boundary only: 0 = free-slip, 1 = no-slip (mirror walls are free-slip)
  const MAX_SPEED = 14; // m/s safety clamp
  const MAX_THROW = 8; // m/s clamp on mouse-throw release velocity
  const DAM_HEIGHT_FRAC = 0.58; // dam-break column height / tank height
  const DROP_RADIUS_CELLS = 5; // held drop radius in cells (~250 particles)

  const SOLID = 0;
  const FLUID = 1;
  const AIR = 2;

  // Slider-controlled parameters (defaults match index.html).
  const params = {
    fill: 0.22, // water amount: fraction of the tank area in the dam column
    gravity: 9.81, // m/s^2
    flipRatio: 0.8, // 0 = pure APIC, 1 = APIC transfer + full FLIP delta (0.8: lively, still quiet at rest)
    // Pressure solver: "pcg" = conjugate gradient with MIC(0) preconditioner, stops on
    // the residual (default); "gs" = the old red-black Gauss-Seidel/SOR (comparison).
    pressureSolver: "pcg",
    pressureIters: 50, // Gauss-Seidel iterations per substep (gs only)
    pcgMaxIters: 500, // PCG safety cap per solve (logged if ever hit)
    // Particle push-apart pass: "neutral" = every pair push is equal and opposite AND
    // limited (per axis) to what both particles can move without leaving the wall
    // bounds, so the pass plus the wall clamp cannot move the centre of mass (default).
    // Step 2e: when the floor blocks a vertical push for a floor-row pair, the lost
    // separation is redirected into x (equal-opposite, re-clipped at the side walls).
    // "legacy" = the old pass (equal/opposite, but the wall clamp afterwards can undo
    // one side of a push, which lifts the fluid). Comparison only.
    separation: "neutral",
    // Box-wall handling:
    //   "separating" (default, step 2g) = free-slip walls that can push but not pull
    //     (Signorini / active-set). Starts from the 2d mirror walls, then releases any
    //     wall face whose adjacent fluid cell has p < 0 and puts it back if the fluid
    //     later moves into the wall. Cap 5 passes per substep.
    //   "mirror" (step 2d) = fixed free-slip walls, normal = 0 both ways. Comparison.
    //   "legacy" = old one-sided walls + ghost faces keeping their own P2G mass.
    boundary: "separating",
    // Pressure face masses (step 2f):
    //   "equal" (default until still-pool passes) = unit-face-mass projection.
    //   "weighted" = mass-weighted projection. Face mass m_f = sum of the same P2G
    //     weights used for momentum; rho_f = m_f / V_f with V_f = h^2; solve
    //     D·(dt/rho_f)·G p = D u* with rho_f clamped to >= 0.1*rho0. Still SPD →
    //     PCG + MIC(0), preconditioner rebuilt every substep.
    pressureMass: "equal",
    // [2a] Drift correction rate (1/s). Softened default 1.0 (was 4.0). Set 4 for the
    // pre-2a baseline. 0 disables drift.
    driftRate: 1.0,
    // [2b] Cell-density sampler: "renorm" (default) skips SOLID corners and renormalizes;
    // "legacy" keeps the old splat (Validation grades density on the legacy sampler when
    // both are posted).
    densitySampler: "renorm",
    viscosity: 0.05, // 0..1 grid velocity diffusion per substep
    // 1/s light linear damping: stands in for boundary-layer losses the grid cannot
    // resolve, so the browser pool settles in ~30 s. Set to 0 for validation runs.
    damping: 0.15,
  };

  // ---- grid state (nx, ny include one solid cell on every side) ------------
  let nx = 0;
  let ny = 0;
  let h = 1;
  let invH = 1;
  let pr = 0.3; // particle radius (m)
  let restDensity = 0; // particles per cell at rest, measured right after reset
  let u, v, uOld, vOld, uW, vW, uTmp, vTmp; // MAC faces: u (nx+1)*ny, v nx*(ny+1)
  let uScale, vScale; // [2f] face scale = rho0/rho_f = (dt/rho_f)·G coeff in velocity-p units
  let cellType, cellDensity, cellDensityLegacy; // nx*ny; legacy = pre-2b sampler for Validation grading
  let uValid, vValid, uValidNext, vValidNext; // face masks for velocity extrapolation
  let fluidList, fluidKu, fluidFlags, fluidDrift; // packed fluid cells for the solver
  let cellDrift; // drift target per cell (same values as fluidDrift)
  // PCG work arrays (per cell, float64) and the 5-point matrix
  let pcgCells, pcgP, pcgR, pcgZ, pcgS, pcgQ, pcgPrecon, pcgDiag, pcgPlusI, pcgPlusJ;
  let pcgCount = 0;
  const lastPcg = { iters: 0, relRes: 0, maxResDt: 0, capHit: false };
  const pcgCapLog = []; // solves that hit pcgMaxIters
  // [2f] face-density clamp stats for the mass-weighted projection (per solve + accum)
  const lastClamp = { hits: 0, faces: 0, fraction: 0 };
  const clampAccum = { hits: 0, faces: 0 };
  // Separating-wall active set (2g): per fluid cell, bits 1=L 2=R 4=T(ceiling) 8=B(floor)
  // mark a solid-neighbour face as released (Dirichlet p=0, face free). Cleared each solve.
  let wallRelease; // Uint8Array nc
  let uSnap, vSnap; // pre-pressure velocity snapshots for active-set restarts
  const lastSep = {
    passes: 0, capHit: false, negPLeft: false,
    worstNegP: 0, worstUInto: 0, compFail: 0,
  };
  const sepPassLog = []; // per-substep pass counts (capped length)
  let sepPassSum = 0, sepPassMax = 0, sepCapHits = 0, sepSubsteps = 0;
  // centre-of-mass shift caused by the last separation pass (free particles, metres;
  // y is down): raw = right after the pass; eff = after the wall clamp is applied to the
  // positions before and after the pass, i.e. what the pass leaves once walls act.
  const sepCom = { rawDx: 0, rawDy: 0, effDx: 0, effDy: 0, n: 0 };
  let pcgWarnFrame = -1e9;
  let fluidRed = 0; // fluidList[0..fluidRed) red cells, then black cells
  let fluidCount = 0;

  // ---- particles (SoA typed arrays) ---------------------------------------
  let maxN = 0;
  let count = 0;
  let px, py, pvx, pvy, c00, c01, c10, c11, held, holdOx, holdOy;

  // ---- separation hash ----------------------------------------------------
  let hashCs = 1;
  let hashNx = 1;
  let hashNy = 1;
  let hashStart, hashIds, hashCellOf;

  // ---- held drop / pointer (world units) ----------------------------------
  let holding = false;
  let pointerX = 0;
  let pointerY = 0;
  let pointerVx = 0;
  let pointerVy = 0;
  let lastSubsteps = 1;
  let frame = 0; // frames stepped since init
  let simTime = 0; // simulated seconds since reset (falls behind real time on cap hits)
  let lastUmax = 0;
  const capLog = []; // { frame, t, umax, behind } for frames that hit MAX_SUBSTEPS
  let capWarnFrame = -1e9;
  let substepHook = null; // optional diagnostics callback (headless validation)
  let nanRecoveries = 0; // particles re-placed after a NaN (the only reseed path)
  let maxCflObserved = 0; // max u_max*dt/h over all substeps since init

  let rngState = 12345;
  function rand() {
    // tiny deterministic LCG so headless runs are repeatable
    rngState = (rngState * 1664525 + 1013904223) >>> 0;
    return rngState / 4294967296;
  }

  /** Allocate a tank of cellsX x cellsY interior cells of size cellSize (m). */
  function init(cellsX, cellsY, cellSize) {
    nx = cellsX + 2;
    ny = cellsY + 2;
    h = cellSize;
    invH = 1 / h;
    pr = PARTICLE_RADIUS * h;
    const nu = (nx + 1) * ny;
    const nv = nx * (ny + 1);
    const nc = nx * ny;
    u = new Float32Array(nu);
    uOld = new Float32Array(nu);
    uW = new Float32Array(nu);
    uTmp = new Float32Array(nu);
    uScale = new Float32Array(nu);
    v = new Float32Array(nv);
    vOld = new Float32Array(nv);
    vW = new Float32Array(nv);
    vTmp = new Float32Array(nv);
    vScale = new Float32Array(nv);
    uValid = new Uint8Array(nu);
    uValidNext = new Uint8Array(nu);
    vValid = new Uint8Array(nv);
    vValidNext = new Uint8Array(nv);
    cellType = new Uint8Array(nc);
    cellDensity = new Float32Array(nc);
    cellDensityLegacy = new Float32Array(nc);
    fluidList = new Int32Array(nc);
    fluidKu = new Int32Array(nc);
    fluidFlags = new Uint8Array(nc);
    fluidDrift = new Float32Array(nc);
    cellDrift = new Float64Array(nc);
    pcgCells = new Int32Array(nc);
    pcgP = new Float64Array(nc);
    pcgR = new Float64Array(nc);
    pcgZ = new Float64Array(nc);
    pcgS = new Float64Array(nc);
    pcgQ = new Float64Array(nc);
    pcgPrecon = new Float64Array(nc);
    pcgDiag = new Float64Array(nc);
    pcgPlusI = new Float64Array(nc);
    pcgPlusJ = new Float64Array(nc);
    wallRelease = new Uint8Array(nc);
    uSnap = new Float32Array(nu);
    vSnap = new Float32Array(nv);
    obsMask = new Uint8Array(nc);
    baseType = new Uint8Array(nc);
    nearMask = new Uint8Array(nc);

    // room for a ~70% full tank (dam + many drops)
    const hexPerCell = (h * h) / (4 * pr * pr * 0.8660254);
    maxN = Math.ceil(cellsX * cellsY * hexPerCell * 0.7) + 1000;
    px = new Float32Array(maxN);
    py = new Float32Array(maxN);
    pvx = new Float32Array(maxN);
    pvy = new Float32Array(maxN);
    c00 = new Float32Array(maxN);
    c01 = new Float32Array(maxN);
    c10 = new Float32Array(maxN);
    c11 = new Float32Array(maxN);
    held = new Uint8Array(maxN);
    holdOx = new Float32Array(maxN);
    holdOy = new Float32Array(maxN);

    hashCs = SEPARATION_DIST * pr;
    hashNx = Math.ceil((nx * h) / hashCs) + 1;
    hashNy = Math.ceil((ny * h) / hashCs) + 1;
    hashStart = new Int32Array(hashNx * hashNy + 1);
    hashIds = new Int32Array(maxN);
    hashCellOf = new Int32Array(maxN);

    count = 0;
    holding = false;
    restDensity = 0;
    rebuildObstacleGrid(); // obstacles (metres) survive a re-init / resize
  }

  // ---- spawning -----------------------------------------------------------
  function addParticle(x, y, vx, vy) {
    if (count >= maxN) return -1;
    const i = count++;
    px[i] = x;
    py[i] = y;
    pvx[i] = vx;
    pvy[i] = vy;
    c00[i] = c01[i] = c10[i] = c11[i] = 0;
    held[i] = 0;
    return i;
  }

  function resetObstacleCounters() {
    obstacleDeleted = 0; lastObstacleDeleted = 0; dropClipped = 0; pushStuck = 0;
  }
  /** Spawn helper for reset/resetColumn: never place a particle inside or touching an obstacle. */
  function addParticleFree(x, y) {
    if (obsActive && obstacleSdfAt(x, y, pr, null) < pr) return -1;
    return addParticle(x, y, 0, 0);
  }

  /**
   * "dam": hex-packed column on the left (DAM_HEIGHT_FRAC of tank height). "pool": flat layer.
   * Obstacles are kept (no particles spawn inside them) unless opts.clearObstacles.
   */
  function reset(mode, opts) {
    if (opts && opts.clearObstacles) clearObstacles();
    resetObstacleCounters();
    count = 0;
    holding = false;
    simTime = 0;
    restDensity = 0; // re-measured on the first substep
    rngState = 12345;
    sepPassLog.length = 0; sepPassSum = 0; sepPassMax = 0; sepCapHits = 0; sepSubsteps = 0;
    clampAccum.hits = 0; clampAccum.faces = 0;
    lastClamp.hits = 0; lastClamp.faces = 0; lastClamp.fraction = 0;
    if (wallRelease) wallRelease.fill(0);
    const innerW = (nx - 2) * h;
    const innerH = (ny - 2) * h;
    const dx = 2 * pr;
    const dy = dx * 0.8660254;
    const pool = mode === "pool";
    const colH = innerH * (pool ? params.fill : DAM_HEIGHT_FRAC);
    const colW = pool ? innerW : Math.min(innerW * 0.92, (params.fill * innerW) / DAM_HEIGHT_FRAC);
    const x0 = h + pr;
    const yFloor = (ny - 1) * h - pr;
    const rows = Math.max(1, Math.floor(colH / dy));
    const cols = Math.max(1, Math.floor((colW - pr) / dx));
    for (let r = 0; r < rows; r++) {
      const ox = (r & 1) * pr;
      for (let c = 0; c < cols; c++) {
        const jx = (rand() - 0.5) * 0.02 * h;
        addParticleFree(x0 + c * dx + ox + jx, yFloor - r * dy);
      }
    }
  }

  /**
   * Exact validation column: cellsA x cellsH cells against the left wall on a dry
   * bed, 2x2 particles per cell on a regular lattice (centres at 1/4 and 3/4 of
   * each cell), so the column edges coincide with cell edges and the particle
   * volume is exactly cellsA*cellsH*h^2. Particles start at rest. Obstacles are
   * kept; lattice sites inside (or touching) one are skipped.
   */
  function resetColumn(cellsA, cellsH) {
    resetObstacleCounters();
    count = 0;
    holding = false;
    simTime = 0;
    restDensity = 0;
    rngState = 12345; // same seed as reset(): coincident-particle jitter must not leak across cases
    sepPassLog.length = 0; sepPassSum = 0; sepPassMax = 0; sepCapHits = 0; sepSubsteps = 0;
    clampAccum.hits = 0; clampAccum.faces = 0;
    lastClamp.hits = 0; lastClamp.faces = 0; lastClamp.fraction = 0;
    if (wallRelease) wallRelease.fill(0);
    for (let j = 0; j < cellsH; j++) {
      for (let i = 0; i < cellsA; i++) {
        const x0 = (1 + i) * h; // interior starts after one wall cell
        const y0 = (ny - 2 - j) * h; // row j counted up from the floor cell
        addParticleFree(x0 + 0.25 * h, y0 + 0.25 * h);
        addParticleFree(x0 + 0.75 * h, y0 + 0.25 * h);
        addParticleFree(x0 + 0.25 * h, y0 + 0.75 * h);
        addParticleFree(x0 + 0.75 * h, y0 + 0.75 * h);
      }
    }
  }

  function clampDropCenter(c, lo, hi) {
    return Math.min(Math.max(c, lo), hi);
  }

  /** Spawn a round hex-packed drop at (cx, cy) and pin it to the pointer until release. */
  function spawnDrop(cx, cy) {
    if (holding) return;
    const R = DROP_RADIUS_CELLS * h;
    const dx = 2 * pr;
    const dy = dx * 0.8660254;
    pointerX = cx;
    pointerY = cy;
    const n = Math.ceil(R / dy);
    const m = Math.ceil(R / dx) + 1;
    let any = false;
    for (let r = -n; r <= n; r++) {
      const oy = r * dy;
      for (let c = -m; c <= m; c++) {
        const ox = c * dx + (r & 1) * pr;
        if (ox * ox + oy * oy > R * R) continue;
        const i = addParticle(cx, cy, 0, 0);
        if (i < 0) break;
        held[i] = 1;
        holdOx[i] = ox;
        holdOy[i] = oy;
        any = true;
      }
    }
    holding = any;
    pinHeld();
  }

  /** Free the held drop with the (clamped) pointer velocity. */
  function releaseDrop() {
    if (!holding) return;
    // drop particles released inside an obstacle would be shot out by the push-out: drop them
    if (obsActive) dropClipped += compactParticles((i) => !(held[i] && obstacleSdfAt(px[i], py[i], pr, null) < 0));
    let tvx = pointerVx;
    let tvy = pointerVy;
    const sp = Math.hypot(tvx, tvy);
    if (sp > MAX_THROW) {
      tvx *= MAX_THROW / sp;
      tvy *= MAX_THROW / sp;
    }
    for (let i = 0; i < count; i++) {
      if (!held[i]) continue;
      held[i] = 0;
      pvx[i] = tvx;
      pvy[i] = tvy;
      c00[i] = c01[i] = c10[i] = c11[i] = 0;
    }
    holding = false;
  }

  /** Pointer position (m) and velocity (m/s) in world units. */
  function setPointer(x, y, vx, vy) {
    pointerX = x;
    pointerY = y;
    pointerVx = vx;
    pointerVy = vy;
  }

  /** Glue the held drop to the pointer (kept fully inside the tank). */
  function pinHeld() {
    if (!holding) return;
    const R = DROP_RADIUS_CELLS * h;
    const cx = clampDropCenter(pointerX, h + R + pr, (nx - 1) * h - R - pr);
    const cy = clampDropCenter(pointerY, h + R + pr, (ny - 1) * h - R - pr);
    for (let i = 0; i < count; i++) {
      if (!held[i]) continue;
      px[i] = cx + holdOx[i];
      py[i] = cy + holdOy[i];
      pvx[i] = pointerVx;
      pvy[i] = pointerVy;
    }
  }

  // ---- particle advection, separation, walls ------------------------------
  function advectParticles(dt) {
    const X = px, Y = py, VX = pvx, VY = pvy, HD = held;
    for (let i = 0; i < count; i++) {
      if (HD[i]) continue;
      X[i] += VX[i] * dt;
      Y[i] += VY[i] * dt;
    }
  }

  /** Gauss-Seidel push-apart so particles keep a minimum spacing (no clumping). */
  function separateParticles(iters) {
    const X = px, Y = py, HD = held;
    const start = hashStart, ids = hashIds, cellOf = hashCellOf;
    const hnx = hashNx, hny = hashNy;
    const minD = SEPARATION_DIST * pr;
    const minD2 = minD * minD;
    const invCs = 1 / hashCs;
    const nCells = hnx * hny;
    const n = count;
    const neutral = params.separation !== "legacy";
    // the same bounds collideWalls() clamps to
    const bx0 = h + pr, bx1 = (nx - 1) * h - pr, by0 = h + pr, by1 = (ny - 1) * h - pr;
    const clampX = (x) => (x < bx0 ? bx0 : x > bx1 ? bx1 : x);
    const clampY = (y) => (y < by0 ? by0 : y > by1 ? by1 : y);
    let rx0 = 0, ry0 = 0, ex0 = 0, ey0 = 0, nFree = 0;
    for (let i = 0; i < n; i++) {
      if (HD[i]) continue;
      rx0 += X[i]; ry0 += Y[i]; ex0 += clampX(X[i]); ey0 += clampY(Y[i]); nFree++;
    }
    // counting sort into hash cells
    start.fill(0);
    for (let i = 0; i < n; i++) {
      let cx = (X[i] * invCs) | 0;
      let cy = (Y[i] * invCs) | 0;
      if (cx < 0) cx = 0;
      else if (cx >= hnx) cx = hnx - 1;
      if (cy < 0) cy = 0;
      else if (cy >= hny) cy = hny - 1;
      const c = cx + cy * hnx;
      cellOf[i] = c;
      start[c]++;
    }
    let acc = 0;
    for (let c = 0; c < nCells; c++) {
      acc += start[c];
      start[c] = acc;
    }
    start[nCells] = acc;
    for (let i = 0; i < n; i++) ids[--start[cellOf[i]]] = i;

    for (let it = 0; it < iters; it++) {
      for (let i = 0; i < n; i++) {
        const c = cellOf[i];
        const cy = (c / hnx) | 0;
        const cx = c - cy * hnx;
        const x0 = cx > 0 ? cx - 1 : 0;
        const x1 = cx < hnx - 1 ? cx + 1 : cx;
        const y0 = cy > 0 ? cy - 1 : 0;
        const y1 = cy < hny - 1 ? cy + 1 : cy;
        const hi = HD[i];
        for (let yy = y0; yy <= y1; yy++) {
          const rowEndBase = yy * hnx;
          const kEnd = start[rowEndBase + x1 + 1];
          // cells x0..x1 of one row are contiguous in the sorted id list
          for (let k = start[rowEndBase + x0]; k < kEnd; k++) {
            const j = ids[k];
            if (j <= i) continue;
            const hj = HD[j];
            if (hi && hj) continue;
            let dx = X[j] - X[i];
            let dy = Y[j] - Y[i];
            const d2 = dx * dx + dy * dy;
            if (d2 >= minD2) continue;
            if (d2 < 1e-18) {
              const a = rand() * 6.2831853;
              dx = Math.cos(a) * 0.5 * minD;
              dy = Math.sin(a) * 0.5 * minD;
            } else {
              const d = Math.sqrt(d2);
              const s = (0.5 * (minD - d)) / d;
              dx *= s;
              dy *= s;
            }
            // a held (pinned) particle is immovable; the free one takes the full push
            // (the user's hand is an external force, so this is not momentum-neutral)
            if (neutral && !hi && !hj) {
              // i moves by -(dx,dy), j by +(dx,dy). Per axis, shrink the shared push to
              // what BOTH can move inside [lo, hi] (a particle already outside the bounds
              // may move inward but not further out), so it stays equal and opposite.
              const dx0 = dx, dy0 = dy;
              if (dx > 0) {
                const li = X[i] - (X[i] < bx0 ? X[i] : bx0), lj = (X[j] > bx1 ? X[j] : bx1) - X[j];
                if (dx > li) dx = li;
                if (dx > lj) dx = lj;
              } else if (dx < 0) {
                const li = (X[i] > bx1 ? X[i] : bx1) - X[i], lj = X[j] - (X[j] < bx0 ? X[j] : bx0);
                if (-dx > li) dx = -li;
                if (-dx > lj) dx = -lj;
              }
              if (dy > 0) {
                const li = Y[i] - (Y[i] < by0 ? Y[i] : by0), lj = (Y[j] > by1 ? Y[j] : by1) - Y[j];
                if (dy > li) dy = li;
                if (dy > lj) dy = lj;
              } else if (dy < 0) {
                const li = (Y[i] > by1 ? Y[i] : by1) - Y[i], lj = Y[j] - (Y[j] < by0 ? Y[j] : by0);
                if (-dy > li) dy = -li;
                if (-dy > lj) dy = -lj;
              }
              // [2e] Floor fully blocked a vertical push (|dy|≈0 after clip but |dy0|>0):
              // replace the lost vertical separation with an equal-opposite x push, then
              // re-clip at the side walls. Magnitude capped to original |δ|. Floor band =
              // one particle radius above the floor bound.
              const floorBand = by1 - pr;
              if (Y[i] >= floorBand && Y[j] >= floorBand) {
                const absDy0 = dy0 < 0 ? -dy0 : dy0;
                const absDy = dy < 0 ? -dy : dy;
                // only when the floor ate essentially all of the vertical half-gap
                if (absDy0 > 1e-15 && absDy <= 1e-15 * absDy0) {
                  const add = (X[j] > X[i] || (X[j] === X[i] && (i & 1) === 0)) ? absDy0 : -absDy0;
                  dx = dx0 + add;
                  dy = 0;
                  if (dx > 0) {
                    const li = X[i] - (X[i] < bx0 ? X[i] : bx0), lj = (X[j] > bx1 ? X[j] : bx1) - X[j];
                    if (dx > li) dx = li;
                    if (dx > lj) dx = lj;
                  } else if (dx < 0) {
                    const li = (X[i] > bx1 ? X[i] : bx1) - X[i], lj = X[j] - (X[j] < bx0 ? X[j] : bx0);
                    if (-dx > li) dx = -li;
                    if (-dx > lj) dx = -lj;
                  }
                  const mag02 = dx0 * dx0 + dy0 * dy0;
                  const mag2 = dx * dx;
                  if (mag2 > mag02 && mag2 > 0) dx *= Math.sqrt(mag02 / mag2);
                }
              }
              if (obsActive) { clipSepPush(X[i], Y[i], X[j], Y[j], dx, dy); dx = sepClipOut.x; dy = sepClipOut.y; }
              X[i] -= dx;
              Y[i] -= dy;
              X[j] += dx;
              Y[j] += dy;
            } else if (hi) {
              X[j] += 2 * dx;
              Y[j] += 2 * dy;
            } else if (hj) {
              X[i] -= 2 * dx;
              Y[i] -= 2 * dy;
            } else {
              X[i] -= dx;
              Y[i] -= dy;
              X[j] += dx;
              Y[j] += dy;
            }
          }
        }
      }
    }
    let rx1 = 0, ry1 = 0, ex1 = 0, ey1 = 0;
    for (let i = 0; i < n; i++) {
      if (HD[i]) continue;
      rx1 += X[i]; ry1 += Y[i]; ex1 += clampX(X[i]); ey1 += clampY(Y[i]);
    }
    const inv = nFree > 0 ? 1 / nFree : 0;
    sepCom.rawDx = (rx1 - rx0) * inv; sepCom.rawDy = (ry1 - ry0) * inv;
    sepCom.effDx = (ex1 - ex0) * inv; sepCom.effDy = (ey1 - ey0) * inv;
    sepCom.n = nFree;
  }

  // [item 2 EXPERIMENT, not in sim.js] obstacle-aware separation: the shared, equal-and-opposite
  // push of a pair is clipped along each particle's obstacle normal so that neither centre is
  // pushed closer than r + m to an obstacle (the same idea as the wall clipping above).
  const sepClipOut = { x: 0, y: 0 }, sepN = { nx: 0, ny: 0 };
  function nearObs(x, y) { let ci = (x * invH) | 0, cj = (y * invH) | 0; if (ci < 0) ci = 0; else if (ci >= nx) ci = nx - 1; if (cj < 0) cj = 0; else if (cj >= ny) cj = ny - 1; return nearMask[ci + cj * nx]; }
  function clipSepPush(xi, yi, xj, yj, dx, dy) {
    const rc = pr + APPROACH_MARGIN_CELLS * h;
    if (nearObs(xi, yi)) { const d = obstacleSdfAt(xi, yi, pr + h, sepN), A = Math.max(d - rc, 0), dn = dx * sepN.nx + dy * sepN.ny; if (dn > A) { dx -= (dn - A) * sepN.nx; dy -= (dn - A) * sepN.ny; } }
    if (nearObs(xj, yj)) { const d = obstacleSdfAt(xj, yj, pr + h, sepN), A = Math.max(d - rc, 0), dn = dx * sepN.nx + dy * sepN.ny; if (dn < -A) { dx += (-A - dn) * sepN.nx; dy += (-A - dn) * sepN.ny; } }
    sepClipOut.x = dx; sepClipOut.y = dy;
  }

  /** Keep particles inside the box; remove velocity into the wall. NaN-safe. */
  function collideWalls() {
    const X = px, Y = py, VX = pvx, VY = pvy, HD = held;
    const minX = h + pr;
    const maxX = (nx - 1) * h - pr;
    const minY = h + pr;
    const maxY = (ny - 1) * h - pr;
    for (let i = 0; i < count; i++) {
      if (HD[i]) continue;
      let x = X[i];
      let y = Y[i];
      if (x !== x || y !== y || VX[i] !== VX[i] || VY[i] !== VY[i]) {
        // should never happen; recover instead of poisoning the grid (counted)
        nanRecoveries++;
        x = 0.5 * (minX + maxX);
        y = minY;
        VX[i] = VY[i] = 0;
        c00[i] = c01[i] = c10[i] = c11[i] = 0;
      }
      if (x < minX) {
        x = minX;
        if (VX[i] < 0) VX[i] = 0;
      } else if (x > maxX) {
        x = maxX;
        if (VX[i] > 0) VX[i] = 0;
      }
      if (y < minY) {
        y = minY;
        if (VY[i] < 0) VY[i] = 0;
      } else if (y > maxY) {
        y = maxY;
        if (VY[i] > 0) VY[i] = 0;
      }
      X[i] = x;
      Y[i] = y;
    }
  }

  // ---- grid transfers (APIC, quadratic B-spline) ---------------------------
  // For node row starting at `base`, with d = f - base in [0.5, 1.5):
  //   w0 = 0.5(1.5-d)^2, w1 = 0.75-(d-1)^2, w2 = 0.5(d-0.5)^2
  // The APIC inertia tensor for quadratic B-splines is D = h^2/4 * I.

  /** Scatter one velocity component (faces offset by (ox0, oy0) cells) with APIC. */
  function p2gComponent(G, GW, stride, offX, offY, VEL, CA, CB) {
    const X = px, Y = py, HD = held;
    const hh = h, ih = invH;
    for (let p = 0; p < count; p++) {
      if (HD[p]) continue;
      const fx = X[p] * ih - offX;
      const fy = Y[p] * ih - offY;
      const bx = Math.floor(fx - 0.5);
      const by = Math.floor(fy - 0.5);
      const dx = fx - bx;
      const dy = fy - by;
      const wx0 = 0.5 * (1.5 - dx) * (1.5 - dx);
      const wx1 = 0.75 - (dx - 1) * (dx - 1);
      const wx2 = 0.5 * (dx - 0.5) * (dx - 0.5);
      const wy0 = 0.5 * (1.5 - dy) * (1.5 - dy);
      const wy1 = 0.75 - (dy - 1) * (dy - 1);
      const wy2 = 0.5 * (dy - 0.5) * (dy - 0.5);
      const vp = VEL[p];
      const ca = CA[p] * hh;
      const cb = CB[p] * hh;
      // node offsets (in cells) from particle: ox = a - dx, oy = b - dy
      const ax0 = vp + ca * (0 - dx), ax1 = vp + ca * (1 - dx), ax2 = vp + ca * (2 - dx);
      let k = by * stride + bx;
      for (let b = 0; b < 3; b++, k += stride) {
        const wy = b === 0 ? wy0 : b === 1 ? wy1 : wy2;
        const cy = cb * (b - dy);
        let w = wx0 * wy;
        GW[k] += w;
        G[k] += w * (ax0 + cy);
        w = wx1 * wy;
        GW[k + 1] += w;
        G[k + 1] += w * (ax1 + cy);
        w = wx2 * wy;
        GW[k + 2] += w;
        G[k + 2] += w * (ax2 + cy);
      }
    }
  }

  function particlesToGrid() {
    u.fill(0);
    v.fill(0);
    uW.fill(0);
    vW.fill(0);
    p2gComponent(u, uW, nx + 1, 0, 0.5, pvx, c00, c01); // u faces at (i h, (j+.5) h)
    p2gComponent(v, vW, nx, 0.5, 0, pvy, c10, c11); // v faces at ((i+.5) h, j h)
    const mirror = params.boundary !== "legacy";
    if (mirror) foldGhostFaces();
    if (obsActive) clearObstacleGhostMass();
    for (let i = 0; i < u.length; i++) if (uW[i] > 0) u[i] /= uW[i];
    for (let i = 0; i < v.length; i++) if (vW[i] > 0) v[i] /= vW[i];
    // mirror walls: apply the wall conditions before saving the FLIP reference, so wall
    // faces never show up as a FLIP delta (no interior normal velocity reaches G2P)
    if (mirror) enforceBoundary();
    uOld.set(u);
    vOld.set(v);
  }

  /**
   * [mirror / separating walls] Fold ghost-face mass+momentum onto the fluid-side
   * image BEFORE dividing by weight (dig 2 / 2f):
   *   - Tangential ghosts (u in ceiling/floor rows, v in side columns) → free-slip image.
   *   - Normal ghosts (deep u beyond left/right walls, deep v beyond ceiling/floor) →
   *     the interface normal face. Keeps m_f / rho_f used by the mass-weighted projection
   *     consistent with total P2G mass; released and constrained faces then share that
   *     same folded mass (no orphan mass on zeroed ghosts).
   */
  function foldGhostFaces() {
    const su = nx + 1, top = 0, bot = (ny - 1) * su;
    // Tangential: ceiling/floor rows of u
    for (let i = 0; i <= nx; i++) {
      u[top + su + i] += u[top + i]; uW[top + su + i] += uW[top + i]; u[top + i] = 0; uW[top + i] = 0;
      u[bot - su + i] += u[bot + i]; uW[bot - su + i] += uW[bot + i]; u[bot + i] = 0; uW[bot + i] = 0;
    }
    // Tangential: side columns of v
    for (let j = 0; j <= ny; j++) {
      const r = j * nx;
      v[r + 1] += v[r]; vW[r + 1] += vW[r]; v[r] = 0; vW[r] = 0;
      v[r + nx - 2] += v[r + nx - 1]; vW[r + nx - 2] += vW[r + nx - 1]; v[r + nx - 1] = 0; vW[r + nx - 1] = 0;
    }
    // Normal: deep left/right u ghosts → interface faces (i=1 and i=nx-1)
    for (let j = 0; j < ny; j++) {
      const r = j * su;
      u[r + 1] += u[r]; uW[r + 1] += uW[r]; u[r] = 0; uW[r] = 0;
      u[r + nx - 1] += u[r + nx]; uW[r + nx - 1] += uW[r + nx]; u[r + nx] = 0; uW[r + nx] = 0;
    }
    // Normal: deep ceiling/floor v ghosts → interface faces (j=1 and j=ny-1)
    const c1 = ny * nx;
    for (let i = 0; i < nx; i++) {
      v[i + nx] += v[i]; vW[i + nx] += vW[i]; v[i] = 0; vW[i] = 0;
      v[i + (ny - 1) * nx] += v[i + c1]; vW[i + (ny - 1) * nx] += vW[i + c1]; v[i + c1] = 0; vW[i + c1] = 0;
    }
  }

  /** Gather one component: APIC velocity + C row, blended with the FLIP delta. */
  function g2pComponent(G, GOLD, stride, offX, offY, VEL, CA, CB) {
    const X = px, Y = py, HD = held;
    const ih = invH;
    const flip = params.flipRatio;
    const pic = 1 - flip;
    const cScale = 4 * ih; // D^-1 * h (offsets are in cells): 4/h^2 * h
    for (let p = 0; p < count; p++) {
      if (HD[p]) continue;
      const fx = X[p] * ih - offX;
      const fy = Y[p] * ih - offY;
      const bx = Math.floor(fx - 0.5);
      const by = Math.floor(fy - 0.5);
      const dx = fx - bx;
      const dy = fy - by;
      const wx0 = 0.5 * (1.5 - dx) * (1.5 - dx);
      const wx1 = 0.75 - (dx - 1) * (dx - 1);
      const wx2 = 0.5 * (dx - 0.5) * (dx - 0.5);
      const wy0 = 0.5 * (1.5 - dy) * (1.5 - dy);
      const wy1 = 0.75 - (dy - 1) * (dy - 1);
      const wy2 = 0.5 * (dy - 0.5) * (dy - 0.5);
      let vel = 0;
      let dv = 0;
      let cx = 0;
      let cy = 0;
      let k = by * stride + bx;
      for (let b = 0; b < 3; b++, k += stride) {
        const wy = b === 0 ? wy0 : b === 1 ? wy1 : wy2;
        const oy = b - dy;
        const g0 = G[k], g1 = G[k + 1], g2 = G[k + 2];
        const w0 = wx0 * wy, w1 = wx1 * wy, w2 = wx2 * wy;
        const s = w0 * g0 + w1 * g1 + w2 * g2;
        vel += s;
        dv += w0 * (g0 - GOLD[k]) + w1 * (g1 - GOLD[k + 1]) + w2 * (g2 - GOLD[k + 2]);
        cx += w0 * g0 * -dx + w1 * g1 * (1 - dx) + w2 * g2 * (2 - dx);
        cy += s * oy;
      }
      VEL[p] = pic * vel + flip * (VEL[p] + dv);
      CA[p] = cx * cScale;
      CB[p] = cy * cScale;
    }
  }

  function gridToParticles() {
    g2pComponent(u, uOld, nx + 1, 0, 0.5, pvx, c00, c01);
    g2pComponent(v, vOld, nx, 0.5, 0, pvy, c10, c11);
    const VX = pvx, VY = pvy;
    const m2 = MAX_SPEED * MAX_SPEED;
    for (let p = 0; p < count; p++) {
      const s2 = VX[p] * VX[p] + VY[p] * VY[p];
      if (s2 > m2) {
        const s = MAX_SPEED / Math.sqrt(s2);
        VX[p] *= s;
        VY[p] *= s;
      }
    }
  }

  // ---- grid: cells, density, forces, boundary, pressure --------------------
  function markCells() {
    const T = cellType, D = cellDensity, X = px, Y = py, HD = held;
    T.set(baseType); // box walls + rasterized obstacles are SOLID, the rest starts as AIR
    D.fill(0);
    cellDensityLegacy.fill(0);
    const ih = invH;
    for (let p = 0; p < count; p++) {
      if (HD[p]) continue;
      const x = X[p] * ih;
      const y = Y[p] * ih;
      const c = (x | 0) + (y | 0) * nx;
      if (T[c] === AIR) T[c] = FLUID;
      // bilinear splat to cell centres -> particle density per cell.
      // [2b] Always write the legacy splat into cellDensityLegacy (Validation grades on
      // it). Physics uses cellDensity: "renorm" (default) drops SOLID corners and
      // renormalizes; "legacy" copies the old splat.
      const fx = x - 0.5;
      const fy = y - 0.5;
      const i0 = fx | 0; // fx > 0 inside the tank
      const j0 = fy | 0;
      const tx = fx - i0;
      const ty = fy - j0;
      const b = i0 + j0 * nx;
      const w00 = (1 - tx) * (1 - ty);
      const w10 = tx * (1 - ty);
      const w01 = (1 - tx) * ty;
      const w11 = tx * ty;
      const DL = cellDensityLegacy;
      DL[b] += w00; DL[b + 1] += w10; DL[b + nx] += w01; DL[b + nx + 1] += w11;
      if (params.densitySampler === "legacy") {
        D[b] += w00; D[b + 1] += w10; D[b + nx] += w01; D[b + nx + 1] += w11;
      } else {
        let r00 = T[b] === SOLID ? 0 : w00;
        let r10 = T[b + 1] === SOLID ? 0 : w10;
        let r01 = T[b + nx] === SOLID ? 0 : w01;
        let r11 = T[b + nx + 1] === SOLID ? 0 : w11;
        const ws = r00 + r10 + r01 + r11;
        if (ws > 0 && ws < 1 - 1e-12) {
          const inv = 1 / ws;
          r00 *= inv; r10 *= inv; r01 *= inv; r11 *= inv;
        }
        D[b] += r00; D[b + 1] += r10; D[b + nx] += r01; D[b + nx + 1] += r11;
      }
    }
    if (restDensity === 0) measureRestDensity();
  }

  /** Rest density = mean density of fluid cells with no air neighbour, right after reset. */
  function measureRestDensity() {
    let sum = 0;
    let n = 0;
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx;
        if (cellType[c] !== FLUID) continue;
        if (cellType[c - 1] === AIR || cellType[c + 1] === AIR) continue;
        if (cellType[c - nx] === AIR || cellType[c + nx] === AIR) continue;
        sum += cellDensity[c];
        n++;
      }
    }
    restDensity = n > 0 ? sum / n : (h * h) / (4 * pr * pr * 0.8660254);
  }

  function applyGravity(dt) {
    const dv = params.gravity * dt;
    for (let i = 0; i < v.length; i++) if (vW[i] > 0) v[i] += dv;
    if (params.damping > 0) {
      const k = Math.max(0, 1 - params.damping * dt);
      for (let i = 0; i < u.length; i++) u[i] *= k;
      for (let i = 0; i < v.length; i++) v[i] *= k;
    }
  }

  /** Light viscosity: one damped Jacobi diffusion pass over faces that carry fluid. */
  function diffuse(G, GW, T, stride, i0, i1, j0, j1, a) {
    T.set(G);
    for (let j = j0; j < j1; j++) {
      for (let i = i0; i < i1; i++) {
        const k = i + j * stride;
        if (GW[k] <= 0) continue;
        let s = 0;
        let n = 0;
        if (GW[k - 1] > 0) { s += T[k - 1]; n++; }
        if (GW[k + 1] > 0) { s += T[k + 1]; n++; }
        if (GW[k - stride] > 0) { s += T[k - stride]; n++; }
        if (GW[k + stride] > 0) { s += T[k + stride]; n++; }
        if (n > 0) G[k] += a * (s / n - T[k]);
      }
    }
  }
  function applyViscosity() {
    const a = params.viscosity * 0.5;
    if (a <= 0) return;
    diffuse(u, uW, uTmp, nx + 1, 1, nx, 1, ny - 1, a);
    diffuse(v, vW, vTmp, nx, 1, nx - 1, 1, ny, a);
  }

  /**
   * Box walls. "separating" (default) / "mirror" / "legacy": normal faces on and beyond the walls are zero, and
   * tangential faces inside the wall cells copy their mirror (fluid-side) face; those
   * faces carry no mass (see foldGhostFaces), so the copy creates no momentum.
   * "legacy": only flow INTO a wall is removed, and the ghost faces keep their own P2G
   * mass. Either way the quadratic stencil of particles next to a wall never reads
   * stale, unprojected values from inside the solid.
   */
  function enforceBoundary() {
    enforceWalls();
    if (obsActive) enforceObstacles(); // obstacle faces follow the same rule as wall faces
  }
  function enforceWalls() {
    const su = nx + 1;
    if (params.boundary === "separating") {
      // [2g] free-slip separating walls: deep ghost normals always 0; interface normals
      // zero only while constrained; released faces keep their post-solve velocity.
      for (let j = 0; j < ny; j++) {
        const r = j * su;
        u[r] = 0; u[r + nx] = 0; // deep ghosts
        const cL = 1 + j * nx, cR = (nx - 2) + j * nx;
        if (j >= 1 && j <= ny - 2) {
          if (!(wallRelease[cL] & 1)) u[r + 1] = 0;
          if (!(wallRelease[cR] & 2)) u[r + nx - 1] = 0;
        } else {
          u[r + 1] = 0; u[r + nx - 1] = 0;
        }
      }
      const c0 = (ny - 1) * nx, c1 = ny * nx;
      for (let i = 0; i < nx; i++) {
        v[i] = 0; v[i + c1] = 0; // deep ghosts
        if (i >= 1 && i <= nx - 2) {
          const cT = i + 1 * nx, cB = i + (ny - 2) * nx;
          if (!(wallRelease[cT] & 4)) v[i + nx] = 0;
          if (!(wallRelease[cB] & 8)) v[i + c0] = 0;
        } else {
          v[i + nx] = 0; v[i + c0] = 0;
        }
      }
      const top = 0, bot = (ny - 1) * su;
      for (let i = 0; i <= nx; i++) {
        u[top + i] = u[top + su + i];
        u[bot + i] = u[bot - su + i];
      }
      for (let j = 0; j <= ny; j++) {
        const r = j * nx;
        v[r] = v[r + 1];
        v[r + nx - 1] = v[r + nx - 2];
      }
      return;
    }
    if (params.boundary !== "legacy") {
      // [mirror / 2d] fixed free-slip walls: normal faces on and beyond every wall are zero
      // (both directions); tangential faces inside the wall copy their mirror face
      for (let j = 0; j < ny; j++) {
        const r = j * su;
        u[r] = 0; u[r + 1] = 0; u[r + nx - 1] = 0; u[r + nx] = 0;
      }
      const c0 = (ny - 1) * nx, c1 = ny * nx;
      for (let i = 0; i < nx; i++) {
        v[i] = 0; v[i + nx] = 0; v[i + c0] = 0; v[i + c1] = 0;
      }
      const top = 0, bot = (ny - 1) * su;
      for (let i = 0; i <= nx; i++) {
        u[top + i] = u[top + su + i];
        u[bot + i] = u[bot - su + i];
      }
      for (let j = 0; j <= ny; j++) {
        const r = j * nx;
        v[r] = v[r + 1];
        v[r + nx - 1] = v[r + nx - 2];
      }
      return;
    }
    // [legacy] One-sided walls: flow INTO a wall is removed, flow AWAY from it is kept,
    // so water can peel off the ceiling / side walls instead of sticking to them.
    for (let j = 0; j < ny; j++) {
      const r = j * su;
      if (u[r] < 0) u[r] = 0; // left wall: only +x allowed
      if (u[r + 1] < 0) u[r + 1] = 0;
      if (u[r + nx - 1] > 0) u[r + nx - 1] = 0; // right wall: only -x
      if (u[r + nx] > 0) u[r + nx] = 0;
    }
    const b0 = (ny - 1) * nx;
    const b1 = ny * nx;
    for (let i = 0; i < nx; i++) {
      if (v[i] < 0) v[i] = 0; // ceiling: only +y (down) allowed
      if (v[i + nx] < 0) v[i + nx] = 0;
      if (v[i + b0] > 0) v[i + b0] = 0; // floor: only -y (up)
      if (v[i + b1] > 0) v[i + b1] = 0;
    }
    const slip = 1 - WALL_FRICTION;
    const top = 0, bot = (ny - 1) * su;
    for (let i = 0; i <= nx; i++) {
      u[top + i] = slip * u[top + su + i];
      u[bot + i] = slip * u[bot - su + i];
    }
    for (let j = 0; j <= ny; j++) {
      const r = j * nx;
      v[r] = slip * v[r + 1];
      v[r + nx - 1] = slip * v[r + nx - 2];
    }
  }

  /** Pack fluid cells (red first, then black) with solid-neighbour flags + drift target. */
  function buildFluidList() {
    const T = cellType, D = cellDensity;
    const rest = restDensity;
    const driftRate = params.driftRate != null ? params.driftRate : DRIFT_RATE;
    const driftScale = rest > 0 ? (driftRate * h) / rest : 0;
    cellDrift.fill(0);
    let n = 0;
    for (let color = 0; color < 2; color++) {
      for (let j = 1; j < ny - 1; j++) {
        for (let i = 1 + ((j + color + 1) & 1); i < nx - 1; i += 2) {
          const c = i + j * nx;
          if (T[c] !== FLUID) continue;
          const f =
            (T[c - 1] !== SOLID ? 1 : 0) |
            (T[c + 1] !== SOLID ? 2 : 0) |
            (T[c - nx] !== SOLID ? 4 : 0) |
            (T[c + nx] !== SOLID ? 8 : 0);
          if (f === 0) continue;
          fluidList[n] = c;
          fluidKu[n] = c + j; // u index of the cell's left face
          fluidFlags[n] = f;
          // Over-dense cells push apart. Under-dense cells are only pulled together
          // when fully submerged (no air neighbour), so the free surface never clumps.
          let comp = D[c] - rest;
          if (comp < 0 && (T[c - 1] === AIR || T[c + 1] === AIR || T[c - nx] === AIR || T[c + nx] === AIR)) comp = 0;
          fluidDrift[n] = driftScale * comp;
          cellDrift[c] = fluidDrift[n];
          n++;
        }
      }
      if (color === 0) fluidRed = n;
    }
    fluidCount = n;
  }

  /** Pressure projection: PCG (default) or the old Gauss-Seidel, per params.pressureSolver. */
  function solvePressure(dt) {
    buildFluidList();
    if (params.pressureSolver === "gs") {
      wallRelease.fill(0);
      solvePressureGS();
      lastSep.passes = 1; lastSep.capHit = false; lastSep.negPLeft = false;
      lastSep.worstNegP = 0; lastSep.worstUInto = 0; lastSep.compFail = 0;
      return;
    }
    if (params.boundary !== "separating") {
      wallRelease.fill(0);
      solvePressurePCG(dt);
      lastSep.passes = 1; lastSep.capHit = false; lastSep.negPLeft = false;
      lastSep.worstNegP = 0; lastSep.worstUInto = 0; lastSep.compFail = 0;
      return;
    }
    solvePressureSeparating(dt);
  }

  /**
   * [2g] Separating (Signorini) wall contact on top of the 2d mirror walls.
   * A wall can push but cannot pull. Active set, NOT a post-solve pressure clamp:
   *   1. Solve with all wall faces constrained (normal = 0).
   *   2. Release any wall face whose adjacent fluid cell has p < 0; re-solve.
   *   3. Put back any released face the fluid is now moving into; re-solve.
   * Cap 5 passes. Hitting the cap with a negative wall pressure left is a fail.
   */
  function solvePressureSeparating(dt) {
    wallRelease.fill(0);
    uSnap.set(u);
    vSnap.set(v);
    // NM active set, interleaved, release-priority on conflict (wall cannot pull):
    //   into > 0 && p >= 0  -> constrain (put back)
    //   p < 0               -> release
    // Cap 5. Never clamp pressure.
    let passes = 0;
    let settled = false;
    const nextRel = new Uint8Array(wallRelease.length);
    const su = nx + 1;
    while (passes < 5) {
      passes++;
      u.set(uSnap); v.set(vSnap); enforceBoundary(); solvePressurePCG(dt);
      nextRel.fill(0);
      let changed = false;
      for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx; if (cellType[c] !== FLUID) continue;
        const rel = wallRelease[c], p = pcgP[c]; let nrel = 0;
        const face = (bit, into) => {
          const was = (rel & bit) !== 0;
          if (p < 0) { nrel |= bit; if (!was) changed = true; }
          else if (into > 0) { if (was) changed = true; /* constrain */ }
          else if (was) nrel |= bit;
        };
        if (cellType[c - 1] === SOLID) face(1, -u[i + j * su]);
        if (cellType[c + 1] === SOLID) face(2, u[(i + 1) + j * su]);
        if (cellType[c - nx] === SOLID) face(4, -v[c]);
        if (cellType[c + nx] === SOLID) face(8, v[c + nx]);
        nextRel[c] = nrel;
      }
      if (!changed) { settled = true; wallRelease.set(nextRel); break; }
      if (passes < 5) wallRelease.set(nextRel);
      else break; // keep wallRelease used by last solve for the check
    }
    // Residual put-back without snap-restore: re-constrain released faces that still
    // show into > 0, zero them, and project once from the current field (keeps the
    // already-solved interior; avoids the chatter of a full snap restart).
    {
      let pb = false;
      for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx; if (cellType[c] !== FLUID) continue;
        let rel = wallRelease[c];
        if (cellType[c - 1] === SOLID && (rel & 1) && -u[i + j * su] > 0) { wallRelease[c] &= ~1; u[i + j * su] = 0; pb = true; }
        if (cellType[c + 1] === SOLID && (rel & 2) && u[(i + 1) + j * su] > 0) { wallRelease[c] &= ~2; u[(i + 1) + j * su] = 0; pb = true; }
        if (cellType[c - nx] === SOLID && (rel & 4) && -v[c] > 0) { wallRelease[c] &= ~4; v[c] = 0; pb = true; }
        if (cellType[c + nx] === SOLID && (rel & 8) && v[c + nx] > 0) { wallRelease[c] &= ~8; v[c + nx] = 0; pb = true; }
      }
      if (pb) { solvePressurePCG(dt); if (passes < 5) passes++; }
      // If put-back reintroduced tension, release those faces once more from the current field.
      let rel2 = false;
      for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx; if (cellType[c] !== FLUID) continue;
        const p = pcgP[c], rel = wallRelease[c];
        if (cellType[c - 1] === SOLID && !(rel & 1) && p < 0) { wallRelease[c] |= 1; rel2 = true; }
        if (cellType[c + 1] === SOLID && !(rel & 2) && p < 0) { wallRelease[c] |= 2; rel2 = true; }
        if (cellType[c - nx] === SOLID && !(rel & 4) && p < 0) { wallRelease[c] |= 4; rel2 = true; }
        if (cellType[c + nx] === SOLID && !(rel & 8) && p < 0) { wallRelease[c] |= 8; rel2 = true; }
      }
      if (rel2) { solvePressurePCG(dt); if (passes < 5) passes++; }
    }
    let worstNegP = 0, worstUInto = 0, compFail = 0, negPLeft = false;
    for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
      const c = i + j * nx; if (cellType[c] !== FLUID) continue;
      const rel = wallRelease[c], pCell = pcgP[c];
      const check = (bit, into) => {
        const released = (rel & bit) !== 0; const p = released ? 0 : pCell;
        if (p < worstNegP) worstNegP = p; if (into > worstUInto) worstUInto = into;
        if (!released && pCell < 0) negPLeft = true;
        if (!(p >= -1e-12 && into <= 1e-6 && (Math.abs(p) <= 1e-6 || Math.abs(into) <= 1e-6))) compFail++;
      };
      if (cellType[c - 1] === SOLID) check(1, -u[i + j * su]);
      if (cellType[c + 1] === SOLID) check(2, u[(i + 1) + j * su]);
      if (cellType[c - nx] === SOLID) check(4, -v[c]);
      if (cellType[c + nx] === SOLID) check(8, v[c + nx]);
    }
    lastSep.passes = passes; lastSep.capHit = !settled; lastSep.negPLeft = negPLeft;
    lastSep.worstNegP = worstNegP; lastSep.worstUInto = worstUInto; lastSep.compFail = compFail;
    sepSubsteps++; sepPassSum += passes; if (passes > sepPassMax) sepPassMax = passes;
    if (!settled) sepCapHits++;
    if (sepPassLog.length < 20000) sepPassLog.push(passes);
  }

  /**
   * [gs] Velocity-form Gauss-Seidel (Ten Minute Physics style), red-black ordered
   * SOR, fixed iteration count. Air cells are p = 0, solid faces stay fixed. The
   * drift term asks over-dense cells for a small positive divergence.
   */
  function solvePressureGS() {
    const U = u, V = v, L = fluidList, KU = fluidKu, F = fluidFlags, DR = fluidDrift;
    const stride = nx;
    const n = fluidCount;
    const omega = OVER_RELAX;
    const iters = params.pressureIters;
    for (let it = 0; it < iters; it++) {
      for (let q = 0; q < n; q++) {
        const c = L[q];
        const ku = KU[q];
        const f = F[q];
        const sl = f & 1;
        const sr = (f >> 1) & 1;
        const st = (f >> 2) & 1;
        const sb = (f >> 3) & 1;
        const div = U[ku + 1] - U[ku] + V[c + stride] - V[c] - DR[q];
        const p = (-div / (sl + sr + st + sb)) * omega;
        U[ku] -= sl * p;
        U[ku + 1] += sr * p;
        V[c] -= st * p;
        V[c + stride] += sb * p;
      }
    }
  }

  /**
   * [pcg] Conjugate gradient + MIC(0) (Bridson). Air p = 0.
   * Equal-mass (params.pressureMass = "equal"): A_cc = # non-solid neighbours,
   *   A_cn = -1, face update u -= p_right - p_left (velocity units).
   * Mass-weighted (params.pressureMass = "weighted", step 2f): face mass m_f from
   *   P2G weights (uW/vW), V_f = h^2, rho_f = m_f/V_f clamped to >= 0.1*rho0.
   *   NM: solve D·(dt/rho_f)·G p = D u*, then u -= (dt/rho_f)·G p (D volumetric).
   *   In velocity-unit p (p absorbs dt/(rho0 h)): scale_f = rho0/rho_f on BOTH the
   *   Laplacian and the face update (same array, same clamp). Still SPD → PCG+MIC(0).
   * Stops on ||r||/||b|| <= 1e-6 OR max|r|/h * dt <= 1e-4; params.pcgMaxIters cap.
   */
  function solvePressurePCG(dt) {
    const T = cellType, U = u, V = v, su = nx + 1;
    const C = pcgCells, P = pcgP, R = pcgR, Z = pcgZ, S = pcgS, Q = pcgQ;
    const M = pcgPrecon, Ad = pcgDiag, Ai = pcgPlusI, Aj = pcgPlusJ, DR = cellDrift;
    P.fill(0); R.fill(0); Z.fill(0); S.fill(0); Q.fill(0); M.fill(0); Ad.fill(0); Ai.fill(0); Aj.fill(0);
    // fluid cells in natural order + matrix + right-hand side.
    // Released wall faces (2g) are treated as air: count in the diagonal (Dirichlet
    // p=0) and update the face after the solve.
    const WR = wallRelease;
    const weighted = params.pressureMass !== "equal";
    // Dig 1+3: shared face scale; rho_f = m_f/V_f with V_f = h^2 (not m_f alone,
    // not a second volume factor). rho0 = restDensity/V_f. Clamp rho_f >= 0.1*rho0.
    // scale_f = rho0/rho_f used identically in A and in u -= scale_f*Δp (velocity-p units).
    const Vf = h * h; // face dual volume in 2D
    const restM = restDensity > 0 ? restDensity : 1; // particle mass per cell at rest = rho0*Vf
    const rho0 = restM / Vf;
    const rhoMin = 0.1 * rho0;
    const mMin = rhoMin * Vf; // == 0.1*restM
    let cHits = 0, cFaces = 0;
    if (weighted) {
      uValid.fill(0); vValid.fill(0);
      for (let k = 0; k < uW.length; k++) {
        const mf = uW[k];
        const rhoF = mf / Vf; // dig 3: density from mass/volume, once
        const rho = rhoF >= rhoMin ? rhoF : rhoMin;
        uScale[k] = rho0 / rho; // == restM/max(mf,mMin); no extra Vf
      }
      for (let k = 0; k < vW.length; k++) {
        const mf = vW[k];
        const rhoF = mf / Vf;
        const rho = rhoF >= rhoMin ? rhoF : rhoMin;
        vScale[k] = rho0 / rho;
      }
    }
    const markU = (idx) => {
      if (uValid[idx]) return uScale[idx];
      uValid[idx] = 1; cFaces++;
      if (uW[idx] < mMin) cHits++;
      return uScale[idx];
    };
    const markV = (idx) => {
      if (vValid[idx]) return vScale[idx];
      vValid[idx] = 1; cFaces++;
      if (vW[idx] < mMin) cHits++;
      return vScale[idx];
    };
    const sigU = (idx) => (weighted ? markU(idx) : 1);
    const sigV = (idx) => (weighted ? markV(idx) : 1);
    let m = 0;
    let bNorm2 = 0;
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx;
        if (T[c] !== FLUID) continue;
        const rel = WR[c];
        const L = T[c - 1] !== SOLID || (rel & 1);
        const Rgt = T[c + 1] !== SOLID || (rel & 2);
        const Tp = T[c - nx] !== SOLID || (rel & 4);
        const Bt = T[c + nx] !== SOLID || (rel & 8);
        const ku = c + j;
        let d;
        if (!weighted) {
          d = (L ? 1 : 0) + (Rgt ? 1 : 0) + (Tp ? 1 : 0) + (Bt ? 1 : 0);
          if (d === 0) continue;
          C[m++] = c;
          Ad[c] = d;
          if (T[c + 1] === FLUID) Ai[c] = -1;
          if (T[c + nx] === FLUID) Aj[c] = -1;
        } else {
          const sL = L ? sigU(ku) : 0;
          const sR = Rgt ? sigU(ku + 1) : 0;
          const sT = Tp ? sigV(c) : 0;
          const sB = Bt ? sigV(c + nx) : 0;
          d = sL + sR + sT + sB;
          if (d === 0) continue;
          C[m++] = c;
          Ad[c] = d;
          if (T[c + 1] === FLUID) Ai[c] = -uScale[ku + 1]; // same scale as face update
          if (T[c + nx] === FLUID) Aj[c] = -vScale[c + nx];
        }
        const b = DR[c] - (U[ku + 1] - U[ku] + V[c + nx] - V[c]);
        R[c] = b;
        bNorm2 += b * b;
      }
    }
    lastClamp.hits = cHits; lastClamp.faces = cFaces;
    lastClamp.fraction = cFaces ? cHits / cFaces : 0;
    clampAccum.hits += cHits; clampAccum.faces += cFaces;
    pcgCount = m;
    const resScale = dt / h; // |r| (= h * divergence) -> divergence * dt
    lastPcg.iters = 0;
    lastPcg.capHit = false;
    let maxR = 0;
    for (let k = 0; k < m; k++) maxR = Math.max(maxR, Math.abs(R[C[k]]));
    lastPcg.relRes = bNorm2 > 0 ? 1 : 0;
    lastPcg.maxResDt = maxR * resScale;
    if (m === 0 || bNorm2 === 0 || maxR * resScale <= PCG_ABS_TOL) return;
    const bNorm = Math.sqrt(bNorm2);

    // MIC(0) factor
    for (let k = 0; k < m; k++) {
      const c = C[k];
      const a = Ad[c];
      const ml = M[c - 1], mt = M[c - nx];
      const ail = Ai[c - 1] * ml; // coupling to the left neighbour (0 if not fluid)
      const ajt = Aj[c - nx] * mt; // coupling to the neighbour above
      let e = a - ail * ail - ajt * ajt - MIC_TAU * (Ai[c - 1] * Aj[c - 1] * ml * ml + Aj[c - nx] * Ai[c - nx] * mt * mt);
      if (e < MIC_SIGMA * a) e = a;
      M[c] = 1 / Math.sqrt(e);
    }
    function applyPrecon() {
      // solve L q = r (forward), then L^T z = q (backward); Q is scratch here
      for (let k = 0; k < m; k++) {
        const c = C[k];
        const t = R[c] - Ai[c - 1] * M[c - 1] * Q[c - 1] - Aj[c - nx] * M[c - nx] * Q[c - nx];
        Q[c] = t * M[c];
      }
      for (let k = m - 1; k >= 0; k--) {
        const c = C[k];
        const t = Q[c] - Ai[c] * M[c] * Z[c + 1] - Aj[c] * M[c] * Z[c + nx];
        Z[c] = t * M[c];
      }
    }
    applyPrecon();
    let sigma = 0;
    for (let k = 0; k < m; k++) {
      const c = C[k];
      S[c] = Z[c];
      sigma += Z[c] * R[c];
    }
    const maxIt = Math.max(1, params.pcgMaxIters | 0);
    let it = 0;
    let converged = false;
    let relRes = 1;
    while (it < maxIt) {
      it++;
      // Q = A S
      let sAs = 0;
      for (let k = 0; k < m; k++) {
        const c = C[k];
        const q = Ad[c] * S[c] + Ai[c] * S[c + 1] + Ai[c - 1] * S[c - 1] + Aj[c] * S[c + nx] + Aj[c - nx] * S[c - nx];
        Q[c] = q;
        sAs += S[c] * q;
      }
      if (!(sAs > 0)) break; // breakdown guard (should not happen: A is SPD with air present)
      const alpha = sigma / sAs;
      let r2 = 0;
      maxR = 0;
      for (let k = 0; k < m; k++) {
        const c = C[k];
        P[c] += alpha * S[c];
        const r = R[c] - alpha * Q[c];
        R[c] = r;
        r2 += r * r;
        const ar = r < 0 ? -r : r;
        if (ar > maxR) maxR = ar;
      }
      relRes = Math.sqrt(r2) / bNorm;
      if (relRes <= PCG_REL_TOL || maxR * resScale <= PCG_ABS_TOL) {
        converged = true;
        break;
      }
      applyPrecon();
      let sigmaNew = 0;
      for (let k = 0; k < m; k++) sigmaNew += Z[C[k]] * R[C[k]];
      const beta = sigmaNew / sigma;
      sigma = sigmaNew;
      for (let k = 0; k < m; k++) {
        const c = C[k];
        S[c] = Z[c] + beta * S[c];
      }
    }
    lastPcg.iters = it;
    lastPcg.relRes = relRes;
    lastPcg.maxResDt = maxR * resScale;
    if (!converged) {
      lastPcg.capHit = true;
      const rec = { frame, t: simTime, iters: it, relRes, maxResDt: maxR * resScale, fluidCells: m };
      if (pcgCapLog.length < 10000) pcgCapLog.push(rec);
      if (typeof console !== "undefined" && frame - pcgWarnFrame >= 60) {
        pcgWarnFrame = frame;
        console.warn("PCG safety cap hit: frame " + frame + ", " + it + " iterations, rel. residual " + relRes.toExponential(2) + ", max|div|*dt " + (maxR * resScale).toExponential(2));
      }
    }
    // apply pressure: u -= scale_f * Δp  with the SAME scale_f used in A (dig 1 / NM).
    const WR2 = wallRelease;
    const useW = weighted;
    for (let j = 0; j < ny; j++) {
      for (let i = 1; i < nx; i++) {
        const cr = i + j * nx;
        const cl = cr - 1;
        const tl = T[cl], tr = T[cr];
        if (tl !== FLUID && tr !== FLUID) continue;
        const f = i + j * su;
        const sf = useW ? uScale[f] : 1; // identical to Ad/Ai entries for this face
        if (tl === SOLID || tr === SOLID) {
          if (tr === FLUID && tl === SOLID && (WR2[cr] & 1)) U[f] -= sf * P[cr];
          else if (tl === FLUID && tr === SOLID && (WR2[cl] & 2)) U[f] -= sf * (-P[cl]);
          continue;
        }
        U[f] -= sf * (P[cr] - P[cl]);
      }
    }
    for (let j = 1; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const cb = i + j * nx;
        const ct = cb - nx;
        const tt = T[ct], tb = T[cb];
        if (tt !== FLUID && tb !== FLUID) continue;
        const sf = useW ? vScale[cb] : 1;
        if (tt === SOLID || tb === SOLID) {
          if (tb === FLUID && tt === SOLID && (WR2[cb] & 4)) V[cb] -= sf * P[cb];
          else if (tt === FLUID && tb === SOLID && (WR2[ct] & 8)) V[cb] -= sf * (-P[ct]);
          continue;
        }
        V[cb] -= sf * (P[cb] - P[ct]);
      }
    }
  }

  /**
   * After projection only faces touching a fluid cell hold divergence-free
   * velocities. Faces between two air cells still carry the raw P2G value plus
   * gravity, and the quadratic stencil of surface particles reaches them. Without
   * this step that unprojected gravity leaks into the surface every substep and a
   * resting pool slowly heats up. Fix: overwrite those faces with the average of
   * valid neighbours (2 layers covers the 1.5-cell stencil reach).
   */
  function extrapolateComponent(G, valid, next, stride, rows, isU) {
    const T = cellType;
    const cols = stride;
    valid.fill(0);
    for (let j = 0; j < rows; j++) {
      for (let i = 0; i < cols; i++) {
        const k = i + j * stride;
        let a, b; // the two cells sharing this face
        if (isU) {
          if (i === 0 || i === cols - 1) continue;
          a = i - 1 + j * nx;
          b = i + j * nx;
        } else {
          if (j === 0 || j === rows - 1) continue;
          a = i + (j - 1) * nx;
          b = i + j * nx;
        }
        if (T[a] === FLUID || T[b] === FLUID) valid[k] = 1;
      }
    }
    for (let layer = 0; layer < 2; layer++) {
      next.set(valid);
      for (let j = 1; j < rows - 1; j++) {
        for (let i = 1; i < cols - 1; i++) {
          const k = i + j * stride;
          if (valid[k]) continue;
          let s = 0;
          let n = 0;
          if (valid[k - 1]) { s += G[k - 1]; n++; }
          if (valid[k + 1]) { s += G[k + 1]; n++; }
          if (valid[k - stride]) { s += G[k - stride]; n++; }
          if (valid[k + stride]) { s += G[k + stride]; n++; }
          if (n > 0) {
            G[k] = s / n;
            next[k] = 1;
          }
        }
      }
      valid.set(next);
    }
  }
  function extrapolateVelocities() {
    extrapolateComponent(u, uValid, uValidNext, nx + 1, ny, true);
    extrapolateComponent(v, vValid, vValidNext, nx, ny + 1, false);
  }

  // ---- fixed solid obstacles ------------------------------------------------
  // Obstacles are rasterized into the solid-cell mask (a cell is solid iff its centre
  // is inside a shape) and joined with the box walls, so every grid step that already
  // handles wall cells (P2G masses, density sampler, drift, PCG, separating active set,
  // extrapolation) treats them the same way:
  //   - normal faces between an obstacle cell and a non-solid cell are walls: 0, or
  //     released by the separating active set exactly like a box-wall face;
  //   - faces with an obstacle cell on both sides are massless ghosts (their P2G mass
  //     is dropped) refilled from the fluid side: normal neighbours contribute 0 and
  //     tangential neighbours their value, which reproduces the box walls' free-slip
  //     mirror on flat sides (deep normal ghost = 0, tangential ghost = copy).
  // Particles collide with the exact shape (analytic SDF), like the wall clamp:
  // pushed out along the gradient to distance >= radius, the velocity component into
  // the solid removed, tangential kept.
  let obstacles = []; // { id, shape (normalized), comp (compiled) }
  let nextObstacleId = 1;
  let obstacleVersion = 0;
  let obstacleView = [];
  let obstacleViewVersion = -1;
  let obsActive = false;
  let obsMask = null; // Uint8Array nc: 1 = obstacle cell
  let baseType = null; // Uint8Array nc: SOLID for wall/obstacle cells, else AIR (markCells starts here)
  let nearMask = null; // Uint8Array nc: cells within radius + 1 cell of an obstacle's AABB
  let obsIfU = new Int32Array(0), obsIfUCell = new Int32Array(0), obsIfUBit = new Uint8Array(0);
  let obsIfV = new Int32Array(0), obsIfVCell = new Int32Array(0), obsIfVBit = new Uint8Array(0);
  const emptyGhosts = { faces: new Int32Array(0), src: new Int32Array(0), n: 0, thin: new Uint8Array(0) };
  let ghostU = emptyGhosts, ghostV = emptyGhosts;
  let obstacleDeleted = 0; // particles removed because an obstacle was placed on them (since reset)
  let lastObstacleDeleted = 0; // ... by the last addObstacle / updateObstacle
  let dropClipped = 0; // held-drop particles released inside an obstacle (removed)
  let pushStuck = 0; // particles wedged by obstacle + wall that needed the free-spot search
  const sdfOut = { nx: 0, ny: 0 };
  const sdfScratch = { nx: 0, ny: 0 };

  /** Min signed distance over all obstacles; shapes whose AABB is farther than pad are skipped. */
  function obstacleSdfAt(x, y, pad, out) {
    let best = Infinity, bnx = 0, bny = -1;
    for (let k = 0; k < obstacles.length; k++) {
      const c = obstacles[k].comp;
      if (x < c.minX - pad || x > c.maxX + pad || y < c.minY - pad || y > c.maxY + pad) continue;
      const d = ObstacleShapes.sdfCompiled(c, x, y, sdfScratch);
      if (d < best) { best = d; bnx = sdfScratch.nx; bny = sdfScratch.ny; }
    }
    if (out) { out.nx = bnx; out.ny = bny; }
    return best;
  }

  /** Ghost faces of one MAC component, ordered by distance (levels) from the fluid side. */
  function buildGhosts(cols, rows, isGhost, solidCount) {
    const nF = cols * rows;
    const level = new Uint8Array(nF); // 0 = not a ghost
    const list = [];
    for (let j = 0; j < rows; j++) for (let i = 0; i < cols; i++) {
      if (isGhost(i, j)) { level[i + j * cols] = 255; list.push(i + j * cols); }
    }
    if (!list.length) return emptyGhosts;
    const nbrs = (k) => {
      const i = k % cols, j = (k / cols) | 0, out = [];
      if (i > 0) out.push(k - 1);
      if (i < cols - 1) out.push(k + 1);
      if (j > 0) out.push(k - cols);
      if (j < rows - 1) out.push(k + cols);
      return out;
    };
    for (const k of list) if (nbrs(k).some((m) => level[m] === 0 && solidCount(m) < 2)) level[k] = 1;
    for (const k of list) if (level[k] === 255 && nbrs(k).some((m) => level[m] === 1)) level[k] = 2;
    const ordered = list.filter((k) => level[k] === 1).concat(list.filter((k) => level[k] === 2), list.filter((k) => level[k] === 255));
    const faces = Int32Array.from(ordered);
    const src = new Int32Array(4 * faces.length).fill(-1);
    for (let q = 0; q < faces.length; q++) {
      const k = faces[q], lv = level[k];
      let t = 0;
      for (const m of nbrs(k)) {
        if (lv === 1) {
          if (level[m] !== 0) continue;
          const sc = solidCount(m);
          if (sc === 0) src[4 * q + t++] = m; // fluid-side face: copy (tangential)
          else if (sc === 1) src[4 * q + t++] = -2; // a wall/obstacle normal face: counts as 0
        } else if (lv === 2 && level[m] === 1) src[4 * q + t++] = m;
      }
    }
    // [item 1] "thin" ghosts: a level-1 ghost with fluid-side (copy) faces on BOTH opposite
    // sides, i.e. inside a solid exactly one cell thick. Its single value is read by the
    // stencils of both sides, so it cannot mirror either side alone (see fillGhosts).
    const thin = new Uint8Array(faces.length);
    for (let q = 0; q < faces.length; q++) {
      if (level[faces[q]] !== 1) continue;
      let nc = 0;
      for (let t = 0; t < 4; t++) if (src[4 * q + t] >= 0) nc++;
      thin[q] = nc >= 2 ? 1 : 0;
    }
    return { faces, src, n: faces.length, thin };
  }

  /**
   * Refill ghost faces. Thin ghosts [item 1] weight their fluid-side faces by P2G face mass
   * (W = uW / vW, fixed for the whole substep) instead of equally: a dry side has mass 0 and
   * never contributes, so the ghost holds the same function of the wet side at P2G time (the
   * FLIP reference) and after extrapolation. With equal weights the dry-side face is 0 at P2G
   * time but a copy of the wet side after extrapolation, so the FLIP delta at the ghost gained
   * +v/2 every substep (exponential energy blow-up next to one-cell-thick obstacles).
   * Normal-face sources keep their 0 contribution (factor nCopy/count). Both sides dry: the
   * old equal-weight rule. Non-thin ghosts (every obstacle >= 2 cells thick) are unchanged.
   */
  function fillGhosts(G, g, W) {
    const F = g.faces, S = g.src, TH = g.thin;
    for (let q = 0; q < g.n; q++) {
      if (TH[q]) {
        let ws = 0, wv = 0, c = 0, nc = 0;
        for (let t = 4 * q, e = t + 4; t < e; t++) {
          const k = S[t];
          if (k === -1) continue;
          c++;
          if (k >= 0) { nc++; ws += W[k]; wv += W[k] * G[k]; }
        }
        if (ws > 0) { G[F[q]] = (wv / ws) * (nc / c); continue; }
      }
      let s = 0, c = 0;
      for (let t = 4 * q, e = t + 4; t < e; t++) {
        const k = S[t];
        if (k === -1) continue;
        c++;
        if (k >= 0) s += G[k];
      }
      G[F[q]] = c ? s / c : 0;
    }
  }

  /** Re-rasterize obstacles into the solid mask and rebuild the face lists. */
  function rebuildObstacleGrid() {
    if (!baseType) return; // before init()
    obstacleVersion++;
    obsMask.fill(0);
    nearMask.fill(0);
    for (const o of obstacles) {
      const c = o.comp;
      const i0 = Math.max(1, Math.floor(c.minX * invH - 0.5)), i1 = Math.min(nx - 2, Math.ceil(c.maxX * invH - 0.5));
      const j0 = Math.max(1, Math.floor(c.minY * invH - 0.5)), j1 = Math.min(ny - 2, Math.ceil(c.maxY * invH - 0.5));
      for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
        if (ObstacleShapes.sdfCompiled(c, (i + 0.5) * h, (j + 0.5) * h, null) < 0) obsMask[i + j * nx] = 1;
      }
      const pad = pr + h;
      const a0 = Math.max(0, Math.floor((c.minX - pad) * invH)), a1 = Math.min(nx - 1, Math.floor((c.maxX + pad) * invH));
      const b0 = Math.max(0, Math.floor((c.minY - pad) * invH)), b1 = Math.min(ny - 1, Math.floor((c.maxY + pad) * invH));
      for (let j = b0; j <= b1; j++) for (let i = a0; i <= a1; i++) nearMask[i + j * nx] = 1;
    }
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const c = i + j * nx;
      baseType[c] = i === 0 || j === 0 || i === nx - 1 || j === ny - 1 || obsMask[c] ? SOLID : AIR;
    }
    // wall-like normal faces between an obstacle cell and a non-solid cell
    const isSolid = (c) => baseType[c] === SOLID;
    const fu = [], fuc = [], fub = [], fv = [], fvc = [], fvb = [];
    for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx; i++) {
      const cl = i - 1 + j * nx, cr = cl + 1, k = i + j * (nx + 1);
      if (obsMask[cl] && !isSolid(cr)) { fu.push(k); fuc.push(cr); fub.push(1); }
      else if (obsMask[cr] && !isSolid(cl)) { fu.push(k); fuc.push(cl); fub.push(2); }
    }
    for (let j = 1; j < ny; j++) for (let i = 1; i < nx - 1; i++) {
      const ct = i + (j - 1) * nx, cb = ct + nx, k = cb;
      if (obsMask[ct] && !isSolid(cb)) { fv.push(k); fvc.push(cb); fvb.push(4); }
      else if (obsMask[cb] && !isSolid(ct)) { fv.push(k); fvc.push(ct); fvb.push(8); }
    }
    obsIfU = Int32Array.from(fu); obsIfUCell = Int32Array.from(fuc); obsIfUBit = Uint8Array.from(fub);
    obsIfV = Int32Array.from(fv); obsIfVCell = Int32Array.from(fvc); obsIfVBit = Uint8Array.from(fvb);
    // ghost faces strictly inside obstacles (box-wall lines are left to enforceWalls)
    const su = nx + 1;
    const cellSolid = (i, j) => (i < 0 || j < 0 || i >= nx || j >= ny ? 1 : isSolid(i + j * nx) ? 1 : 0);
    ghostU = buildGhosts(su, ny,
      (i, j) => i >= 2 && i <= nx - 2 && j >= 1 && j <= ny - 2 && obsMask[i - 1 + j * nx] && obsMask[i + j * nx],
      (k) => { const i = k % su, j = (k / su) | 0; return cellSolid(i - 1, j) + cellSolid(i, j); });
    ghostV = buildGhosts(nx, ny + 1,
      (i, j) => j >= 2 && j <= ny - 2 && i >= 1 && i <= nx - 2 && obsMask[i + (j - 1) * nx] && obsMask[i + j * nx],
      (k) => { const i = k % nx, j = (k / nx) | 0; return cellSolid(i, j - 1) + cellSolid(i, j); });
    if (wallRelease) wallRelease.fill(0);
    obsActive = obstacles.length > 0;
  }

  /** Obstacle part of enforceBoundary (only called when obstacles exist). */
  function enforceObstacles() {
    const mode = params.boundary;
    const sep = mode === "separating", legacy = mode === "legacy";
    const WR = wallRelease;
    for (let q = 0; q < obsIfU.length; q++) {
      const k = obsIfU[q], bit = obsIfUBit[q];
      if (sep) { if (!(WR[obsIfUCell[q]] & bit)) u[k] = 0; }
      else if (legacy) { if (bit === 1 ? u[k] < 0 : u[k] > 0) u[k] = 0; } // one-sided: only flow into the solid removed
      else u[k] = 0;
    }
    for (let q = 0; q < obsIfV.length; q++) {
      const k = obsIfV[q], bit = obsIfVBit[q];
      if (sep) { if (!(WR[obsIfVCell[q]] & bit)) v[k] = 0; }
      else if (legacy) { if (bit === 4 ? v[k] < 0 : v[k] > 0) v[k] = 0; }
      else v[k] = 0;
    }
    fillGhosts(u, ghostU, uW);
    fillGhosts(v, ghostV, vW);
  }

  /** P2G: faces inside obstacles carry no mass (like the folded box-wall ghosts). */
  function clearObstacleGhostMass() {
    let F = ghostU.faces;
    for (let q = 0; q < ghostU.n; q++) { u[F[q]] = 0; uW[F[q]] = 0; }
    F = ghostV.faces;
    for (let q = 0; q < ghostV.n; q++) { v[F[q]] = 0; vW[F[q]] = 0; }
  }

  /** Free spot (sdf >= radius, inside the box) near (x, y): rings of 16 directions. */
  function findFreeSpot(x, y, out) {
    const minX = h + pr, maxX = (nx - 1) * h - pr, minY = h + pr, maxY = (ny - 1) * h - pr;
    for (let ring = 1; ring <= 10; ring++) {
      const r = ring * 0.4 * h;
      for (let a = 0; a < 16; a++) {
        const qx = x + r * Math.cos(a * 0.39269908), qy = y + r * Math.sin(a * 0.39269908);
        if (qx < minX || qx > maxX || qy < minY || qy > maxY) continue;
        if (obstacleSdfAt(qx, qy, pr, null) >= pr) { out.x = qx; out.y = qy; return true; }
      }
    }
    return false;
  }

  /**
   * After advection + wall clamp: any free particle closer than its radius to an
   * obstacle is projected out along the SDF gradient (up to 4 passes, re-clamped to
   * the box each pass) and its velocity into the solid is removed (tangential kept),
   * exactly like collideWalls does for the box.
   */
  function pushOutObstacles() {
    const X = px, Y = py, VX = pvx, VY = pvy, HD = held, NM = nearMask, ih = invH;
    const minX = h + pr, maxX = (nx - 1) * h - pr, minY = h + pr, maxY = (ny - 1) * h - pr;
    const o = sdfOut;
    const spot = { x: 0, y: 0 };
    for (let i = 0; i < count; i++) {
      if (HD[i]) continue;
      let x = X[i], y = Y[i];
      let ci = (x * ih) | 0, cj = (y * ih) | 0;
      if (ci < 0) ci = 0; else if (ci >= nx) ci = nx - 1;
      if (cj < 0) cj = 0; else if (cj >= ny) cj = ny - 1;
      if (!NM[ci + cj * nx]) continue;
      let d = obstacleSdfAt(x, y, pr, o);
      if (d >= pr) continue;
      let vx = VX[i], vy = VY[i];
      for (let it = 0; it < 4 && d < pr; it++) {
        const push = pr - d;
        x += o.nx * push;
        y += o.ny * push;
        const vn = vx * o.nx + vy * o.ny;
        if (vn < 0) { vx -= vn * o.nx; vy -= vn * o.ny; }
        if (x < minX) { x = minX; if (vx < 0) vx = 0; } else if (x > maxX) { x = maxX; if (vx > 0) vx = 0; }
        if (y < minY) { y = minY; if (vy < 0) vy = 0; } else if (y > maxY) { y = maxY; if (vy > 0) vy = 0; }
        d = obstacleSdfAt(x, y, pr, o);
      }
      if (d < 0) {
        // wedged between an obstacle and a wall (or two obstacles): jump to the nearest free spot
        pushStuck++;
        if (findFreeSpot(x, y, spot)) { x = spot.x; y = spot.y; vx = 0; vy = 0; }
      }
      X[i] = x; Y[i] = y; VX[i] = vx; VY[i] = vy;
    }
  }

  /**
   * [item 2] Obstacle contact at the source (runs before advection). For every obstacle a
   * particle could reach this substep, with n = the obstacle's SDF gradient (outward unit
   * normal) and d = its signed distance:
   *   - d <= r + m (touching; m = APPROACH_MARGIN = 4e-5 h, a few float32 ulps of a position,
   *     so rounding of the stored position cannot carry a touching particle inside): strip the
   *     whole into-solid normal component, v -= (v·n) n if v·n < 0;
   *   - d > r + m: limit the normal approach to what keeps the centre at least r + m from the
   *     surface after advection, v·n >= -(d - r - m)/dt (only the excess into-solid part is removed).
   * Only an into-solid normal component is ever reduced (Δv is parallel to n and v·n stays <= 0):
   * outward normal velocity and tangential velocity are never touched (frictionless free-slip).
   * For a convex obstacle the SDF outside it is convex, so d(x + v dt) >= d + (v·n) dt >= r + m:
   * advection cannot carry a centre into the contact band, and pushOutObstacles() (which lifts
   * particles along n and so raised potential energy) is left as a rare backstop for the
   * separation pass, concave freehand shapes and wedges.
   */
  const approachOut = { nx: 0, ny: 0 };
  const APPROACH_MARGIN_CELLS = 4e-5;
  function limitObstacleApproach(dt) {
    const X = px, Y = py, VX = pvx, VY = pvy, HD = held, o = approachOut, nObs = obstacles.length;
    const rc = pr + APPROACH_MARGIN_CELLS * h; // contact distance incl. the rounding margin
    for (let i = 0; i < count; i++) {
      if (HD[i]) continue;
      const x = X[i], y = Y[i];
      let vx = VX[i], vy = VY[i];
      const reach = rc + Math.sqrt(vx * vx + vy * vy) * dt; // farther obstacles cannot be reached
      let changed = false;
      for (let k = 0; k < nObs; k++) {
        const c = obstacles[k].comp;
        if (x < c.minX - reach || x > c.maxX + reach || y < c.minY - reach || y > c.maxY + reach) continue;
        const d = ObstacleShapes.sdfCompiled(c, x, y, o);
        if (d >= reach) continue;
        const vn = vx * o.nx + vy * o.ny;
        if (vn >= 0) continue;
        const lim = d > rc ? -(d - rc) / dt : 0;
        if (vn < lim) { const dv = lim - vn; vx += dv * o.nx; vy += dv * o.ny; changed = true; }
      }
      if (changed) { VX[i] = vx; VY[i] = vy; }
    }
  }

  /** Stable in-place removal of particles where keep(i) is false. Returns the count removed. */
  function compactParticles(keep) {
    let w = 0;
    for (let r = 0; r < count; r++) {
      if (!keep(r)) continue;
      if (w !== r) {
        px[w] = px[r]; py[w] = py[r]; pvx[w] = pvx[r]; pvy[w] = pvy[r];
        c00[w] = c00[r]; c01[w] = c01[r]; c10[w] = c10[r]; c11[w] = c11[r];
        held[w] = held[r]; holdOx[w] = holdOx[r]; holdOy[w] = holdOy[r];
      }
      w++;
    }
    const removed = count - w;
    count = w;
    return removed;
  }
  function removeParticlesInside(comp) {
    if (!count) return 0;
    return compactParticles((i) => {
      const x = px[i], y = py[i];
      if (x < comp.minX || x > comp.maxX || y < comp.minY || y > comp.maxY) return true;
      return !(ObstacleShapes.sdfCompiled(comp, x, y, null) < 0);
    });
  }
  function findObstacle(id) {
    for (let k = 0; k < obstacles.length; k++) if (obstacles[k].id === id) return k;
    return -1;
  }
  /**
   * Add a fixed obstacle (see ObstacleShapes for the shape format). Particles whose
   * centre is inside it are deleted (count in lastRemovedParticles / stats). Returns id.
   * opts.removeParticles = false keeps them (they are then pushed out by the solver).
   */
  function addObstacle(shape, opts) {
    const s = ObstacleShapes.normalize(shape);
    const rec = { id: nextObstacleId++, shape: s, comp: ObstacleShapes.compile(s) };
    obstacles.push(rec);
    rebuildObstacleGrid();
    lastObstacleDeleted = opts && opts.removeParticles === false ? 0 : removeParticlesInside(rec.comp);
    obstacleDeleted += lastObstacleDeleted;
    return rec.id;
  }
  /** Replace an obstacle's shape (move / rotate). Same particle rule as addObstacle. */
  function updateObstacle(id, shape, opts) {
    const k = findObstacle(id);
    if (k < 0) return false;
    const s = ObstacleShapes.normalize(shape);
    obstacles[k] = { id, shape: s, comp: ObstacleShapes.compile(s) };
    rebuildObstacleGrid();
    lastObstacleDeleted = opts && opts.removeParticles === false ? 0 : removeParticlesInside(obstacles[k].comp);
    obstacleDeleted += lastObstacleDeleted;
    return true;
  }
  function removeObstacle(id) {
    const k = findObstacle(id);
    if (k < 0) return false;
    obstacles.splice(k, 1);
    rebuildObstacleGrid();
    return true;
  }
  function clearObstacles() {
    obstacles = [];
    rebuildObstacleGrid();
  }
  /** Topmost obstacle id whose shape contains (x, y) (or is within tol m), else -1. */
  function hitTest(x, y, tol) {
    const t = tol > 0 ? tol : 0;
    for (let k = obstacles.length - 1; k >= 0; k--) {
      const c = obstacles[k].comp;
      if (x < c.minX - t || x > c.maxX + t || y < c.minY - t || y > c.maxY + t) continue;
      if (ObstacleShapes.sdfCompiled(c, x, y, null) <= t) return obstacles[k].id;
    }
    return -1;
  }

  // ---- time stepping ------------------------------------------------------
  function substep(dt) {
    pinHeld();
    if (obsActive) limitObstacleApproach(dt);
    advectParticles(dt);
    separateParticles(SEPARATION_ITERS);
    collideWalls();
    if (obsActive) pushOutObstacles();
    particlesToGrid();
    markCells();
    applyGravity(dt);
    applyViscosity();
    enforceBoundary();
    solvePressure(dt);
    extrapolateVelocities();
    enforceBoundary();
    if (substepHook) substepHook(divergenceStats(dt), dt);
    gridToParticles();
  }

  /**
   * Post-projection divergence over fluid cells (diagnostics only):
   * |div u|*dt (raw, includes the intended drift divergence) and the solver
   * residual |div u - drift target|*dt.
   */
  function divergenceStats(dt) {
    const su = nx + 1;
    let maxDiv = 0, sumDiv = 0, maxRes = 0, sumRes = 0, n = 0;
    const rest = restDensity;
    const driftRate = params.driftRate != null ? params.driftRate : DRIFT_RATE;
    const driftScale = rest > 0 ? (driftRate * h) / rest : 0;
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx;
        if (cellType[c] !== FLUID) continue;
        const ku = i + j * su;
        const d = u[ku + 1] - u[ku] + v[c + nx] - v[c]; // = h * div u
        let comp = cellDensity[c] - rest;
        if (comp < 0 && (cellType[c - 1] === AIR || cellType[c + 1] === AIR || cellType[c - nx] === AIR || cellType[c + nx] === AIR)) comp = 0;
        const a = (Math.abs(d) * dt) / h;
        const r = (Math.abs(d - driftScale * comp) * dt) / h;
        if (a > maxDiv) maxDiv = a;
        if (r > maxRes) maxRes = r;
        sumDiv += a;
        sumRes += r;
        n++;
      }
    }
    return {
      t: simTime + dt, fluidCells: n, maxDiv, meanDiv: n ? sumDiv / n : 0, maxRes, meanRes: n ? sumRes / n : 0,
      solver: params.pressureSolver, omega: OVER_RELAX,
      pcgIters: lastPcg.iters, pcgRelRes: lastPcg.relRes, pcgMaxResDt: lastPcg.maxResDt, pcgCapHit: lastPcg.capHit,
      clampHits: lastClamp.hits, clampFaces: lastClamp.faces, clampFraction: lastClamp.fraction,
      sepPasses: lastSep.passes, sepCapHit: lastSep.capHit, sepNegPLeft: lastSep.negPLeft,
      sepWorstNegP: lastSep.worstNegP, sepWorstUInto: lastSep.worstUInto, sepCompFail: lastSep.compFail,
      separation: params.separation, sepRawDx: sepCom.rawDx, sepRawDy: sepCom.rawDy, sepEffDx: sepCom.effDx, sepEffDy: sepCom.effDy,
      umax: lastUmax, cfl: (lastUmax * dt) / h, cflFull: ((lastUmax + Math.sqrt(Math.max(0, params.gravity) * h)) * dt) / h,
    };
  }

  function maxFreeSpeed() {
    let m2 = 0;
    for (let i = 0; i < count; i++) {
      if (held[i]) continue; // pinned particles do not advect
      const s2 = pvx[i] * pvx[i] + pvy[i] * pvy[i];
      if (s2 > m2) m2 = s2;
    }
    return Math.sqrt(m2);
  }

  /**
   * Advance one frame (default FIXED_DT = 1/60 s). Each substep uses
   *   dt = CFL * h / (u_max + sqrt(g h)),  CFL = 1
   * (spread evenly over the rest of the frame). At most MAX_SUBSTEPS per frame:
   * if the frame is not finished by then it is logged in capLog and the sim
   * simply falls behind real time; CFL 1 is never exceeded.
   */
  function step(frameDt) {
    const target = frameDt > 0 ? frameDt : FIXED_DT;
    frame++;
    const cRef = Math.sqrt(Math.max(0, params.gravity) * h);
    let remaining = target;
    let n = 0;
    let frameMaxCfl = 0;
    while (remaining > 1e-12 && n < MAX_SUBSTEPS) {
      const umax = maxFreeSpeed();
      lastUmax = umax;
      const dtCfl = (CFL * h) / (umax + cRef + 1e-12);
      let dt;
      if (dtCfl >= remaining) dt = remaining;
      else if (MAX_SUBSTEPS - n > 1) dt = remaining / Math.ceil(remaining / dtCfl); // even split
      else dt = dtCfl;
      if (dt > dtCfl) dt = dtCfl;
      const cfl = (umax * dt) / h; // advective CFL of this substep (particles move u*dt)
      if (cfl > frameMaxCfl) frameMaxCfl = cfl;
      if (cfl > maxCflObserved) maxCflObserved = cfl;
      substep(dt);
      simTime += dt;
      remaining -= dt;
      n++;
    }
    lastSubsteps = n;
    if (remaining > 1e-12) {
      const rec = { frame, t: simTime, umax: lastUmax, behind: remaining, covered: target - remaining, substeps: n, maxCfl: frameMaxCfl };
      if (capLog.length < 10000) capLog.push(rec);
      if (typeof console !== "undefined" && frame - capWarnFrame >= 60) {
        capWarnFrame = frame; // at most one console line per ~second
        console.warn("FLIP substep cap hit: frame " + frame + ", u_max " + lastUmax.toFixed(2) + " m/s, sim behind by " + (remaining * 1000).toFixed(2) + " ms");
      }
    }
  }

  /** Diagnostics (used by test/headless.js). */
  function stats() {
    let nan = 0;
    let outside = 0;
    let vmax = 0;
    let ke = 0;
    const eps = 1e-6;
    const minX = h + pr - eps, maxX = (nx - 1) * h - pr + eps;
    const minY = h + pr - eps, maxY = (ny - 1) * h - pr + eps;
    for (let i = 0; i < count; i++) {
      const a = px[i], b = py[i], c = pvx[i], d = pvy[i];
      if (!Number.isFinite(a) || !Number.isFinite(b) || !Number.isFinite(c) || !Number.isFinite(d)) nan++;
      if (!held[i] && (a < minX || a > maxX || b < minY || b > maxY)) outside++;
      const s2 = c * c + d * d;
      if (s2 > vmax * vmax) vmax = Math.sqrt(s2);
      ke += 0.5 * s2;
    }
    // particles whose centre is inside an obstacle (should stay 0) and the closest approach
    let insideObstacle = 0;
    let minObstacleSdf = Infinity;
    if (obsActive) {
      for (let i = 0; i < count; i++) {
        if (held[i]) continue;
        const d = obstacleSdfAt(px[i], py[i], 2 * pr, null);
        if (d < 0) insideObstacle++;
        if (d < minObstacleSdf) minObstacleSdf = d;
      }
    }
    let obstacleCells = 0;
    if (obsActive) for (let c = 0; c < obsMask.length; c++) obstacleCells += obsMask[c];
    let fluidCells = 0;
    let interior = 0;
    let interiorDensity = 0;
    let interiorDensityLegacy = 0;
    for (let j = 1; j < ny - 1; j++) {
      for (let i = 1; i < nx - 1; i++) {
        const c = i + j * nx;
        if (cellType[c] !== FLUID) continue;
        fluidCells++;
        if (cellType[c - 1] === AIR || cellType[c + 1] === AIR || cellType[c - nx] === AIR || cellType[c + nx] === AIR) continue;
        interior++;
        interiorDensity += cellDensity[c];
        interiorDensityLegacy += cellDensityLegacy[c];
      }
    }
    return {
      count,
      nan,
      outside,
      vmax,
      ke: count ? ke / count : 0,
      fluidCells,
      interiorDensity: interior ? interiorDensity / interior : 0,
      interiorDensityLegacy: interior ? interiorDensityLegacy / interior : 0,
      restDensity,
      substeps: lastSubsteps,
      simTime,
      frame,
      umax: lastUmax,
      obstacles: obstacles.length,
      obstacleCells,
      insideObstacle,
      minObstacleSdf, // closest particle-centre distance to an obstacle within 2 radii (Infinity = none)
      obstacleDeleted,
      lastObstacleDeleted,
      dropClipped,
      pushStuck,
    };
  }

  // ---- public API: renderers read the typed arrays, never write them --------
  return {
    params,
    FIXED_DT,
    DROP_RADIUS_CELLS,
    MAX_SUBSTEPS,
    CFL,
    OVER_RELAX,
    capLog,
    pcgCapLog,
    lastPcg,
    lastClamp,
    clampAccum,
    resetClampAccum() { clampAccum.hits = 0; clampAccum.faces = 0; },
    lastSep,
    sepPassLog,
    get sepStats() {
      const a = Float64Array.from(sepPassLog).sort();
      const p95 = a.length ? a[Math.min(a.length - 1, Math.ceil(0.95 * a.length) - 1)] : 0;
      return {
        mean: sepSubsteps ? sepPassSum / sepSubsteps : 0,
        p95, max: sepPassMax, capHits: sepCapHits, substeps: sepSubsteps,
      };
    },
    PCG_REL_TOL,
    PCG_ABS_TOL,
    init,
    reset,
    resetColumn,
    step,
    spawnDrop,
    releaseDrop,
    setPointer,
    pinHeld,
    stats,
    get count() { return count; },
    get px() { return px; }, // Float32Array positions (m), y down
    get py() { return py; },
    get vx() { return pvx; }, // Float32Array velocities (m/s)
    get vy() { return pvy; },
    get held() { return held; },
    get capacity() { return maxN; },
    get nx() { return nx; }, // grid size incl. 1 solid wall cell per side
    get ny() { return ny; },
    get cellSize() { return h; },
    get radius() { return pr; },
    get holding() { return holding; },
    get simTime() { return simTime; },
    get frame() { return frame; },
    get nanRecoveries() { return nanRecoveries; },
    get maxCflObserved() { return maxCflObserved; },
    get restDensity() { return restDensity; },
    get cellType() { return cellType; }, // Uint8Array nx*ny: 0 solid, 1 fluid, 2 air
    get cellDensity() { return cellDensity; },
    get cellDensityLegacy() { return cellDensityLegacy; },
    setSubstepHook(fn) { substepHook = fn || null; },
    get pointerX() { return pointerX; },
    get pointerY() { return pointerY; },
    // ---- fixed obstacles ----
    addObstacle,
    updateObstacle,
    removeObstacle,
    clearObstacles,
    hitTest,
    /** Min signed distance (m) from (x, y) to any obstacle; Infinity when there are none. */
    obstacleSdf(x, y) { return obstacleSdfAt(x, y, Infinity, null); },
    /** Read-only copies: [{ id, type, cx, cy, angle, ...params }], rebuilt only on change. */
    get obstacles() {
      if (obstacleViewVersion !== obstacleVersion) {
        obstacleView = obstacles.map((o) => Object.freeze(Object.assign({ id: o.id }, ObstacleShapes.clone(o.shape))));
        obstacleViewVersion = obstacleVersion;
      }
      return obstacleView;
    },
    get obstacleVersion() { return obstacleVersion; },
    get obstacleMask() { return obsMask; }, // Uint8Array nx*ny, 1 = rasterized obstacle cell
    /** Editor rule check on this grid (no side effects): see ObstacleShapes.thicknessCheck. */
    checkObstacleThickness(shape) { return ObstacleShapes.thicknessCheck(shape, nx, ny, h); },
    get lastRemovedParticles() { return lastObstacleDeleted; },
    get obstacleDeleted() { return obstacleDeleted; },
  };
})();

/**
 * Map a canvas (CSS px) to a tank: canvas width = TANK_WIDTH_M, cells ~10 px,
 * tank floor on the canvas bottom. Pure function; shared by the app and tests.
 */
function computeLayout(wPx, hPx) {
  const TANK_WIDTH_M = 3.0;
  const cellsX = Math.max(64, Math.min(128, Math.round(wPx / 10)));
  const cellSize = TANK_WIDTH_M / cellsX;
  const scale = wPx / TANK_WIDTH_M; // px per meter
  const cellsY = Math.max(8, Math.floor(hPx / (cellSize * scale)));
  return {
    cellsX,
    cellsY,
    cellSize,
    scale,
    // world (m) -> screen (px): sx = (x - originX) * scale + offX
    originX: cellSize, // first interior cell starts one wall cell in
    originY: cellSize,
    offX: 0,
    offY: hPx - cellsY * cellSize * scale,
  };
}

// ===========================================================================
// 2. Canvas2DRenderer — reads solver arrays, draws speed-coloured discs
//    Interface (keep for a future WebGL2 renderer): resize(w, h, dpr), draw(solver, view[, overlay])
// ===========================================================================
const Canvas2DRenderer = (function () {
  const BUCKETS = 16;
  const SPEED_WHITE = 4.0; // m/s that maps to white foam
  const DISC_SCALE = 1.35; // draw slightly fat discs so the body reads as continuous
  const OBSTACLE_FILL = "#2a323d";
  const OBSTACLE_EDGE = "rgba(200, 215, 230, 0.35)";

  function create(canvas) {
    const ctx = canvas.getContext("2d", { alpha: false });
    const colors = [];
    for (let b = 0; b < BUCKETS; b++) {
      const t = b / (BUCKETS - 1);
      const e = t * t * (3 - 2 * t);
      const r = Math.round(28 + (235 - 28) * e);
      const g = Math.round(110 + (246 - 110) * Math.sqrt(t));
      const bl = Math.round(205 + (255 - 205) * t);
      colors.push("rgb(" + r + "," + g + "," + bl + ")");
    }
    let bucketOf = new Uint8Array(0);
    let cssW = 1;
    let cssH = 1;

    function resize(w, hh, dpr) {
      cssW = w;
      cssH = hh;
      canvas.width = Math.floor(w * dpr);
      canvas.height = Math.floor(hh * dpr);
      canvas.style.width = w + "px";
      canvas.style.height = hh + "px";
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    /**
     * overlay (optional, App only): { editMode, pending, editingId, selectedId,
     * stroke, handle, hint, status } for Object Edit Mode. Obstacles are always drawn.
     */
    function draw(solver, view, overlay) {
      ctx.fillStyle = "#07090d";
      ctx.fillRect(0, 0, cssW, cssH);
      const n = solver.count;
      const X = solver.px, Y = solver.py, VX = solver.vx, VY = solver.vy;
      if (bucketOf.length < solver.capacity) bucketOf = new Uint8Array(solver.capacity);
      const s = view.scale;
      const ox = view.offX - view.originX * s;
      const oy = view.offY - view.originY * s;
      const rad = solver.radius * s * DISC_SCALE;
      const invSp = (BUCKETS - 1) / SPEED_WHITE;
      for (let i = 0; i < n; i++) {
        const b = (Math.sqrt(VX[i] * VX[i] + VY[i] * VY[i]) * invSp) | 0;
        bucketOf[i] = b > BUCKETS - 1 ? BUCKETS - 1 : b;
      }
      const tau = Math.PI * 2;
      for (let b = 0; b < BUCKETS; b++) {
        ctx.fillStyle = colors[b];
        ctx.beginPath();
        let any = false;
        for (let i = 0; i < n; i++) {
          if (bucketOf[i] !== b) continue;
          const sx = X[i] * s + ox;
          const sy = Y[i] * s + oy;
          ctx.moveTo(sx + rad, sy);
          ctx.arc(sx, sy, rad, 0, tau);
          any = true;
        }
        if (any) ctx.fill();
      }
      if (solver.holding) {
        ctx.strokeStyle = "rgba(180, 210, 230, 0.35)";
        ctx.lineWidth = 1;
        ctx.beginPath();
        const R = (solver.DROP_RADIUS_CELLS + 1) * solver.cellSize * s;
        ctx.arc(solver.pointerX * s + ox, solver.pointerY * s + oy, R, 0, tau);
        ctx.stroke();
      }
      drawObstacles(solver, s, ox, oy, overlay);
      if (overlay && overlay.editMode) drawEditOverlay(solver, s, ox, oy, overlay);
    }

    /** Trace a shape (solver units) as a closed path in screen space. */
    function traceShape(shape, s, ox, oy) {
      if (shape.type === "circle") {
        const cx = shape.cx * s + ox, cy = shape.cy * s + oy, r = shape.r * s;
        ctx.moveTo(cx + r, cy);
        ctx.arc(cx, cy, r, 0, Math.PI * 2);
        return;
      }
      const V = ObstacleShapes.worldVertices(shape);
      ctx.moveTo(V[0][0] * s + ox, V[0][1] * s + oy);
      for (let k = 1; k < V.length; k++) ctx.lineTo(V[k][0] * s + ox, V[k][1] * s + oy);
      ctx.closePath();
    }

    function drawObstacles(solver, s, ox, oy, overlay) {
      const list = solver.obstacles;
      if (!list || !list.length) return;
      const editingId = overlay ? overlay.editingId : -1;
      const selectedId = overlay ? overlay.selectedId : -1;
      for (const o of list) {
        const ghosted = o.id === editingId; // being moved: original shown faint
        ctx.beginPath();
        traceShape(o, s, ox, oy);
        ctx.globalAlpha = ghosted ? 0.3 : 1;
        ctx.fillStyle = OBSTACLE_FILL;
        ctx.fill();
        ctx.lineWidth = 1;
        ctx.strokeStyle = OBSTACLE_EDGE;
        ctx.stroke();
        ctx.globalAlpha = 1;
        if (o.id === selectedId && !ghosted) {
          ctx.lineWidth = 2;
          ctx.strokeStyle = "rgba(90, 167, 212, 0.95)";
          ctx.stroke();
        }
      }
    }

    function drawEditOverlay(solver, s, ox, oy, overlay) {
      // rasterized solid cells the grid actually sees (staircase), very faint
      const M = solver.obstacleMask;
      if (M) {
        const nx = solver.nx, ny = solver.ny, hh = solver.cellSize * s;
        ctx.fillStyle = "rgba(255, 255, 255, 0.06)";
        ctx.beginPath();
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
          if (M[i + j * nx]) ctx.rect(i * solver.cellSize * s + ox, j * solver.cellSize * s + oy, hh, hh);
        }
        ctx.fill();
      }
      const p = overlay.pending;
      if (p) {
        ctx.beginPath();
        traceShape(p, s, ox, oy);
        ctx.fillStyle = "rgba(90, 167, 212, 0.28)";
        ctx.fill();
        ctx.setLineDash([6, 4]);
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = "rgba(150, 205, 240, 0.95)";
        ctx.stroke();
        ctx.setLineDash([]);
        if (overlay.handle) {
          const cx = p.cx * s + ox, cy = p.cy * s + oy;
          const hx = overlay.handle.x * s + ox, hy = overlay.handle.y * s + oy;
          ctx.beginPath();
          ctx.moveTo(cx, cy);
          ctx.lineTo(hx, hy);
          ctx.strokeStyle = "rgba(150, 205, 240, 0.6)";
          ctx.lineWidth = 1;
          ctx.stroke();
          ctx.beginPath();
          ctx.arc(hx, hy, 7, 0, Math.PI * 2);
          ctx.fillStyle = "rgba(90, 167, 212, 0.9)";
          ctx.fill();
          ctx.strokeStyle = "#e8eef4";
          ctx.stroke();
        }
      }
      if (overlay.thin && overlay.thin.length) {
        // cells that fail the 2-cell thickness rule (rejected Confirm)
        const cs = solver.cellSize * s;
        ctx.beginPath();
        for (const c of overlay.thin) ctx.rect(c[0] * cs + ox, c[1] * cs + oy, cs, cs);
        ctx.fillStyle = "rgba(235, 80, 70, 0.55)";
        ctx.fill();
        ctx.lineWidth = 1;
        ctx.strokeStyle = "rgba(255, 140, 120, 0.95)";
        ctx.stroke();
      }
      const st = overlay.stroke;
      if (st && st.length > 1) {
        ctx.beginPath();
        ctx.moveTo(st[0][0] * s + ox, st[0][1] * s + oy);
        for (let k = 1; k < st.length; k++) ctx.lineTo(st[k][0] * s + ox, st[k][1] * s + oy);
        ctx.strokeStyle = "rgba(150, 205, 240, 0.95)";
        ctx.lineWidth = 2;
        ctx.stroke();
      }
      const lines = [overlay.hint, overlay.status].filter(Boolean);
      if (lines.length) drawHint(lines);
    }

    /** Word-wrapped hint box at the bottom centre of the canvas. */
    function drawHint(lines) {
      ctx.font = "12px system-ui, -apple-system, Segoe UI, sans-serif";
      const maxW = Math.min(cssW - 24, 760);
      const rows = [];
      for (const text of lines) {
        let cur = "";
        for (const word of text.split(" ")) {
          const t = cur ? cur + " " + word : word;
          if (cur && ctx.measureText(t).width > maxW - 20) { rows.push(cur); cur = word; } else cur = t;
        }
        if (cur) rows.push(cur);
      }
      let w = 0;
      for (const r of rows) w = Math.max(w, ctx.measureText(r).width);
      const lh = 16, bw = w + 20, bh = rows.length * lh + 12;
      const bx = (cssW - bw) / 2, by = cssH - bh - 14;
      ctx.fillStyle = "rgba(10, 13, 18, 0.82)";
      ctx.strokeStyle = "rgba(90, 167, 212, 0.45)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(bx, by, bw, bh, 8); else ctx.rect(bx, by, bw, bh);
      ctx.fill();
      ctx.stroke();
      ctx.fillStyle = "#d7e2ec";
      ctx.textBaseline = "top";
      ctx.textAlign = "center";
      rows.forEach((r, k) => ctx.fillText(r, cssW / 2, by + 6 + k * lh + 1));
      ctx.textAlign = "start";
    }

    return { resize, draw };
  }

  return { create };
})();

// ===========================================================================
// 3. App — DOM, HUD, sliders, mouse, Object Edit Mode, main loop (skipped in Node)
// ===========================================================================
if (typeof document !== "undefined") {
  (function () {
    const solver = FlipSolver;
    const canvas = document.getElementById("view");
    const renderer = Canvas2DRenderer.create(canvas);
    const hudN = document.getElementById("hud-n");
    const hudFps = document.getElementById("hud-fps");
    const btnPause = document.getElementById("btn-pause");
    const panel = document.getElementById("panel");
    const MAX_FRAMES_PER_TICK = 2; // catch-up limit so a slow tab cannot spiral

    let view = computeLayout(1, 1);
    let paused = false;
    let mouseDown = false;
    let lastMx = 0;
    let lastMy = 0;
    let lastMouseStamp = 0;
    let mvx = 0; // px/s
    let mvy = 0;
    let msx = 0; // px
    let msy = 0;
    let accTime = 0;
    let lastStamp = 0;
    let fpsEma = 60;
    let fpsFrames = 0;
    let fpsStamp = 0;

    const toWorldX = (sx) => (sx - view.offX) / view.scale + view.originX;
    const toWorldY = (sy) => (sy - view.offY) / view.scale + view.originY;

    // ---- Object Edit Mode state ----------------------------------------------
    const TOOLS = ["circle", "box", "ramp", "triangle", "freehand"];
    const TOOL_LABEL = { circle: "Circle", box: "Box", ramp: "Ramp", triangle: "Triangle", freehand: "Freehand" };
    const ROT_STEP = Math.PI / 12; // 15° per [ / ] press (Shift: 3°)
    const edit = {
      on: false,
      prevPaused: false, // pause state before entering, restored on exit
      tool: "circle",
      pending: null, // shape (solver units) shown as a ghost, not yet in the solver
      editingId: -1, // committed obstacle that `pending` is a moved/rotated copy of
      selectedId: -1,
      drag: null, // { kind: "place" | "move" | "rotate" | "grab" | "freehand", ... }
      stroke: null, // freehand points while drawing
      status: "",
      thin: null, // { shape, cells } of the last shape rejected by the thickness rule (red overlay)
    };

    function setupCanvas() {
      const dpr = Math.max(1, window.devicePixelRatio || 1);
      const w = window.innerWidth;
      const hh = window.innerHeight;
      renderer.resize(w, hh, dpr);
      view = computeLayout(w, hh);
      solver.init(view.cellsX, view.cellsY, view.cellSize);
    }

    function syncPointer() {
      solver.setPointer(toWorldX(msx), toWorldY(msy), mvx / view.scale, mvy / view.scale);
    }

    function tick(stamp) {
      requestAnimationFrame(tick);
      if (!lastStamp) lastStamp = stamp;
      const frameDt = Math.min(0.05, (stamp - lastStamp) / 1000); // clamp the raw gap
      lastStamp = stamp;
      if (lastMouseStamp && stamp - lastMouseStamp > 50) {
        mvx = 0; // stale pointer -> no throw velocity
        mvy = 0;
      }
      syncPointer();

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
        let frames = 0;
        while (accTime >= solver.FIXED_DT && frames < MAX_FRAMES_PER_TICK) {
          solver.step();
          accTime -= solver.FIXED_DT;
          frames++;
        }
        if (frames === MAX_FRAMES_PER_TICK) accTime = 0;
      } else if (solver.holding) {
        solver.pinHeld(); // follow the cursor while paused
      }
      renderer.draw(solver, view, edit.on ? editOverlay() : null);
      hudN.textContent = String(solver.count);
    }

    function pointerFromEvent(e) {
      const rect = canvas.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;
      const now = performance.now();
      if (lastMouseStamp > 0) {
        const dt = (now - lastMouseStamp) / 1000;
        if (dt > 1e-4 && dt < 0.1) {
          // light smoothing so one jittery event does not decide the throw
          mvx = 0.5 * mvx + 0.5 * ((x - lastMx) / dt);
          mvy = 0.5 * mvy + 0.5 * ((y - lastMy) / dt);
        }
      }
      lastMx = x;
      lastMy = y;
      lastMouseStamp = now;
      msx = x;
      msy = y;
      syncPointer();
    }

    function setPaused(next) {
      paused = next;
      btnPause.textContent = paused ? "Resume" : "Pause";
      updateEditPanel();
    }

    function doReset() {
      accTime = 0;
      solver.reset("dam"); // keeps obstacles; no particles spawn inside them
    }

    // ---- Object Edit Mode ------------------------------------------------------
    const el = (id) => document.getElementById(id);

    function rotateHandle(shape) {
      if (!shape || shape.type === "circle") return null;
      const r = ObstacleShapes.boundingRadius(shape) + 26 / view.scale;
      const a = shape.angle || 0;
      return { x: shape.cx + Math.sin(a) * r, y: shape.cy - Math.cos(a) * r };
    }

    function editHint() {
      const t = TOOL_LABEL[edit.tool];
      if (window.innerWidth < 600) {
        return "EDIT · " + t + " · " + (edit.tool === "freehand" ? "draw a loop" : "drag to place") + " · drag object to move · handle rotates · Confirm / Cancel / Done in the panel";
      }
      const place = edit.tool === "freehand" ? "hold & draw, release to close" : "drag to place/size (tap = default size)";
      return "EDIT MODE · " + t + " (1-5) · " + place + " · drag an object to move · [ ] / wheel / handle rotate · Enter confirm · Esc cancel · Del delete · E done";
    }

    function editOverlay() {
      return {
        editMode: true,
        pending: edit.pending,
        editingId: edit.editingId,
        selectedId: edit.selectedId,
        stroke: edit.stroke,
        handle: edit.pending ? rotateHandle(edit.pending) : null,
        hint: editHint(),
        status: edit.status,
        thin: edit.thin && edit.thin.shape === edit.pending ? edit.thin.cells : null,
      };
    }

    function setStatus(msg) {
      edit.status = msg || "";
      const st = el("edit-status");
      if (st) st.textContent = edit.status;
    }

    function updateEditPanel() {
      if (!panel) return;
      panel.classList.toggle("editing", edit.on);
      const ep = el("edit-panel");
      if (ep) ep.hidden = !edit.on;
      const b = el("btn-edit");
      if (b) b.textContent = edit.on ? "Done editing (E)" : "Edit objects (E)";
      for (const tb of document.querySelectorAll("[data-tool]")) tb.classList.toggle("active", tb.dataset.tool === edit.tool);
      const pz = el("edit-paused");
      if (pz) pz.textContent = paused ? "sim paused" : "sim running";
      const has = !!edit.pending;
      const cf = el("btn-confirm");
      if (cf) cf.disabled = !has;
      const cc = el("btn-cancel");
      if (cc) cc.disabled = !has && !edit.stroke;
      const del = el("btn-delete");
      if (del) del.disabled = !has && edit.selectedId < 0;
      const n = el("edit-count");
      if (n) n.textContent = String(solver.obstacles.length);
      canvas.style.cursor = edit.on ? "default" : "";
    }

    function setEditMode(on) {
      if (on === edit.on) return;
      if (on) {
        if (mouseDown) releaseDrop(); // never leave a drop pinned to the pointer
        edit.on = true;
        edit.prevPaused = paused;
        setPaused(true);
        setStatus("");
      } else {
        // Done = accept the shape on screen (Esc first to discard); a shape rejected by the
        // thickness rule keeps edit mode open so the message and the marked cells stay visible
        if (edit.pending && !confirmPending()) return;
        edit.drag = null;
        edit.stroke = null;
        edit.selectedId = -1;
        edit.on = false;
        setPaused(edit.prevPaused);
      }
      updateEditPanel();
    }

    function setTool(t) {
      if (TOOLS.indexOf(t) < 0) return;
      edit.tool = t;
      updateEditPanel();
    }

    /** Shape for the current tool from a drag (x0,y0)->(x1,y1) in metres; isTap = default size. */
    function toolShape(x0, y0, x1, y1, isTap) {
      const h = solver.cellSize;
      const minS = 3 * h; // keep shapes a few cells across so they show up on the grid
      const span = (a, b, def) => {
        if (isTap) return [a - def / 2, a + def / 2];
        const d = b - a;
        const m = Math.max(Math.abs(d), minS);
        return d >= 0 ? [a, a + m] : [a - m, a];
      };
      if (edit.tool === "circle") {
        const r = isTap ? 0.12 : Math.max(Math.hypot(x1 - x0, y1 - y0), minS / 2);
        return ObstacleShapes.circle(x0, y0, r);
      }
      if (edit.tool === "box") {
        const [xa, xb] = span(x0, x1, 0.36), [ya, yb] = span(y0, y1, 0.18);
        return ObstacleShapes.box(0.5 * (xa + xb), 0.5 * (ya + yb), xb - xa, yb - ya, 0);
      }
      if (edit.tool === "ramp") {
        if (isTap) return ObstacleShapes.ramp(x0 + 0.3, y0 - 0.15, x0 - 0.3, y0 + 0.15); // rises to the right
        const [xa, xb] = span(x0, x1, 0), [ya, yb] = span(y0, y1, 0);
        // the drag start is the ramp's high end if it is above the release point
        const sx = x1 >= x0 ? xa : xb, ex = x1 >= x0 ? xb : xa, sy = y1 >= y0 ? ya : yb, ey = y1 >= y0 ? yb : ya;
        return ObstacleShapes.ramp(sx, sy, ex, ey);
      }
      if (edit.tool === "triangle") {
        const [xa, xb] = span(x0, x1, 0.4), [ya, yb] = span(y0, y1, 0.32);
        return ObstacleShapes.triangle(xa, ya, xb, yb);
      }
      return null;
    }

    function startEditFromObstacle(id) {
      const o = solver.obstacles.find((q) => q.id === id);
      if (!o) return false;
      const c = ObstacleShapes.clone(o);
      delete c.id;
      edit.pending = c;
      edit.editingId = id;
      edit.selectedId = id;
      return true;
    }

    /**
     * Thickness rule (item 1): every part of a shape must be >= 2 grid cells thick
     * (ObstacleShapes.thicknessCheck). Polygons (ramp / triangle / freehand) whose only
     * problem is a sharp corner get that corner clipped (cut 2 cells, widened in half-cell
     * steps up to 6 cells until the rule holds); anything else that fails is rejected with
     * a message and the failing cells marked in red, and stays pending to be fixed.
     * Returns { shape, clipped } or { reject: message, cells }.
     */
    function applyThicknessRule(shape, label) {
      const h = solver.cellSize, cm = (2 * h * 100).toFixed(1);
      let chk;
      try { chk = solver.checkObstacleThickness(shape); } catch (err) { return { reject: "Could not place shape: " + err.message, cells: [] }; }
      if (chk.ok) return { shape, clipped: 0 };
      if (shape.type === "polygon" && chk.cells > 0) {
        for (let cut = 2 * h; cut <= 6 * h + 1e-9; cut += 0.5 * h) {
          const c = ObstacleShapes.clipAcuteTips(shape, cut);
          if (c.clipped && solver.checkObstacleThickness(c.shape).ok) return { shape: c.shape, clipped: c.clipped };
        }
      }
      if (chk.cells === 0) return { reject: "Not placed: this " + label + " is smaller than one grid cell, so the water would not see it. Make it at least 2 cells (" + cm + " cm) thick.", cells: [] };
      return {
        reject: "Not placed: parts of this " + label + " (red) are thinner than 2 grid cells (" + cm + " cm). Make it thicker, or Esc to cancel.",
        cells: chk.thin,
      };
    }

    function confirmPending() {
      if (!edit.pending) {
        setStatus("Nothing to confirm: drag on the canvas to place a shape.");
        return true;
      }
      const label = edit.pending.type === "polygon" ? edit.pending.kind : edit.pending.type === "rect" ? "box" : "circle";
      const rule = applyThicknessRule(edit.pending, label);
      if (rule.reject) {
        edit.thin = { shape: edit.pending, cells: rule.cells };
        setStatus(rule.reject);
        updateEditPanel();
        return false;
      }
      let id;
      try {
        if (edit.editingId > 0 && solver.updateObstacle(edit.editingId, rule.shape)) id = edit.editingId;
        else id = solver.addObstacle(rule.shape);
      } catch (err) {
        setStatus("Could not place shape: " + err.message);
        return false;
      }
      const removed = solver.lastRemovedParticles;
      const clip = rule.clipped ? " · sharp tip" + (rule.clipped > 1 ? "s" : "") + " clipped so every part is at least 2 grid cells thick" : "";
      setStatus((edit.editingId > 0 ? "Moved " : "Placed ") + label + clip + (removed ? " · removed " + removed + " particles that were inside it" : ""));
      edit.thin = null;
      edit.pending = null;
      edit.editingId = -1;
      edit.selectedId = id;
      updateEditPanel();
      return true;
    }

    function cancelPending() {
      const had = edit.pending || edit.stroke;
      edit.pending = null;
      edit.editingId = -1;
      edit.stroke = null;
      edit.drag = null;
      if (!had) edit.selectedId = -1; // second Esc clears the selection
      setStatus(had ? "Cancelled" : "");
      updateEditPanel();
    }

    function deleteSelected() {
      if (edit.editingId > 0) {
        solver.removeObstacle(edit.editingId);
        edit.pending = null;
        edit.editingId = -1;
        edit.selectedId = -1;
        setStatus("Deleted object");
      } else if (edit.pending) {
        cancelPending();
      } else if (edit.selectedId > 0) {
        solver.removeObstacle(edit.selectedId);
        edit.selectedId = -1;
        setStatus("Deleted object");
      }
      updateEditPanel();
    }

    function rotateBy(da) {
      if (!edit.pending && edit.selectedId > 0) startEditFromObstacle(edit.selectedId);
      if (!edit.pending) {
        setStatus("Select or place an object to rotate it.");
        return;
      }
      edit.pending = ObstacleShapes.rotate(edit.pending, da);
      updateEditPanel();
    }

    function clearAll() {
      solver.clearObstacles();
      edit.pending = null;
      edit.editingId = -1;
      edit.selectedId = -1;
      edit.stroke = null;
      setStatus("All objects removed");
      updateEditPanel();
    }

    function eventWorld(e) {
      const rect = canvas.getBoundingClientRect();
      const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
      return { sx, sy, x: toWorldX(sx), y: toWorldY(sy) };
    }

    function editPointerDown(e) {
      const p = eventWorld(e);
      const touch = e.pointerType === "touch" || e.pointerType === "pen";
      const tol = (touch ? 14 : 5) / view.scale;
      const hnd = edit.pending ? rotateHandle(edit.pending) : null;
      if (hnd && Math.hypot(p.x - hnd.x, p.y - hnd.y) * view.scale <= (touch ? 24 : 12)) {
        edit.drag = { kind: "rotate" };
      } else if (edit.pending && ObstacleShapes.sdf(edit.pending, p.x, p.y) <= tol) {
        edit.drag = { kind: "move", lx: p.x, ly: p.y };
      } else {
        const id = solver.hitTest(p.x, p.y, tol);
        if (id > 0) {
          if (edit.pending && edit.editingId !== id) {
            edit.pending = null; // an unconfirmed ghost is replaced by the new selection
            edit.editingId = -1;
          }
          edit.selectedId = id;
          edit.drag = { kind: "grab", id, sx0: p.sx, sy0: p.sy, lx: p.x, ly: p.y };
        } else {
          edit.pending = null;
          edit.editingId = -1;
          edit.selectedId = -1;
          if (edit.tool === "freehand") {
            edit.stroke = [[p.x, p.y]];
            edit.drag = { kind: "freehand" };
          } else {
            edit.drag = { kind: "place", x0: p.x, y0: p.y, sx0: p.sx, sy0: p.sy };
            edit.pending = toolShape(p.x, p.y, p.x, p.y, true);
          }
        }
      }
      setStatus(edit.status && edit.drag && edit.drag.kind === "place" ? "" : edit.status);
      updateEditPanel();
    }

    function editPointerMove(e) {
      const p = eventWorld(e);
      const d = edit.drag;
      if (!d) {
        const tol = 5 / view.scale;
        const hnd = edit.pending ? rotateHandle(edit.pending) : null;
        let cur = "crosshair";
        if (hnd && Math.hypot(p.x - hnd.x, p.y - hnd.y) * view.scale <= 12) cur = "grab";
        else if ((edit.pending && ObstacleShapes.sdf(edit.pending, p.x, p.y) <= tol) || solver.hitTest(p.x, p.y, tol) > 0) cur = "move";
        canvas.style.cursor = cur;
        return;
      }
      if (d.kind === "rotate" && edit.pending) {
        let a = Math.atan2(p.y - edit.pending.cy, p.x - edit.pending.cx) + Math.PI / 2;
        if (e.shiftKey) a = Math.round(a / ROT_STEP) * ROT_STEP;
        edit.pending = ObstacleShapes.rotate(edit.pending, a - (edit.pending.angle || 0));
      } else if (d.kind === "move" && edit.pending) {
        edit.pending = ObstacleShapes.translate(edit.pending, p.x - d.lx, p.y - d.ly);
        d.lx = p.x;
        d.ly = p.y;
      } else if (d.kind === "grab") {
        if (Math.hypot(p.sx - d.sx0, p.sy - d.sy0) > 4 && startEditFromObstacle(d.id)) {
          edit.pending = ObstacleShapes.translate(edit.pending, p.x - d.lx, p.y - d.ly);
          edit.drag = { kind: "move", lx: p.x, ly: p.y };
          updateEditPanel();
        }
      } else if (d.kind === "place") {
        const tap = Math.hypot(p.sx - d.sx0, p.sy - d.sy0) < 6;
        edit.pending = toolShape(d.x0, d.y0, p.x, p.y, tap);
      } else if (d.kind === "freehand" && edit.stroke) {
        const q = edit.stroke[edit.stroke.length - 1];
        if (Math.hypot(p.x - q[0], p.y - q[1]) > 0.35 * solver.cellSize) edit.stroke.push([p.x, p.y]);
      }
    }

    function editPointerUp() {
      const d = edit.drag;
      edit.drag = null;
      if (d && d.kind === "freehand" && edit.stroke) {
        const res = ObstacleShapes.freehand(edit.stroke, solver.cellSize);
        edit.stroke = null;
        if (res.shape) {
          edit.pending = res.shape;
          setStatus(res.note ? "Freehand: " + res.note + ". Enter to confirm." : "Freehand shape ready. Enter to confirm.");
        } else setStatus("Freehand rejected: " + res.error);
      } else if (d && (d.kind === "place" || d.kind === "move" || d.kind === "rotate") && edit.pending) {
        setStatus("Enter / Confirm to commit · Esc / Cancel to discard");
      }
      updateEditPanel();
    }

    function releaseDrop() {
      if (!mouseDown) return;
      mouseDown = false;
      syncPointer();
      solver.releaseDrop();
    }

    function wireUi() {
      const slN = el("sl-n");
      const slG = el("sl-g");
      const slFlip = el("sl-flip");
      const slIters = el("sl-iters");
      const slMu = el("sl-mu");
      const p = solver.params;

      function sync() {
        p.fill = parseInt(slN.value, 10) / 100;
        p.gravity = parseFloat(slG.value);
        p.flipRatio = parseFloat(slFlip.value);
        // PCG stops on the residual; the slider is only its safety cap (GS: fixed count)
        if (p.pressureSolver === "gs") p.pressureIters = parseInt(slIters.value, 10);
        else p.pcgMaxIters = parseInt(slIters.value, 10);
        p.viscosity = parseFloat(slMu.value);
        el("val-n").textContent = slN.value + "%";
        el("val-g").textContent = p.gravity.toFixed(2);
        el("val-flip").textContent = p.flipRatio.toFixed(2);
        el("val-iters").textContent = slIters.value;
        el("val-mu").textContent = p.viscosity.toFixed(2);
      }
      slN.addEventListener("input", function () {
        sync();
        doReset();
      });
      for (const s of [slG, slFlip, slIters, slMu]) s.addEventListener("input", sync);
      sync();

      btnPause.addEventListener("click", () => setPaused(!paused));
      el("btn-reset").addEventListener("click", doReset);
      el("btn-dam").addEventListener("click", doReset);

      // Object Edit Mode buttons (touch: Chromebook / iPhone)
      const on = (id, fn) => { const b = el(id); if (b) b.addEventListener("click", fn); };
      on("btn-edit", () => setEditMode(!edit.on));
      on("btn-done", () => setEditMode(false));
      on("btn-confirm", confirmPending);
      on("btn-cancel", cancelPending);
      on("btn-delete", deleteSelected);
      on("btn-rot-l", () => rotateBy(-ROT_STEP));
      on("btn-rot-r", () => rotateBy(ROT_STEP));
      on("btn-clear-obj", clearAll);
      for (const tb of document.querySelectorAll("[data-tool]")) tb.addEventListener("click", () => setTool(tb.dataset.tool));

      canvas.addEventListener("pointerdown", function (e) {
        if (e.button !== 0) return;
        pointerFromEvent(e);
        mvx = mvy = 0;
        try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* synthetic / stale pointer */ }
        if (edit.on) {
          editPointerDown(e); // no water drops while editing objects
          return;
        }
        mouseDown = true;
        syncPointer();
        if (!solver.holding) solver.spawnDrop(toWorldX(msx), toWorldY(msy));
      });
      canvas.addEventListener("pointermove", function (e) {
        pointerFromEvent(e);
        if (edit.on) editPointerMove(e);
      });
      const release = function () {
        if (edit.on && edit.drag) {
          editPointerUp();
          return;
        }
        releaseDrop();
      };
      canvas.addEventListener("pointerup", release);
      canvas.addEventListener("pointercancel", release);
      canvas.addEventListener("wheel", function (e) {
        if (!edit.on) return;
        e.preventDefault();
        if (!edit.pending && edit.selectedId < 0) {
          const p = eventWorld(e);
          const id = solver.hitTest(p.x, p.y, 5 / view.scale);
          if (id > 0) edit.selectedId = id;
        }
        rotateBy((e.deltaY > 0 ? 1 : -1) * (e.shiftKey ? ROT_STEP / 5 : ROT_STEP / 3));
      }, { passive: false });

      window.addEventListener("keydown", function (e) {
        if (e.ctrlKey || e.metaKey || e.altKey) return;
        const code = e.code;
        if (code === "Space") {
          e.preventDefault();
          setPaused(!paused);
        } else if (code === "KeyR") {
          doReset();
        } else if (code === "KeyE") {
          if (!e.repeat) setEditMode(!edit.on);
        } else if (edit.on) {
          const digit = /^(Digit|Numpad)([1-5])$/.exec(code);
          if (digit) setTool(TOOLS[parseInt(digit[2], 10) - 1]);
          else if (code === "BracketLeft") rotateBy(-(e.shiftKey ? ROT_STEP / 5 : ROT_STEP));
          else if (code === "BracketRight") rotateBy(e.shiftKey ? ROT_STEP / 5 : ROT_STEP);
          else if (code === "Delete" || code === "Backspace") {
            e.preventDefault();
            deleteSelected();
          } else if (code === "Enter" || code === "NumpadEnter") {
            e.preventDefault();
            confirmPending();
          } else if (code === "Escape") cancelPending();
        }
      });
      let resizeTimer = 0;
      window.addEventListener("resize", function () {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(function () {
          setupCanvas(); // obstacles (metres) are kept and re-rasterized
          doReset();
        }, 150);
      });
      updateEditPanel();
    }

    setupCanvas();
    wireUi();
    const query = new URLSearchParams(window.location.search);
    // Dev option: index.html?demoObstacles=1 preloads a circle and a ramp in the
    // dam-break path (screenshots / quick look). No effect without the parameter.
    if (query.get("demoObstacles") === "1") for (const s of demoObstacleShapes(solver)) solver.addObstacle(s);
    doReset();
    // Dev option: index.html?warm=3 pre-runs 3 s of simulation before the first
    // paint (handy for screenshots / thumbnails). No effect without the parameter.
    const warm = parseFloat(query.get("warm"));
    if (warm > 0) for (let f = 0; f < Math.min(warm, 60) / solver.FIXED_DT; f++) solver.step();
    requestAnimationFrame(tick);
  })();
}

// Node / headless harness hook (ignored by the browser).
if (typeof module !== "undefined" && module.exports) {
  module.exports = { FlipSolver, computeLayout, Canvas2DRenderer, ObstacleShapes, demoObstacleShapes };
}
