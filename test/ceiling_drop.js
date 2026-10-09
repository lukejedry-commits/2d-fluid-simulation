// Ceiling-drop release timing: 2g separating vs legacy walls.
// Contact = first 0.02 s sample with any particle centre within one cell of the lid.
// Release = first later sample with none. t_release = release - contact.
// Cases match the earlier diagnostic (fresh process per boundary).
const { FlipSolver: S } = require("../sim.js");
const h = 3 / 128;
const CASES = [
  { name: "center", x: 1.5, y: 0.6, vx: 0, vy: -8 },
  { name: "left", x: 0.4, y: 0.6, vx: -3, vy: -8 },
  { name: "right", x: 2.6, y: 0.6, vx: 5, vy: -6 },
];
function run(boundary) {
  const out = [];
  for (const C of CASES) {
    S.init(128, 72, h);
    S.params.boundary = boundary;
    S.params.separation = "neutral";
    S.params.damping = 0.15; // browser default; drop is a UI stress case
    S.reset("pool");
    S.spawnDrop(C.x, C.y);
    S.setPointer(C.x, C.y, C.vx, C.vy);
    for (let k = 0; k < 3; k++) S.step();
    S.releaseDrop();
    const lid = h; // one cell of the lid (y down from top wall face at y=h... particle near ceiling has py small)
    // Ceiling solid occupies j=0; first fluid row centres near y ~= 1.5 h. "within one cell of the lid"
    // = particle centre within one cell size of the inner lid face at y = h, i.e. py < h + h = 2h.
    const near = () => { for (let i = 0; i < S.count; i++) if (S.py[i] < 2 * h) return true; return false; };
    let contact = null, release = null;
    const dt = 0.02;
    let tSample = 0;
    while (S.simTime < 4) {
      S.step();
      while (tSample <= S.simTime + 1e-12) {
        if (near()) { if (contact === null) contact = tSample; release = null; }
        else if (contact !== null && release === null) release = tSample;
        tSample = +(tSample + dt).toFixed(6);
      }
    }
    const tRel = contact !== null && release !== null ? release - contact : (contact !== null ? Infinity : 0);
    out.push({ ...C, boundary, contact, release, t_release: tRel, DROP_RADIUS_CELLS: S.DROP_RADIUS_CELLS, R_m: S.DROP_RADIUS_CELLS * h });
  }
  return out;
}
function damLid(boundary) {
  S.init(128, 72, h);
  S.params.boundary = boundary; S.params.separation = "neutral"; S.params.damping = 0;
  S.reset("dam");
  const near = () => { for (let i = 0; i < S.count; i++) if (S.py[i] < 2 * h) return true; return false; };
  let contact = null, release = null; const dt = 0.02; let tSample = 0;
  while (S.simTime < 6) {
    S.step();
    while (tSample <= S.simTime + 1e-12) {
      if (near()) { if (contact === null) contact = tSample; release = null; }
      else if (contact !== null && release === null) release = tSample;
      tSample = +(tSample + dt).toFixed(6);
    }
  }
  return { boundary, contact, release, t_release: contact != null && release != null ? release - contact : (contact != null ? Infinity : 0) };
}
const sep = run("separating");
const leg = run("legacy");
const damS = damLid("separating");
const damL = damLid("legacy");
const rows = sep.map((s, i) => {
  const l = leg[i];
  const pass = s.t_release <= l.t_release + 0.10 && s.t_release <= 0.55;
  return { case: s.name, start: [s.x, s.y], v: [s.vx, s.vy], R_m: s.R_m, legacy_t: l.t_release, sep_t: s.t_release, delta: s.t_release - l.t_release, pass,
    legacy_contact: l.contact, sep_contact: s.contact };
});
const damPass = damS.t_release <= damL.t_release + 0.10 && damS.t_release <= 0.55;
const out = { dropParams: { DROP_RADIUS_CELLS: 5, R_m: 5 * h, holdStepsBeforeRelease: 3, sampleDt: 0.02, contactDef: "any particle centre within one cell (0.0234375 m) of the lid (py < 2h)", damping: 0.15 },
  drops: rows, allDropsPass: rows.every((r) => r.pass),
  damBreak: { legacy: damL, separating: damS, pass: damPass, damping: 0 } };
console.log(JSON.stringify(out, null, 2));
require("fs").writeFileSync(require("path").join(__dirname, "../validation/ceiling_drop_2g.json"), JSON.stringify(out, null, 2));
