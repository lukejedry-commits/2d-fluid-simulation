// Browser smoke test for index.html + Object Edit Mode, driven over the Chrome
// DevTools Protocol with REAL input events (mouse, keyboard, touch). No npm deps.
// Usage: node --experimental-websocket test/browser_smoke.js [chromePath]
// Writes obstacles_preview.png (?demoObstacles=1&warm=0.6), obstacles_editmode.png and
// obstacles_editmode_phone.png (390x844 touch viewport). --preview-only: just the preview.
// Not referenced by index.html.
const { spawn } = require("child_process");
const fs = require("fs");
const os = require("os");
const path = require("path");

const ROOT = path.join(__dirname, "..");
const CHROME = (process.argv[2] && !process.argv[2].startsWith("--") ? process.argv[2] : "") || ["google-chrome", "chromium", "chromium-browser"].find((c) => {
  try { require("child_process").execSync("command -v " + c, { stdio: "ignore" }); return true; } catch (e) { return false; }
});
if (typeof WebSocket === "undefined") {
  console.error("needs WebSocket: run with node --experimental-websocket (Node 20) or Node >= 22");
  process.exit(2);
}
if (!CHROME) { console.error("no chrome/chromium found"); process.exit(2); }
const PREVIEW_ONLY = process.argv.includes("--preview-only");
const PREVIEW_WARM = process.env.PREVIEW_WARM || "0.6";
const PREVIEW_OUT = process.env.PREVIEW_OUT || path.join(ROOT, "obstacles_preview.png");
const W = 1280, H = 800, PORT = 9333 + Math.floor(Math.random() * 500);
const profile = fs.mkdtempSync(path.join(os.tmpdir(), "flip-smoke-"));
const chrome = spawn(CHROME, ["--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars", "--no-first-run",
  "--remote-debugging-port=" + PORT, "--user-data-dir=" + profile, "--window-size=" + W + "," + H, "about:blank"], { stdio: "ignore" });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

let ws, nextId = 1;
const pending = new Map();
const listeners = [];
function send(method, params) {
  return new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject, method });
    ws.send(JSON.stringify({ id, method, params: params || {} }));
  });
}
const waitEvent = (name, ms = 15000) => new Promise((resolve, reject) => {
  const t = setTimeout(() => reject(new Error("timeout waiting for " + name)), ms);
  listeners.push({ name, fn: (p) => { clearTimeout(t); resolve(p); }, once: true });
});
async function evaluate(expr) {
  const r = await send("Runtime.evaluate", { expression: expr, returnByValue: true });
  if (r.exceptionDetails) throw new Error("eval failed: " + expr + " → " + JSON.stringify(r.exceptionDetails.exception && r.exceptionDetails.exception.description));
  return r.result.value;
}
const mouse = (type, x, y, extra) => send("Input.dispatchMouseEvent", Object.assign({ type, x, y, button: "left", buttons: type === "mouseReleased" ? 0 : 1, clickCount: 1 }, extra || {}));
async function drag(x0, y0, x1, y1, steps = 8) {
  await mouse("mousePressed", x0, y0);
  for (let k = 1; k <= steps; k++) await mouse("mouseMoved", x0 + ((x1 - x0) * k) / steps, y0 + ((y1 - y0) * k) / steps);
  await mouse("mouseReleased", x1, y1);
  await sleep(50);
}
const click = async (x, y) => { await mouse("mousePressed", x, y); await mouse("mouseReleased", x, y); await sleep(50); };
const KEYS = {
  e: ["KeyE", "e", 69], Enter: ["Enter", "Enter", 13], Escape: ["Escape", "Escape", 27], Delete: ["Delete", "Delete", 46],
  "]": ["BracketRight", "]", 221], "[": ["BracketLeft", "[", 219], 1: ["Digit1", "1", 49], 2: ["Digit2", "2", 50],
  r: ["KeyR", "r", 82], 3: ["Digit3", "3", 51], 4: ["Digit4", "4", 52], 5: ["Digit5", "5", 53], Space: ["Space", " ", 32],
};
async function key(k) {
  const [code, keyName, vk] = KEYS[k];
  const text = keyName.length === 1 ? keyName : undefined;
  await send("Input.dispatchKeyEvent", { type: "keyDown", code, key: keyName, windowsVirtualKeyCode: vk, text });
  await send("Input.dispatchKeyEvent", { type: "keyUp", code, key: keyName, windowsVirtualKeyCode: vk });
  await sleep(30);
}
async function touchDrag(x0, y0, x1, y1, steps = 8) {
  await send("Input.dispatchTouchEvent", { type: "touchStart", touchPoints: [{ x: x0, y: y0, id: 1 }] });
  for (let k = 1; k <= steps; k++) await send("Input.dispatchTouchEvent", { type: "touchMove", touchPoints: [{ x: x0 + ((x1 - x0) * k) / steps, y: y0 + ((y1 - y0) * k) / steps, id: 1 }] });
  await send("Input.dispatchTouchEvent", { type: "touchEnd", touchPoints: [] });
  await sleep(300); // a new touch sequence right after touchEnd can get folded into the gesture
}
async function screenshot(file) {
  const r = await send("Page.captureScreenshot", { format: "png" });
  fs.writeFileSync(file, Buffer.from(r.data, "base64"));
}

const results = [];
function check(name, ok, info) {
  results.push({ name, ok: !!ok });
  console.log((ok ? "ok   " : "FAIL ") + name + (info !== undefined ? "  (" + (typeof info === "string" ? info : JSON.stringify(info)) + ")" : ""));
}

async function main() {
  let target;
  for (let t = 0; t < 100 && !target; t++) {
    try {
      const list = await (await fetch("http://127.0.0.1:" + PORT + "/json/list")).json();
      target = list.find((x) => x.type === "page");
    } catch (e) { await sleep(100); }
  }
  if (!target) throw new Error("chrome did not start");
  ws = new WebSocket(target.webSocketDebuggerUrl);
  await new Promise((r, j) => { ws.onopen = r; ws.onerror = j; });
  const errors = [];
  ws.onmessage = (m) => {
    const msg = JSON.parse(m.data);
    if (msg.id && pending.has(msg.id)) {
      const p = pending.get(msg.id);
      pending.delete(msg.id);
      if (msg.error) p.reject(new Error(p.method + ": " + msg.error.message)); else p.resolve(msg.result);
    } else if (msg.method) {
      if (msg.method === "Runtime.exceptionThrown") errors.push(msg.params.exceptionDetails.exception ? msg.params.exceptionDetails.exception.description : msg.params.exceptionDetails.text);
      if (msg.method === "Runtime.consoleAPICalled" && msg.params.type === "error") errors.push(msg.params.args.map((a) => a.value || a.description).join(" "));
      for (let k = listeners.length - 1; k >= 0; k--) if (listeners[k].name === msg.method) { const l = listeners.splice(k, 1)[0]; l.fn(msg.params); }
    }
  };
  await send("Page.enable");
  await send("Runtime.enable");
  await send("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  await send("Emulation.setTouchEmulationEnabled", { enabled: true, maxTouchPoints: 5 });

  const url = "file://" + path.join(ROOT, "index.html");
  let loaded;
  if (!PREVIEW_ONLY) await uiChecks();
  async function uiChecks() {
  loaded = waitEvent("Page.loadEventFired");
  await send("Page.navigate", { url });
  await loaded;
  await sleep(600);

  check("page loads, solver running", await evaluate("FlipSolver.count > 0 && FlipSolver.obstacles.length === 0 && document.getElementById('btn-pause').textContent === 'Pause'"));
  // --- normal mode: mouse drop still works
  let n0 = await evaluate("FlipSolver.count");
  await mouse("mousePressed", 900, 200);
  await sleep(80);
  check("normal mode: mouse hold spawns a held drop", await evaluate("FlipSolver.holding === true"));
  await mouse("mouseReleased", 900, 200);
  await sleep(50);
  check("normal mode: release frees the drop", await evaluate("FlipSolver.holding === false && FlipSolver.count > " + n0));

  // --- enter edit mode
  await key("e");
  check("E enters edit mode and pauses", await evaluate("document.getElementById('panel').classList.contains('editing') && !document.getElementById('edit-panel').hidden && document.getElementById('btn-pause').textContent === 'Resume'"));
  const scale = await evaluate("window.innerWidth / 3.0"); // px per metre (computeLayout: canvas width = 3 m)
  const st0 = await evaluate("FlipSolver.simTime");
  await sleep(200);
  check("edit mode: simulation paused", (await evaluate("FlipSolver.simTime")) === st0);

  // circle (tool 1): drag 60 px radius, drop must NOT spawn water
  n0 = await evaluate("FlipSolver.count");
  await key("1");
  await drag(700, 330, 760, 330);
  check("edit mode: dragging does not spawn a water drop", await evaluate("FlipSolver.holding === false && FlipSolver.count === " + n0));
  check("pending circle is not committed before Enter", await evaluate("FlipSolver.obstacles.length === 0"));
  await key("Enter");
  let obs = await evaluate("FlipSolver.obstacles");
  check("Enter commits circle", obs.length === 1 && obs[0].type === "circle" && Math.abs(obs[0].r * scale - 60) < 2, obs[0] && { r_px: +(obs[0].r * scale).toFixed(1) });
  const circleId = obs[0] && obs[0].id;

  // box (tool 2)
  await key("2");
  await drag(380, 520, 530, 590); // (keep clear of the panel at the top left)
  await key("Enter");
  obs = await evaluate("FlipSolver.obstacles");
  const box = obs.find((o) => o.type === "rect");
  check("box placed by corner drag", obs.length === 2 && box && Math.abs(box.w * scale - 150) < 2 && Math.abs(box.h * scale - 70) < 2, box && { w_px: +(box.w * scale).toFixed(1), h_px: +(box.h * scale).toFixed(1) });
  const removedTxt = await evaluate("document.getElementById('edit-status').textContent");
  check("status line reports placement", /Placed box/.test(removedTxt), removedTxt);

  // ramp (tool 3) then Esc cancels
  await key("3");
  await drag(1150, 600, 950, 700);
  check("ramp ghost pending", (await evaluate("FlipSolver.obstacles.length")) === 2);
  await key("Escape");
  await key("Enter"); // nothing pending -> no-op
  check("Esc discards pending ramp", (await evaluate("FlipSolver.obstacles.length")) === 2);

  // ramp placed for real
  await drag(1150, 600, 950, 700);
  await key("Enter");
  obs = await evaluate("FlipSolver.obstacles");
  check("ramp committed", obs.length === 3 && obs[2].type === "polygon" && obs[2].kind === "ramp");
  const rampTxt = await evaluate("document.getElementById('edit-status').textContent");
  check("mid-air ramp: sharp tips clipped to >= 2 cells, still a ramp", /clipped/.test(rampTxt) && obs[2].points.length > 3 && (await evaluate("FlipSolver.checkObstacleThickness(FlipSolver.obstacles[2]).ok")), rampTxt);

  // freehand (tool 5): a closed loop
  await key("5");
  const fx = 1000, fy = 250, fr = 45;
  await mouse("mousePressed", fx + fr, fy);
  for (let k = 1; k <= 36; k++) await mouse("mouseMoved", fx + fr * Math.cos((k * 2 * Math.PI) / 36), fy + fr * 0.8 * Math.sin((k * 2 * Math.PI) / 36));
  await mouse("mouseReleased", fx + fr, fy);
  await sleep(50);
  await key("Enter");
  obs = await evaluate("FlipSolver.obstacles");
  const fh = obs.find((o) => o.kind === "freehand");
  check("freehand loop closes into a simplified polygon", obs.length === 4 && fh && fh.points.length >= 6 && fh.points.length <= 64, fh && { vertices: fh.points.length });

  // freehand: a straight line is rejected
  await drag(400, 150, 550, 150, 12);
  const rej = await evaluate("document.getElementById('edit-status').textContent");
  await key("Enter");
  check("freehand straight line rejected", (await evaluate("FlipSolver.obstacles.length")) === 4 && /rejected/.test(rej), rej);

  // thickness rule (item 1): a freehand dumbbell whose neck is ~1.2 cells thick is rejected
  // on Enter with a message, stays pending (red cells), and Esc discards it
  {
    const out = [[780, 430], [840, 430], [840, 454], [920, 454], [920, 430], [980, 430], [980, 490], [920, 490], [920, 466], [840, 466], [840, 490], [780, 490], [780, 432]];
    const pts = [];
    for (let k = 0; k + 1 < out.length; k++) {
      const [ax, ay] = out[k], [bx, by] = out[k + 1], n = Math.max(1, Math.round(Math.hypot(bx - ax, by - ay) / 5));
      for (let q = 0; q < n; q++) pts.push([ax + ((bx - ax) * q) / n, ay + ((by - ay) * q) / n]);
    }
    await mouse("mousePressed", pts[0][0], pts[0][1]);
    for (const [x, y] of pts.slice(1)) await mouse("mouseMoved", x, y);
    await mouse("mouseReleased", 780, 432);
    await sleep(50);
    const pendingTxt = await evaluate("document.getElementById('edit-status').textContent");
    await key("Enter");
    const thinTxt = await evaluate("document.getElementById('edit-status').textContent");
    const stillPending = await evaluate("document.getElementById('btn-confirm').disabled === false");
    check("thin freehand neck (1.2 cells) rejected on Enter, stays pending", (await evaluate("FlipSolver.obstacles.length")) === 4 && /Not placed/.test(thinTxt) && /thinner than 2 grid cells/.test(thinTxt) && stillPending, { before: pendingTxt, after: thinTxt });
    await screenshot(path.join(ROOT, "validation", "item1_thin_reject.png"));
    await key("Escape");
    check("Esc discards the rejected thin shape", (await evaluate("FlipSolver.obstacles.length")) === 4 && (await evaluate("document.getElementById('btn-confirm').disabled")) === true);
  }

  // move the circle by dragging it
  const before = obs.find((o) => o.id === circleId);
  await drag(700, 330, 640, 300);
  await key("Enter");
  obs = await evaluate("FlipSolver.obstacles");
  const after = obs.find((o) => o.id === circleId);
  check("drag moves an existing object (same id)", after && Math.abs((after.cx - before.cx) * scale + 60) < 1.5 && Math.abs((after.cy - before.cy) * scale + 30) < 1.5,
    after && { dx_px: +((after.cx - before.cx) * scale).toFixed(1), dy_px: +((after.cy - before.cy) * scale).toFixed(1) });

  // rotate the box: select by click, ] ] then Enter = +30°
  await click(455, 555);
  await key("]");
  await key("]");
  await key("Enter");
  obs = await evaluate("FlipSolver.obstacles");
  let b2 = obs.find((o) => o.type === "rect");
  check("[ ] keys rotate the selected object", b2 && Math.abs((b2.angle * 180) / Math.PI - 30) < 0.01, b2 && { deg: +((b2.angle * 180) / Math.PI).toFixed(2) });
  // wheel rotates, Esc reverts
  await send("Input.dispatchMouseEvent", { type: "mouseWheel", x: 455, y: 555, deltaX: 0, deltaY: 100 });
  await sleep(50);
  await key("Escape");
  obs = await evaluate("FlipSolver.obstacles");
  b2 = obs.find((o) => o.type === "rect");
  check("wheel rotate + Esc leaves the committed box unchanged", b2 && Math.abs((b2.angle * 180) / Math.PI - 30) < 0.01);
  // rotate handle drag (box selected → ] creates pending; handle is above the centre)
  await click(455, 555);
  await key("[");
  await key("]"); // pending copy, angle back to 30°
  const hp = await evaluate("(() => { const o = FlipSolver.obstacles.find(q => q.type === 'rect'); return [o.cx, o.cy, Math.hypot(o.w, o.h) / 2]; })()");
  const ox0 = 0, oy0 = await evaluate("window.innerHeight - Math.floor(window.innerHeight / ((3 / Math.max(64, Math.min(128, Math.round(window.innerWidth / 10)))) * window.innerWidth / 3)) * (3 / Math.max(64, Math.min(128, Math.round(window.innerWidth / 10)))) * window.innerWidth / 3");
  const h = 3 / Math.max(64, Math.min(128, Math.round(W / 10)));
  const toS = (x, y) => [(x - h) * scale + ox0, (y - h) * scale + oy0];
  const [ccx, ccy] = toS(hp[0], hp[1]);
  const rr = hp[2] * scale + 26, a30 = Math.PI / 6;
  const hx = ccx + Math.sin(a30) * rr, hy = ccy - Math.cos(a30) * rr;
  await drag(hx, hy, ccx + rr, ccy, 6); // drag handle to the right → angle 90°
  await key("Enter");
  b2 = (await evaluate("FlipSolver.obstacles")).find((o) => o.type === "rect");
  check("rotate handle drag sets the angle", b2 && Math.abs((b2.angle * 180) / Math.PI - 90) < 1, b2 && { deg: +((b2.angle * 180) / Math.PI).toFixed(1) });

  // delete the freehand object: click it, Delete
  await click(fx, fy);
  await key("Delete");
  check("Delete removes the selected object", (await evaluate("FlipSolver.obstacles.filter(o => o.kind === 'freehand').length")) === 0 && (await evaluate("FlipSolver.obstacles.length")) === 3);

  // touch: triangle (tool 4) via touch drag, then the Confirm button by touch
  await key("4");
  await touchDrag(600, 600, 700, 690);
  if (process.env.SMOKE_DEBUG) console.log("  after touch drag:", await evaluate("document.getElementById('edit-status').textContent + ' | confirm disabled=' + document.getElementById('btn-confirm').disabled"));
  if (process.env.SMOKE_DEBUG) await evaluate("window.__log = []; for (const t of ['pointerdown','pointerup','click','touchstart','touchend']) document.addEventListener(t, (e) => window.__log.push(t + ':' + (e.target.id || e.target.tagName)), true); 0");
  const btn = await evaluate("(() => { const r = document.getElementById('btn-confirm').getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()");
  if (process.env.SMOKE_DEBUG) console.log("  confirm at", btn);
  await send("Input.dispatchTouchEvent", { type: "touchStart", touchPoints: [{ x: btn[0], y: btn[1], id: 2 }] });
  await send("Input.dispatchTouchEvent", { type: "touchEnd", touchPoints: [] });
  await sleep(400); // allow for a synthetic tap -> click delay
  if (process.env.SMOKE_DEBUG) console.log("  events:", await evaluate("window.__log.join(' ')"), "| status:", await evaluate("document.getElementById('edit-status').textContent"));
  obs = await evaluate("FlipSolver.obstacles");
  check("touch drag + Confirm button places a triangle", obs.length === 4 && obs.some((o) => o.kind === "triangle"));

  // leave a pending ghost on screen for the edit-mode screenshot
  await key("1");
  await drag(860, 470, 900, 470);
  await sleep(120);
  await screenshot(path.join(ROOT, "obstacles_editmode.png"));
  await key("Escape");

  // Done button (mouse click on it)
  const done = await evaluate("(() => { const r = document.getElementById('btn-done').getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()");
  await click(done[0], done[1]);
  check("Done exits edit mode and restores the previous (running) state", await evaluate("!document.getElementById('panel').classList.contains('editing') && document.getElementById('btn-pause').textContent === 'Pause'"));
  n0 = await evaluate("FlipSolver.count");
  await mouse("mousePressed", 1100, 120);
  await sleep(60);
  check("after exit: mouse drop works again", await evaluate("FlipSolver.holding === true"));
  await mouse("mouseReleased", 1100, 120);
  await send("Input.dispatchTouchEvent", { type: "touchStart", touchPoints: [{ x: 600, y: 120, id: 3 }] });
  await sleep(60);
  check("after exit: touch drop works", await evaluate("FlipSolver.holding === true"));
  await send("Input.dispatchTouchEvent", { type: "touchEnd", touchPoints: [] });

  // paused-before-edit is restored as paused
  await key("Space");
  await key("e");
  await key("e");
  check("entering/exiting edit mode keeps a user pause", await evaluate("document.getElementById('btn-pause').textContent === 'Resume'"));
  await key("Space");

  await sleep(2500); // let the sim run in real time with 4 obstacles + 2 drops
  const st = await evaluate("FlipSolver.stats()");
  check("live run: no NaN, no particle inside an obstacle", st.nan === 0 && st.insideObstacle === 0 && st.simTime > 1, { simTime: +st.simTime.toFixed(2), count: st.count, insideObstacle: st.insideObstacle, obstacleDeleted: st.obstacleDeleted });
  await key("r");
  await sleep(50);
  const rs = await evaluate("FlipSolver.stats()");
  check("R reset keeps obstacles, spawns nothing inside them", rs.obstacles === 4 && rs.insideObstacle === 0 && rs.simTime < 0.5, { obstacles: rs.obstacles, count: rs.count, simTime: +rs.simTime.toFixed(2) });
  check("no page errors", errors.length === 0, errors.slice(0, 3).join(" | ") || "none");
  }

  // --- preview screenshot: water around a circle and a ramp
  loaded = waitEvent("Page.loadEventFired", 60000);
  await send("Page.navigate", { url: url + "?demoObstacles=1&warm=" + PREVIEW_WARM });
  await loaded;
  await sleep(300);
  await evaluate("document.getElementById('btn-pause').click()"); // freeze the frame
  await sleep(100);
  const ds = await evaluate("FlipSolver.stats()");
  check("demo preload: circle + ramp, no leaks", ds.obstacles === 2 && ds.insideObstacle === 0 && ds.nan === 0, { simTime: +ds.simTime.toFixed(2), count: ds.count, obstacleCells: ds.obstacleCells });
  await screenshot(PREVIEW_OUT);
  if (!PREVIEW_ONLY) await phoneChecks();

  // --- iPhone-sized touch viewport: edit via buttons + touch only
  async function phoneChecks() {
    await send("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: true });
    loaded = waitEvent("Page.loadEventFired", 60000);
    await send("Page.navigate", { url: url + "?demoObstacles=1&warm=1" });
    await loaded;
    await sleep(300);
    const tapEl = async (sel) => {
      const c = await evaluate("(() => { const r = document.querySelector(" + JSON.stringify(sel) + ").getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()");
      await send("Input.dispatchTouchEvent", { type: "touchStart", touchPoints: [{ x: c[0], y: c[1], id: 7 }] });
      await send("Input.dispatchTouchEvent", { type: "touchEnd", touchPoints: [] });
      await sleep(400);
    };
    await tapEl("#btn-edit");
    await tapEl('[data-tool="box"]');
    const panelBottom = await evaluate("document.getElementById('panel').getBoundingClientRect().bottom");
    const y0 = Math.max(panelBottom + 40, 560);
    await touchDrag(120, y0, 260, y0 + 60);
    await tapEl("#btn-confirm");
    const ps = await evaluate("FlipSolver.stats()");
    check("phone: tap Edit, tap Box, touch-drag, tap Confirm", ps.obstacles === 3 && ps.insideObstacle === 0, { panelBottom: Math.round(panelBottom), obstacles: ps.obstacles, removed: ps.lastObstacleDeleted });
    await tapEl('[data-tool="circle"]');
    await touchDrag(300, y0 - 20, 330, y0 - 20);
    await sleep(100);
    await screenshot(path.join(ROOT, "obstacles_editmode_phone.png"));
    await tapEl("#btn-done");
    check("phone: Done commits the pending circle and resumes", (await evaluate("FlipSolver.obstacles.length")) === 4 && (await evaluate("document.getElementById('btn-pause').textContent")) === "Pause");
  }
}

main()
  .catch((e) => { console.error(e); results.push({ name: "exception", ok: false }); })
  .finally(() => {
    try { ws && ws.close(); } catch (e) { /* ignore */ }
    chrome.kill("SIGKILL");
    setTimeout(() => { try { fs.rmSync(profile, { recursive: true, force: true }); } catch (e) { /* ignore */ } }, 300);
    const fail = results.filter((r) => !r.ok).length;
    console.log((fail ? "FAIL" : "PASS") + " — " + (results.length - fail) + "/" + results.length + " browser checks");
    setTimeout(() => process.exit(fail ? 1 : 0), 400);
  });
