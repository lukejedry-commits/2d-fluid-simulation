// Unsmoothed energy-rise metrics shared by test/validate.js and diagnostics.
// rows: arrays t[] and E[] sampled once per frame (60 Hz). E0 = E[0].

/** Largest rise above the running minimum inside [t0, t1]. */
function riseIn(t, E, t0, t1) {
  let mn = Infinity, mnT = 0, best = 0, bT = null, fromT = null;
  for (let i = 0; i < E.length; i++) {
    if (t[i] < t0 || t[i] > t1) continue;
    if (E[i] < mn) { mn = E[i]; mnT = t[i]; }
    if (E[i] - mn > best) { best = E[i] - mn; bT = t[i]; fromT = mnT; }
  }
  return { rise: best, risePctE0: (100 * best) / E[0], fromT, toT: bT };
}

/**
 * Individual upswings ("spikes") in [t0, t1]: zig-zag over the frame series with a
 * reversal threshold of revPct % of E0, so frame-level jitter does not split one
 * spike. Returns every upswing >= minPct % of E0, largest first.
 */
function spikes(t, E, t0, t1, revPct = 0.02, minPct = 0.05) {
  const E0 = E[0], rev = (revPct / 100) * E0;
  const out = [];
  let lo = null, loT = 0, hi = null, hiT = 0, rising = false;
  for (let i = 0; i < E.length; i++) {
    if (t[i] < t0 || t[i] > t1) continue;
    const e = E[i];
    if (lo === null) { lo = hi = e; loT = hiT = t[i]; continue; }
    if (!rising) {
      if (e < lo) { lo = e; loT = t[i]; }
      else if (e - lo > rev) { rising = true; hi = e; hiT = t[i]; }
    } else {
      if (e > hi) { hi = e; hiT = t[i]; }
      else if (hi - e > rev) { out.push({ fromT: loT, toT: hiT, rise: hi - lo }); rising = false; lo = e; loT = t[i]; }
    }
  }
  if (rising) out.push({ fromT: loT, toT: hiT, rise: hi - lo });
  return out
    .map((s) => ({ fromT: +s.fromT.toFixed(3), toT: +s.toT.toFixed(3), rise: s.rise, pctE0: +((100 * s.rise) / E0).toFixed(3) }))
    .filter((s) => s.pctE0 >= minPct)
    .sort((a, b) => b.rise - a.rise);
}

module.exports = { riseIn, spikes };
