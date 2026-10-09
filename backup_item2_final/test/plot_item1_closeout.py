"""Item 1 close-out: plots + analysis (no simulation; reads CSV/JSON written by the harnesses).
Usage: python3 test/plot_item1_closeout.py [summary.json from test/item1_closeout.js]
Writes validation/item1_energy.png, item1_stage_rows.png, item1_api_tilt.png,
       validation/item1_api_tilt.json, validation/item1_dt_check.json
"""
import csv, json, math, sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
V = os.path.join(ROOT, "validation")
H = 3 / 128; G = 9.81; RHO = 1000.0
T_CONTACT = 0.16666666666666669

def rcsv(p):
    with open(p) as f: return list(csv.DictReader(f))

# ---------------------------------------------------------------- series (Phase A/B budget runs)
def series(path):
    rows = list(csv.reader(open(path)))
    hdr = rows[0]; out = {}
    for r in rows[1:]:
        out.setdefault(r[0], []).append(r[1:])
    return hdr[1:], out
hdrA, serPre = series(os.path.join(V, "item1_stage_series.csv"))
hdrB, serFix = series(os.path.join(V, "item1_stage_series_fixed.csv"))
budPre = json.load(open(os.path.join(V, "item1_stage_budget.json")))
budFix = json.load(open(os.path.join(V, "item1_stage_budget_fixed.json")))
E0 = budPre["cases"]["b1_cfl1"]["E0"]

def num(x):
    try: return float(x)
    except: return float("nan")

def arr(ser, case, hdr):
    rows = ser[case]
    t0 = np.array([num(r[0]) for r in rows]); dt = np.array([num(r[1]) for r in rows]); dE = np.array([num(r[2]) for r in rows])
    stages = {h[2:]: np.array([num(r[k]) for r in rows]) for k, h in enumerate(hdr) if h.startswith("A:")}
    return t0, dt, dE, stages

# ---------------------------------------------------------------- 1a energy plot
gE = rcsv(os.path.join(V, "item1_gates_energy.csv"))
tg = np.array([float(r["t"]) for r in gE]); Eb1 = np.array([float(r["E_b1"]) for r in gE]); Eb2 = np.array([float(r["E_b2"]) for r in gE])
bE = rcsv(os.path.join(V, "item1_baseline_energy.csv"))
tb = np.array([float(r["t"]) for r in bE]); Eb0 = np.array([float(r["E_none"]) for r in bE])
t0p, dtp, dEp, _ = arr(serPre, "b1_cfl1", hdrA)
ok = np.isfinite(dEp); last = np.argmax(~ok) if (~ok).any() else len(dEp)
tp = np.concatenate([[0], (t0p + dtp)[:last]]); Ep = E0 + np.concatenate([[0], np.cumsum(dEp[:last])])
t_nan = budPre["cases"]["b1_cfl1"]["simTime"]
# consistency: fixed-b1 series cumulative vs gates frames (first 1 s)
t0f, dtf, dEf, stF = arr(serFix, "b1_cfl1", hdrB)
tf = t0f + dtf; Ef = E0 + np.cumsum(dEf)
cons = float(np.max(np.abs(np.interp(tg[(tg > 0) & (tg <= 1)], tf, Ef) - Eb1[(tg > 0) & (tg <= 1)])))

plt.rcParams.update({"font.size": 12})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 9.6), dpi=100, gridspec_kw={"height_ratios": [1, 1.15]})
for ax in (a1,):
    ax.semilogy(tp, Ep / E0, color="#d62728", lw=2.2, label="b1 pre-fix (backup_pre_item1/sim.js), per substep")
    ax.semilogy(tg, Eb1 / E0, color="#1f77b4", lw=2, label="b1 fixed (1-cell barrier)")
    ax.semilogy(tg, Eb2 / E0, color="#2ca02c", lw=1.6, ls="--", label="b2 (2-cell barrier)")
    ax.semilogy(tb, Eb0 / E0, color="#7f7f7f", lw=1.6, label="no obstacle")
    ax.axvline(T_CONTACT, color="k", ls=":", lw=1); ax.text(T_CONTACT, 1e40, " contact\n 0.167 s", fontsize=11)
    ax.axvline(t_nan, color="#d62728", ls=":", lw=1.2); ax.text(t_nan, 1e10, " NaN\n %.3f s" % t_nan, color="#d62728", fontsize=11)
    ax.set_xlim(0, 0.5); ax.set_ylim(0.3, 1e80); ax.set_ylabel("E / E0  (log)"); ax.set_xlabel("t (s)")
    ax.set_title("Dam break, damping 0, speed clamp off, CFL 1: first 0.5 s (log y)")
    ax.legend(loc="center right", fontsize=10.5); ax.grid(alpha=0.3, which="major")
a2.plot(tg, Eb1 / E0, color="#1f77b4", lw=1.8, label="b1 fixed: E(30)/E0 = %.4f" % (Eb1[-1] / E0))
a2.plot(tg, Eb2 / E0, color="#2ca02c", lw=1.4, ls="--", label="b2: E(30)/E0 = %.4f" % (Eb2[-1] / E0))
a2.plot(tb, Eb0 / E0, color="#7f7f7f", lw=1.4, label="no obstacle: E(30)/E0 = %.4f" % (Eb0[-1] / E0))
m = Ep / E0 < 1.12
a2.plot(tp[m], Ep[m] / E0, color="#d62728", lw=2.2, label="b1 pre-fix (leaves the axis at %.3f s, NaN at %.3f s)" % (tp[np.argmax(~m)] if (~m).any() else tp[-1], t_nan))
a2.annotate("pre-fix b1 → ∞", xy=(0.3, 1.1), xytext=(2.5, 1.07), color="#d62728", arrowprops=dict(arrowstyle="->", color="#d62728"), fontsize=11)
a2.axhline(1, color="k", lw=0.8)
a2.set_xlim(0, 30); a2.set_ylim(0.35, 1.12); a2.set_xlabel("t (s)"); a2.set_ylabel("E / E0  (linear)")
a2.set_title("Full 30 s (linear y); E0 = %.4f J/kg (same initial state for all four)" % E0)
a2.legend(loc="lower left", fontsize=10.5); a2.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(V, "item1_energy.png")); plt.close(fig)

# ---------------------------------------------------------------- 1b stage rows
HL = {"extrapolation": ("#d62728", "extrapolation"), "boundary: obstacle ghost-face fill": ("#ff7f0e", "obstacle ghost-face fill")}
def cum_rows(ser, hdr, case):
    t0s, dts, dEs, st = arr(ser, case, hdr)
    te = t0s + dts; good = np.isfinite(dEs)
    for v in st.values(): good &= np.isfinite(v)
    n = np.argmax(~good) if (~good).any() else len(te)
    return te[:n], np.cumsum(dEs[:n]), {k: np.cumsum(v[:n]) for k, v in st.items()}
pre = cum_rows(serPre, hdrA, "b1_cfl1"); fix = cum_rows(serFix, hdrB, "b1_cfl1")
fig = plt.figure(figsize=(12, 10), dpi=100)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1], sharey=axB)
def draw(ax, data, title, tnan, xlim, lin, legend=False):
    te, cE, rows = data
    for k, c in rows.items():
        if k in HL or np.max(np.abs(c)) < 1e-12: continue
        ax.plot(te, c, color="0.62", lw=0.9)
    for k, (col, lab) in HL.items(): ax.plot(te, rows[k], color=col, lw=2.6, label=lab)
    ax.plot(te, cE, color="k", lw=1.6, ls="--", label="measured ΔE (= sum of all rows)")
    ax.plot([], [], color="0.62", lw=0.9, label="other rows (P2G, projection, walls, viscosity, ...)")
    ax.axvline(T_CONTACT, color="k", ls=":", lw=1.1)
    ax.set_yscale("symlog", linthresh=lin); ax.set_xlim(*xlim); ax.set_xlabel("t (s)"); ax.set_title(title, fontsize=12); ax.grid(alpha=0.3)
    yl = ax.get_ylim(); ax.text(T_CONTACT, yl[1] * 0.3 if yl[1] > 0 else 0, " contact %.3f s" % T_CONTACT, fontsize=10, va="top")
    if tnan: ax.axvline(tnan, color="#d62728", ls=":", lw=1.3); ax.text(tnan, yl[0] * 0.3, " NaN %.4f s" % tnan, color="#d62728", fontsize=10)
    if legend: ax.legend(loc="upper left", fontsize=10)
draw(axA, pre, "b1 PRE-FIX, full scale (symlog, linear for |y| < 0.01): rows blow up to ~1e71 J/kg before NaN", t_nan, (0, 0.25), 1e-2, legend=True)
axA.set_ylabel("cumulative ΔE by stage (J/kg)")
for ax, d, ttl, tn in [(axB, pre, "b1 PRE-FIX, first 1 s (y clipped to ±30)", t_nan), (axC, fix, "b1 FIXED, first 1 s (same y scale)", None)]:
    draw(ax, d, ttl, None, (0, 1), 1e-2)
    ax.set_ylim(-30, 30)
    for tx in ax.texts[:]: tx.remove()
    ax.text(T_CONTACT, 12, " contact", fontsize=10)
    if tn: ax.axvline(tn, color="#d62728", ls=":", lw=1.3); ax.text(tn, -12, " NaN %.4f s" % tn, color="#d62728", fontsize=10)
axB.set_ylabel("cumulative ΔE by stage (J/kg), symlog")
fig.suptitle("Cumulative per-stage energy rows (method A), CFL 1, damping 0, speed clamp off", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(V, "item1_stage_rows.png")); plt.close(fig)

# ---------------------------------------------------------------- 3 dt check
def frames_caps(t0, dt, target=1 / 60, maxsub=8, tol=1e-9):
    """Reconstruct frames: a frame ends when its substeps cover 1/60 s; 8 substeps short of it = cap hit."""
    cap = np.zeros(len(dt), bool); capframes = []; acc = 0; k0 = 0; n = 0
    for k in range(len(dt)):
        acc += dt[k]; n += 1
        if abs(acc - target) < tol or acc > target:
            acc = 0; n = 0; k0 = k + 1
        elif n == maxsub:
            cap[k0:k + 1] = True; capframes.append((k0, k)); acc = 0; n = 0; k0 = k + 1
    return cap, capframes
DTX = {}
for r in rcsv(os.path.join(V, "item1_dt_series.csv")):
    DTX.setdefault((r["tag"], r["case"]), ([], []))
    DTX[(r["tag"], r["case"])][0].append(float(r["t0"])); DTX[(r["tag"], r["case"])][1].append(float(r["dt"]))
DTX = {k: (np.array(v[0]), np.array(v[1])) for k, v in DTX.items()}
dtc = {"rule": "p and halving claims stand only if mean dt(CFL1)/mean dt(CFL0.5) is within 10% of 2 (1.8-2.2)",
       "windowMembership": "substep start time t0 in [a, b), exactly as test/item1_budget.js; full-precision dt from validation/item1_dt_series.csv", "capHitInference": "frames reconstructed from the per-substep dt series (8 substeps not covering 1/60 s = cap hit); totals cross-checked against the harness capLog count", "cases": {}}
for tag, ser, hdr, bud, wins in [("prefix", serPre, hdrA, budPre, ["onsetWindow", "growthWindow"]), ("fixed", serFix, hdrB, budFix, ["first1s", "contactWindow"])]:
    for case in ["b1_cfl1", "b1_cfl05"]:
        t0s, dts = DTX[(tag, case)]
        cap, cf = frames_caps(t0s, dts, tol=1e-9)
        dtc["cases"].setdefault(tag, {})[case + "_capFramesTotal_inferred"] = len(cf)
        dtc["cases"][tag][case + "_capFramesTotal_harness"] = bud["cases"][case]["capHits"]
    for w in wins:
        res = {}
        for case in ["b1_cfl1", "b1_cfl05"]:
            t0s, dts, dEs, _ = arr(ser, case, hdr)
            cap, cf = frames_caps(t0s, dts)
            a, b = bud["cases"][case][w]["window"]
            t0s, dts = DTX[(tag, case)]  # full-precision t0, dt (validation/item1_dt_series.csv)
            cap, cf = frames_caps(t0s, dts)
            m = (t0s >= a - 1e-12) & (t0s < b - 1e-12)  # same membership as test/item1_budget.js
            capF = sum(1 for (k0, k1) in cf if np.any(m[k0:k1 + 1]))
            res[case] = {"window": [a, b], "substeps": int(m.sum()), "substeps_budget": bud["cases"][case][w]["substeps"], "mean_dt": float(dts[m].mean()), "min_dt": float(dts[m].min()), "max_dt": float(dts[m].max()),
                         "capHitFramesInWindow": capF, "substepsInCapFrames": int((cap & m).sum())}
        ratio = res["b1_cfl1"]["mean_dt"] / res["b1_cfl05"]["mean_dt"]
        res["dt_ratio_cfl1_over_cfl05"] = ratio; res["p_stands"] = bool(1.8 <= ratio <= 2.2)
        dtc["cases"][tag][w] = res
json.dump(dtc, open(os.path.join(V, "item1_dt_check.json"), "w"), indent=1)

# ---------------------------------------------------------------- 2 API tilt
summ = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "/tmp/item1co_summary.json"))
def fser(suf):
    r = rcsv(os.path.join(V, "item1_api_tilt_force%s.csv" % suf))
    return {k: np.array([float(x[k]) for x in r]) for k in r[0]}
def eser(suf):
    r = rcsv(os.path.join(V, "item1_api_tilt_energy%s.csv" % suf))
    return np.array([float(x["t"]) for x in r]), np.array([float(x["E"]) for x in r])
def period_fft(t, y):
    tt = np.arange(t[0], t[-1], 1 / 120); yy = np.interp(tt, t, y); yy = yy - yy.mean()
    yy = yy * np.hanning(len(yy)); N = 1 << 20
    Y = np.abs(np.fft.rfft(yy, N)); f = np.fft.rfftfreq(N, tt[1] - tt[0])
    band = (f > 0.2) & (f < 3); k = np.argmax(np.where(band, Y, 0))
    return 1 / f[k]
def period_zc(t, y):
    # upward zero crossings of the signal low-passed with a 0.4 s moving average (removes the
    # high-frequency content that would add spurious crossings)
    nw = max(1, int(round(0.4 / np.median(np.diff(t))))); yy = np.convolve(y - y.mean(), np.ones(nw) / nw, mode="same")
    s = np.sign(yy); idx = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    tz = t[idx] - yy[idx] * (t[idx + 1] - t[idx]) / (yy[idx + 1] - yy[idx])
    return float(np.mean(np.diff(tz))) if len(tz) > 2 else None, len(tz)
api = {"setup": "see test/item1_closeout.js header", "cases": {}}
for name, suf in [("tilt_b1", ""), ("tilt_b1_eqN", "_eqN")]:
    S = summ[name]; F = fser(suf); t = F["t"]
    dL, dR = S["d_columnDepth_left_m"], S["d_columnDepth_right_m"]
    sd = S["seed"]; LL, LR = sd[0]["width_m"], sd[1]["width_m"]
    def Tth(L, d): k = math.pi / L; return 2 * math.pi / math.sqrt(G * k * math.tanh(k * d))
    Tm = period_fft(t, F["F_net_N_per_m"]); Tz, nz = period_zc(t, F["F_net_N_per_m"])
    TmL = period_fft(t, F["F_left_N_per_m"]); TmR = period_fft(t, F["F_right_N_per_m"])
    nper = int(math.floor((t[-1] + 1e-9) / Tm)); w = t <= nper * Tm + 1e-9
    dbar = 0.5 * (dL + dR)
    api["cases"][name] = {
        "count_left": S["count_left"], "count_right": S["count_right"], "count_diff_pct": S["count_diff_pct"], "counts_within_1pct": abs(S["count_diff_pct"]) <= 1.0,
        "d_left_m_columnDepth": dL, "d_right_m_columnDepth": dR, "d_left_m_fromCount": sd[0]["d_fromCount_m"], "d_right_m_fromCount": sd[1]["d_fromCount_m"],
        "width_left_m": LL, "width_right_m": LR, "tilt_amplitude_cells_at_ends": [sd[0]["amplitude_cells_at_ends"], sd[1]["amplitude_cells_at_ends"]],
        "T_theory_left_s": Tth(LL, dL), "T_theory_right_s": Tth(LR, dR),
        "T_measured_Fnet_fft_s": Tm, "T_measured_Fnet_zerocross_lowpass_s": Tz, "zero_crossings": nz, "T_measured_Fleft_fft_s": TmL, "T_measured_Fright_fft_s": TmR,
        "mean_window_s": [0, nper * Tm], "periods_in_window": nper,
        "mean_Fnet_N_per_m": float(np.mean(F["F_net_N_per_m"][w])), "mean_Fnet_full30s": float(np.mean(F["F_net_N_per_m"])),
        "bound_rho_g_dbar_h": RHO * G * dbar * H, "mean_Fnet_within_bound": bool(abs(np.mean(F["F_net_N_per_m"][w])) <= RHO * G * dbar * H),
        "Fnet_initial": float(F["F_net_N_per_m"][0]), "Fnet_rms_after_mean_removed": float(np.std(F["F_net_N_per_m"][w])),
        "mean_Fleft_N_per_m": float(np.mean(F["F_left_N_per_m"][w])), "mean_Fright_N_per_m": float(np.mean(F["F_right_N_per_m"][w])),
        "hydrostatic_left_half_rho_g_d2": 0.5 * RHO * G * dL ** 2, "hydrostatic_right_half_rho_g_d2": 0.5 * RHO * G * dR ** 2,
        "E0": S["E0"], "E30": S["E30"], "E30_le_E0": S["E30_le_E0"], "capHits": S["capHits"], "pcgCapHits": S["pcgCapHits"], "pcgSolves": S["pcgSolves"], "pcgMaxIters": S["pcgMaxIters"],
        "maxSpeed": S["maxParticleSpeed"], "leaks_everCrossed": S["leaks_everCrossed"], "minSdfOverR": S["minSdfOverR"], "insideMax": S["insideMaxPerSubstep"], "deleted": S["obstacleDeletedEnd"], "nanGuard": S["nanGuard"],
        "thinGhost": S["thinGhost"]}
nb = summ["tilt_none"]
api["noBarrier"] = {k: nb[k] for k in ["E0", "E30", "E30_le_E0", "capHits", "pcgCapHits", "maxParticleSpeed", "count_left", "count_right", "nanGuard"]}
for name in ["tilt_b1", "tilt_b1_eqN"]:
    c = api["cases"][name]; c["maxSpeed_ratio_vs_noBarrier"] = c["maxSpeed"] / nb["maxParticleSpeed"]
    c["expectations"] = {"E30_le_E0": c["E30_le_E0"], "capHits_0": c["capHits"] == 0, "maxSpeed_le_2x_noBarrier": c["maxSpeed_ratio_vs_noBarrier"] <= 2}
json.dump(api, open(os.path.join(V, "item1_api_tilt.json"), "w"), indent=1)

F = fser(""); Fe = fser("_eqN"); te1, E1 = eser(""); ten, En = eser("_nobarrier"); tee, Ee = eser("_eqN")
c = api["cases"]["tilt_b1"]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 9), dpi=100, sharex=True)
a1.plot(te1, E1, color="#1f77b4", lw=1.8, label="1-cell barrier (API): E(30)/E0 = %.4f" % (E1[-1] / E1[0]))
a1.plot(ten, En, color="#7f7f7f", lw=1.5, label="no barrier, same particles: E(30)/E0 = %.4f" % (En[-1] / En[0]))
a1.plot(tee, Ee, color="#9467bd", lw=1.0, ls="--", label="extra: barrier, equal counts (right d lower): E(30)/E0 = %.4f" % (Ee[-1] / Ee[0]))
a1.set_ylabel("E (J/kg)"); a1.legend(fontsize=10.5, loc="upper right"); a1.grid(alpha=0.3)
a1.set_title("Tilted pools on both sides of a 1-cell solid (addObstacle), damping 0, clamp off, CFL 1")
a2.plot(F["t"], F["F_net_N_per_m"], color="#1f77b4", lw=1.3, label="net force F_left - F_right (equal depth)")
a2.plot(Fe["t"], Fe["F_net_N_per_m"], color="#9467bd", lw=0.9, ls="--", label="extra: equal counts")
a2.axhline(c["mean_Fnet_N_per_m"], color="#d62728", lw=1.6, label="mean over %d periods (%.2f s): %.1f N/m" % (c["periods_in_window"], c["mean_window_s"][1], c["mean_Fnet_N_per_m"]))
for s in (1, -1): a2.axhline(s * c["bound_rho_g_dbar_h"], color="k", lw=0.8, ls=":")
a2.text(0.3, c["bound_rho_g_dbar_h"] + 8, "±ρ g d h = ±%.0f N/m" % c["bound_rho_g_dbar_h"], fontsize=10)
a2.set_xlabel("t (s)"); a2.set_ylabel("net horizontal force on barrier (N/m)")
a2.set_title("T measured %.3f s (FFT), theory %.3f / %.3f s (left / right, ω² = g k tanh kd, k = π/L)" % (c["T_measured_Fnet_fft_s"], c["T_theory_left_s"], c["T_theory_right_s"]), fontsize=11.5)
a2.legend(fontsize=10, loc="upper right"); a2.grid(alpha=0.3); a2.set_xlim(0, 30)
fig.tight_layout(); fig.savefig(os.path.join(V, "item1_api_tilt.png")); plt.close(fig)
print(json.dumps({"consistency_fixed_series_vs_gates_maxAbsDiff_first1s": cons}, indent=0))
print(json.dumps(dtc, indent=1)); print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "thinGhost"} for k, v in api["cases"].items()}, indent=1)); print(api["noBarrier"])
