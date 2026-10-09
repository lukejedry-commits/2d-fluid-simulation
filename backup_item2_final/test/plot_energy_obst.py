"""Energy plots from the existing *_obst validation data (no simulation).
Usage: python3 test/plot_energy_obst.py
  -> validation/energy_contact_obst.png, validation/thin_barrier_energy_obst.png
Inputs: energy_contact_series_*_obst.csv, energy_contact_obst.json (gate numbers, for cross-check),
        energy_series_obst.json (E(t) per frame, from test/energy_series_obst.js),
        thin_barrier_energy_obst.json (E(t), vmax(t) every 1 s, from test/thin_barrier_energy_obst.js)."""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

V = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "validation")
GATE = 1e-4
def jl(n):
    with open(os.path.join(V, n)) as f: return json.load(f)
def series(n):
    with open(os.path.join(V, n)) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}

gates = jl("energy_contact_obst.json")
ES = jl("energy_series_obst.json")
report = {}

# ------------------------------------------------------------------ 1. energy at contact
cases = [("box", "cfl1", "r2_box", "r2_none", "Box, CFL 1"), ("ramp", "cfl1", "r2_ramp", "r2_none", "Ramp, CFL 1"),
         ("box_cfl05", "cfl05", "r2_box_cfl05", "r2_none_cfl05", "Box, CFL 0.5"), ("ramp_cfl05", "cfl05", "r2_ramp_cfl05", "r2_none_cfl05", "Ramp, CFL 0.5")]
fig = plt.figure(figsize=(12, 11), dpi=100)
gs = GridSpec(4, 2, figure=fig, height_ratios=[1, 2.4, 1, 2.4], hspace=0.12, wspace=0.18,
              left=0.075, right=0.985, top=0.925, bottom=0.095)
for idx, (tag, cfl, key, base, title) in enumerate(cases):
    s = series(f"energy_contact_series_{tag}_obst.csv")
    t, pe, ke, vr = s["t_end"], s["pushPE"], s["pushKE_other"], s["velRemoval"]
    tot = pe + ke + vr
    cpe, cvr, ctot = np.cumsum(pe), np.cumsum(vr), np.cumsum(tot)
    # rolling 0.5 s window, closed: substeps with t_end in [t_j - 0.5, t_j] (1e-7 absorbs CSV time rounding).
    # obstacle_gates.js used (t_j - 0.5, t_j] on its float clock, so whether the substep exactly 0.5 s back is
    # included there depends on float drift; the two definitions differ by < 1 % in the maxima.
    pref = np.concatenate([[0.0], ctot])
    i0 = np.searchsorted(t, t - 0.5 - 1e-7, side="right")
    win = pref[1:] - pref[i0]
    jmax = int(np.argmax(win))
    g = gates[cfl][tag.split("_")[0]]
    assert abs(win[jmax] / g["maxWindow05s"] - 1) < 1e-2 and abs(ctot[-1] - g["pushTotal"]) < 1e-6, tag  # CSV has 7 significant digits
    r, c = 2 * (idx // 2), idx % 2
    axE = fig.add_subplot(gs[r, c])
    ax = fig.add_subplot(gs[r + 1, c], sharex=axE)
    e, b = np.array(ES[key]["rows"]), np.array(ES[base]["rows"])
    axE.plot(e[:, 0], e[:, 1], color="C0", lw=1.4, label=f"E(t) {tag.split('_')[0]}")
    axE.plot(b[:, 0], b[:, 1], color="0.45", lw=1.2, ls="--", label="E(t) no obstacle")
    axE.set_ylabel("E [J/kg]"); axE.set_ylim(0, 5.2); axE.grid(alpha=0.3)
    axE.set_title(f"{title} — max 0.5 s window {win[jmax]:+.2e} J/kg ({t[i0[jmax]]:.2f}–{t[jmax]:.2f} s) → FAIL", fontsize=10)
    axE.legend(loc="upper right", fontsize=8.5, ncol=2)
    plt.setp(axE.get_xticklabels(), visible=False)
    ax.plot(t, cpe, color="C3", lw=1.6, label="cumulative push-out PE")
    ax.plot(t, cvr, color="C2", lw=1.6, label="cumulative normal-velocity removal")
    ax.plot(t, ctot, color="k", lw=1.8, label="cumulative push-out (PE+KE) + removal")
    ax.plot(t, win, color="C1", lw=1.0, alpha=0.9, label="rolling 0.5 s window of push-out + removal")
    ax.axhline(GATE, color="m", ls="--", lw=1.4, label="gate +1e-4 J/kg")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.plot([t[jmax]], [win[jmax]], "o", ms=6, mfc="none", mec="C1", mew=1.6)
    ax.annotate(f"max window {win[jmax]:+.2e}\n@ {t[i0[jmax]]:.2f}–{t[jmax]:.2f} s", (t[jmax], win[jmax]),
                xytext=(12, -1.5e-4), textcoords="data", fontsize=8.5, color="C1",
                arrowprops=dict(arrowstyle="->", color="C1", lw=0.9))
    ax.set_yscale("symlog", linthresh=1e-5, linscale=0.6)
    ax.set_ylim(-0.06, 0.06)
    ax.set_yticks([-1e-2, -1e-3, -1e-4, -1e-5, 0, 1e-5, 1e-4, 1e-3, 1e-2])
    ax.grid(alpha=0.3, which="major")
    ax.set_xlim(0, 30)
    ax.set_ylabel("J/kg  (symlog, linear |y|<1e-5)")
    if r + 1 == 3: ax.set_xlabel("t [s]")
    else: plt.setp(ax.get_xticklabels(), visible=False)
    ax.text(0.99, 0.03, f"end: PE {cpe[-1]:+.4f}, removal {cvr[-1]:+.4f}, sum {ctot[-1]:+.4f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.5, bbox=dict(fc="w", ec="0.7", alpha=0.85))
    report[tag] = dict(maxWin=win[jmax], at=(t[i0[jmax]], t[jmax]), cum=ctot[-1], cpe=cpe[-1], cvr=cvr[-1],
                       firstOverGate=float(t[np.argmax(win > GATE)]) if (win > GATE).any() else None,
                       timeOverGate=float(np.sum(np.diff(np.concatenate([[0], t]))[win > GATE])))
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Run 2 — energy at obstacle contact (dam break, damping 0, 30 s): push-out / removal budget rows per substep, cumulated\n"
             "small top panels: total E(t) = KE + g·height per unit mass, case vs no-obstacle baseline at the same CFL (sampled every frame)", fontsize=11)
fig.savefig(os.path.join(V, "energy_contact_obst.png"))
plt.close(fig)

# ------------------------------------------------------------------ 2. thin barriers
T = jl("thin_barrier_energy_obst.json")
by = {s["case"]: np.array([[r["t"], r["e"], r["vmax"]] for r in s["rows"]]) for s in T["series"]}
fr = np.array(ES["r3_b1_cfl05"]["rows"])  # per frame: t, E, KE, vmax, caps
one_s = []
k = 0
for row in fr:
    if row[0] >= k - 1e-9:
        one_s.append([row[0], row[1], row[3]]); k += 1
by["b1_cfl05"] = np.array(one_s)
spec = [("b1", "1 cell [64,65], CFL 1", "C3", "-"), ("b1_cfl05", "1 cell [64,65], CFL 0.5", "C1", "-"),
        ("b05c", "0.5 cell over cell centre [64.25,64.75], CFL 1", "C4", "-"),
        ("b2", "2 cells [64,66], CFL 1", "C0", "-"), ("b05e", "0.5 cell, no cell centre (unmarked), CFL 1", "C2", "-"),
        ("none", "no barrier (baseline), CFL 1", "k", "--")]
fig, (ax, axv) = plt.subplots(2, 1, figsize=(12, 9), dpi=100, sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1], hspace=0.08))
CL = 14.0
for key, lab, col, ls in spec:
    d = by[key]
    ax.plot(d[:, 0], d[:, 1], color=col, ls=ls, lw=1.8, marker="o", ms=3.5, label=lab)
    axv.plot(d[:, 0], d[:, 2], color=col, ls=ls, lw=1.4, marker="o", ms=3, label=lab)
    clamp = d[:, 2] >= CL - 0.01
    if clamp.any():
        t0 = d[np.argmax(clamp), 0]
        plateau = d[(d[:, 0] >= 10) & clamp, 1]
        if len(plateau):
            ax.hlines(np.median(plateau), t0, 30, colors=col, linestyles="dashed", lw=1.1, alpha=0.9)
            ax.text(30.3, np.median(plateau) + {"b1": 1.2, "b05c": -1.2}.get(key, 0), f"{np.median(plateau):.1f}", color=col, va="center", fontsize=9)
        report["clamp_" + key] = dict(first_sample_at_clamp=t0, plateau_median=float(np.median(plateau)), Emax=float(d[:, 1].max()), tEmax=float(d[np.argmax(d[:, 1]), 0]))
ax.plot(fr[:, 0], fr[:, 1], color="C1", lw=0.5, alpha=0.35, zorder=0, label="1 cell CFL 0.5, per frame (faint)")
onset = fr[np.argmax(fr[:, 3] >= CL - 0.01), 0]
e2x = fr[np.argmax(fr[:, 1] >= 2 * fr[0, 1]), 0]
report["b1_cfl05_frame"] = dict(first_frame_at_clamp=float(onset), E_doubles_at=float(e2x))
ax.axvspan(1, 30, color="0.85", alpha=0.35, zorder=-1)
ax.text(1.3, 58, "grey band: 1-cell & centred 0.5-cell cases have particles at the 14 m/s speed clamp in every sample (t ≥ 1 s)\n"
        "dashed horizontals: clamp-limited plateau (median E over clamp samples with t ≥ 10 s)", fontsize=8.8, va="top")
ax.set_ylabel("total E = KE + g·height  [J/kg]"); ax.set_ylim(0, 60); ax.grid(alpha=0.3)
h_, l_ = ax.get_legend_handles_labels()
axv.legend(h_, l_, loc="center right", fontsize=8.6, ncol=2, framealpha=0.95)
axv.axhline(CL, color="r", ls="--", lw=1.3)
axv.text(0.3, CL + 0.4, "speed clamp 14 m/s", color="r", fontsize=8.8)
axv.set_ylabel("max particle speed [m/s]"); axv.set_ylim(0, 17); axv.grid(alpha=0.3)
axv.set_xlabel("t [s]"); axv.set_xlim(0, 31.2)
ax.set_title("Run 3 — full-height barriers, dam break, damping 0: total energy vs time\n"
             "sampled every 1 s of sim time (1-cell CFL 0.5 also shown per frame, faint)", fontsize=11)
fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.07)
fig.savefig(os.path.join(V, "thin_barrier_energy_obst.png"))
plt.close(fig)
print(json.dumps(report, indent=1, default=float))
