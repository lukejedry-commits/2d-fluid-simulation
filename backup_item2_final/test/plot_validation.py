"""Plots for validation/ (needs matplotlib).
Usage: python3 test/plot_validation.py            -> front_ZT.png, energy.png from the round-2 GS files (front_h2.csv ...)
       python3 test/plot_validation.py _pcg       -> front_ZT_pcg.png (release + shifted clock), energy_pcg.png
       python3 test/plot_validation.py _2c        -> front_ZT_2c.png, energy_2c.png (any suffix with a summary<suffix>.json)
"""
import csv, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = os.path.join(os.path.dirname(__file__), "..", "validation")
def rd(name):
    with open(os.path.join(V, name)) as f:
        return list(csv.DictReader(f))

SUF = sys.argv[1] if len(sys.argv) > 1 else ""
ref = rd("mm_reference.csv")
shifts = {}
if SUF:
    with open(os.path.join(V, f"summary{SUF}.json")) as f:
        summ = json.load(f)
    shifts = {k: (summ["comparison"][k]["shift"], summ["comparison"][k]["T144"]) for k in ("h2", "h1")}
clocks = ["release", "shifted"] if SUF else ["release"]
fig, axes = plt.subplots(len(clocks), 2, figsize=(12, 5 * len(clocks)), squeeze=False)
panels = []
for row, clock in enumerate(clocks):
    for col, (tag, n2, label) in enumerate([("h2", 2, "H/a = 2 (24x48 cells)"), ("h1", 1, "H/a = 1 (24x24 cells)")]):
        panels.append((axes[row][col], tag, n2, label, clock))
for ax, tag, n2, label, clock in panels:
    fr = rd(f"front_{tag}{SUF}.csv")
    dT = shifts[tag][0] if clock == "shifted" else 0.0
    if clock == "shifted":
        label += f" — shifted clock: T_1.44 = {shifts[tag][1]:.3f}, shift = {dT:+.3f}"
    elif SUF:
        label += " — release clock"
    ax.plot([float(r["T"]) + dT for r in fr], [float(r["Z"]) for r in fr], "-", color="tab:blue", lw=2,
            label=f"FLIP/APIC sim ({ {'_pcg': 'PCG+MIC(0)', '_2c': 'PCG, neutral separation', '_2c_legacy': 'PCG, legacy separation', '_2d': 'PCG, neutral sep., mirror walls', '_2d_base2c': '2d off (2c state)', '_2d_baseR3': 'all legacy (round 3)', '_2g': 'PCG, separating walls', '_2g_base2d': '2g off (2d mirror)', '_2g_baseR3': 'all legacy'}.get(SUF, 'GS' if SUF in ('', '_gs') else SUF) }, damping 0)")
    for a_in, mk in sorted({(r["a_in"], None) for r in ref if int(r["n2"]) == n2}):
        pts = [r for r in ref if int(r["n2"]) == n2 and r["a_in"] == a_in and float(r["Z"]) <= 6.0]
        ax.plot([float(r["T"]) for r in pts], [float(r["Z"]) for r in pts], "o" if a_in in ("1.125", "2.25") and n2 == 2 and a_in == "1.125" else "s",
                mfc="none", label=f"Martin & Moyce 1952, n²={n2}, a={a_in} in. (table mean)")
    ax.axhline(5.33, color="gray", ls=":", lw=1)
    ax.text(0.1, 5.38, "far wall Z = 5.33", color="gray", fontsize=8)
    ax.axvline(0.5, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("T = t·sqrt(gH)/a")
    ax.set_ylabel("Z = x_front / a")
    ax.set_title(label, fontsize=9)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(0.8, 5.8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
fig.suptitle("Surge front: sim vs Martin & Moyce (1952) Tables 1–2 (reference T normalized by the authors at Z = 1.44)"
             + ("; dashed line = comparison window start T = 0.5" if SUF else ""), fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(V, f"front_ZT{SUF}.png"), dpi=130)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, tag, target, title in [(axes[0], "h2", 1.034, "H/a = 2"), (axes[1], "h1", None, "H/a = 1")]:
    e = rd(f"energy_{tag}{SUF}.csv")
    t = [float(r["t_s"]) for r in e]
    ax.plot(t, [float(r["KE_J_per_kg"]) for r in e], label="KE")
    ax.plot(t, [float(r["PE_J_per_kg"]) for r in e], label="PE (from floor)")
    ax.plot(t, [float(r["total_J_per_kg"]) for r in e], "k", lw=1.5, label="total")
    if target:
        ax.axhline(target, color="gray", ls=":", lw=1)
        ax.text(20, target + 0.12, "g·0.2109/2 ≈ 1.03", color="gray", fontsize=8)
        ax.axhline(5.52, color="gray", ls=":", lw=1)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("J/kg")
    ax.set_title(f"Energy per kg, {title}, damping 0" + (f" ({SUF[1:].upper()})" if SUF else ""))
    ax.grid(alpha=0.3)
    ax.legend()
    ins = ax.inset_axes([0.45, 0.45, 0.5, 0.35])
    ins.plot(t, [float(r["total_J_per_kg"]) for r in e], "k", lw=1)
    ins.set_xlim(0, 1.0)
    e0 = float(e[0]["total_J_per_kg"])
    ins.set_ylim(e0 * 0.97, e0 * 1.012)
    ins.axhline(e0, color="gray", ls=":", lw=0.8)
    ins.set_title("total, first 1 s", fontsize=8)
    ins.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(os.path.join(V, f"energy{SUF}.png"), dpi=130)
print("ok")
