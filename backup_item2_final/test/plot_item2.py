#!/usr/bin/env python3
"""Item 2 aggregation: Run 2 gates from /tmp/item2_<tag>_all.json (test/item2_gates.js).
Writes validation/summary_item2.json, item2_energy.png, item2_window_rows.png."""
import json, os, sys, math, hashlib
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
V = os.path.join(ROOT, "validation")
TOL, CLOSE = 1e-4, 1e-6
BUILDS = [("pre", "backup_pre_item2/sim.js (item-1 final)"), ("final", "sim.js (item 2 final: strip + obstacle-aware separation clip)")]
COL = {"t": 0, "dt": 1, "dE": 2, "sum": 3, "pPE": 4, "pKE": 5, "pVR": 6, "strip": 7, "pushes": 8, "nStrip": 9, "bandAdv": 10, "bandSep": 11, "sep": 12}

def load(tag):
    f = f"/tmp/item2_{tag}_all.json"
    return json.load(open(f)) if os.path.exists(f) else None

def windows(t, x, w=0.5):
    """max over sliding windows (start at every substep start) and over fixed windows [k w, (k+1) w)."""
    c = np.concatenate([[0.0], np.cumsum(x)])
    t0 = np.concatenate([[0.0], t[:-1]])  # substep start times
    j = np.searchsorted(t, t0 + w + 1e-12, side="right")  # rows with t_end <= start + w
    slide = c[j] - c[np.arange(len(t))]
    k = np.minimum((t0 / w + 1e-9).astype(int), int(round(t[-1] / w)) - 1)
    fixed = np.bincount(k, weights=x, minlength=int(round(t[-1] / w)))
    return slide, fixed

def case_summary(r):
    s = np.array(r["series"]); t = s[:, 0]
    rows = r["rows"]
    pushrow = s[:, COL["pPE"]] + s[:, COL["pKE"]]
    removal = s[:, COL["pVR"]] + s[:, COL["strip"]]
    gb = pushrow + removal
    sl_b, fx_b = windows(t, gb); sl_p, fx_p = windows(t, pushrow)
    sep = s[:, 12] if s.shape[1] > 12 else None
    clos = s[:, COL["dE"]] - s[:, COL["sum"]]
    comb = {}
    for k, v in rows.items():
        if k == "P2G" or k.endswith("obstacleGhostFill"): comb["P2G + obstacle ghost fill"] = comb.get("P2G + obstacle ghost fill", 0) + v
        else: comb[k] = v
    out = {
        "rows_cumulative_J_per_kg": comb,
        "sum_rows": float(sum(rows.values())), "dE_particles": float(s[:, COL["dE"]].sum()),
        "closure_cumulative": float(s[:, COL["dE"]].sum() - sum(rows.values())), "closure_max_substep": float(np.abs(clos).max()),
        "gate_b": {"cumulative": float(gb.sum()), "max_sliding_0p5s": float(sl_b.max()), "max_fixed_0p5s": float(fx_b.max()),
                   "fixed_25_30s": [float(v) for v in fx_b[50:60]], "removal_cumulative": float(removal.sum()), "removal_max_substep": float(removal.max())},
        "gate_pushout_alone": {"cumulative": float(pushrow.sum()), "PE_cumulative": float(s[:, COL["pPE"]].sum()), "KE_cumulative": float(s[:, COL["pKE"]].sum()),
                               "max_sliding_0p5s": float(sl_p.max()), "max_fixed_0p5s": float(fx_p.max()), "fixed_25_30s": [float(v) for v in fx_p[50:60]]},
        "containment": {k: r[k] for k in ["insideMaxAtSubstepEnd", "minSdfOverR", "countChanged", "particles", "obstacleDeleted", "nan", "nanGuard", "pushStuck", "capHits", "pcgCapHits"]},
        "push": {"count": r["push"]["pushed"], "mean_depth_over_r": (r["push"]["depthSum"] / r["push"]["pushed"] / RADIUS) if r["push"]["pushed"] else 0.0,
                 "max_depth_over_r": r["push"]["depthMax"] / RADIUS},
        "separation_row": ({"cumulative": float(sep.sum()), "max_sliding_0p5s": float(windows(t, sep)[0].max()), "max_fixed_0p5s": float(windows(t, sep)[1].max()), "fixed_25_30s": [float(v) for v in windows(t, sep)[1][50:60]], "min_sliding_0p5s": float(windows(t, sep)[0].min())} if sep is not None else None),
        "count": r["particles"][1], "gate4b": r.get("gate4b"), "clip": r.get("clip"), "band": r["band"], "limiter": r.get("limiter"),
        "dt_mean": float(s[:, 1].mean()), "substeps": len(t), "E0": r["E0"], "E30": r["E30"], "finalStateHash": r["finalStateHash"],
    }
    out["_series"] = (t, gb, pushrow, removal, sl_b, sl_p, s)
    return out

RADIUS = None
def main():
    global RADIUS
    data = {tag: load(tag) for tag, _ in BUILDS}
    H = 1.0 / 128  # not used for depth; radius from minSdf scale below
    summ = {"run": "Run 2 (dam break past box / ramp), CFL 1, cap 8, damping 0, 30 s", "tolerance_J_per_kg": TOL, "closure_tol": CLOSE, "builds": {}}
    md5 = lambda p: hashlib.md5(open(os.path.join(ROOT, p), "rb").read()).hexdigest()
    for tag, desc in BUILDS:
        D = data[tag]
        if D is None: continue
        b = {"description": desc, "cases": {}}
        b["sim_md5"] = md5({"pre": "backup_pre_item2/sim.js", "nm": "validation/item2_runs/sim_strip_only.js", "final": "sim.js"}[tag])
        for kind in ["box", "ramp"]:
            r = D.get(kind)
            if r is None: continue
            RADIUS = r["push"]["depthMax"] and None
            RADIUS = R_PARTICLE
            cs = case_summary(r)
            if f"{kind}_plain" in D: cs["readonly_hash_match"] = D[f"{kind}_plain"]["finalStateHash"] == r["finalStateHash"]
            c5 = D.get(f"{kind}_cfl05")
            if c5 is not None:
                cs5 = case_summary(c5)
                ratio = cs["dt_mean"] / cs5["dt_mean"]
                p1, p5 = cs["gate_pushout_alone"]["PE_cumulative"], cs5["gate_pushout_alone"]["PE_cumulative"]
                cs["cfl05"] = {"dt_mean_ratio": ratio, "pushPE_cfl1": p1, "pushPE_cfl05": p5, "report_only": "push counts and depths only; no p computed",
                               "gate_b_cfl05": cs5["gate_b"], "pushout_alone_cfl05": cs5["gate_pushout_alone"], "push_cfl05": cs5["push"], "containment_cfl05": cs5["containment"]}
                cs.pop("_series_cfl05", None)
            b["cases"][kind] = cs
        # gates
        G = {}
        C = b["cases"]
        if C:
            g1 = all(c["gate_b"]["cumulative"] <= TOL and c["gate_b"]["max_sliding_0p5s"] <= TOL and c["gate_b"]["max_fixed_0p5s"] <= TOL and c["gate_b"]["removal_max_substep"] <= 0
                     and abs(c["closure_cumulative"]) <= CLOSE and c["closure_max_substep"] <= CLOSE for c in C.values())
            G["1_pushout_plus_removal"] = "PASS" if g1 else "FAIL"
            for k, c in C.items():
                g = c["gate_pushout_alone"]
                G[f"2_pushout_alone_{k}"] = "PASS" if (g["cumulative"] <= TOL and g["max_sliding_0p5s"] <= TOL and g["max_fixed_0p5s"] <= TOL) else "FAIL"
            g3 = all(c["containment"]["insideMaxAtSubstepEnd"] == 0 and not c["containment"]["countChanged"] and c["containment"]["obstacleDeleted"] == 0
                     and c["containment"]["nan"] == 0 and c["containment"]["nanGuard"] is None for c in C.values())
            G["3_containment_run2"] = "PASS" if g3 else "FAIL"
            if tag != "pre":
                L = {k: c["limiter"] for k, c in C.items()}
                G["4_new_code_float64_inside_code (tangential, outward) per case"] = {k: [l["tangIn"], l["outIn"]] for k, l in L.items()}
                G["4_clip_velocity_writes (particles whose velocity changed in separation)"] = {k: c["clip"]["velChanged"] for k, c in C.items()}
                G["4_verdict"] = "PASS" if all(l["tangIn"] == 0 and l["outIn"] == 0 for l in L.values()) and all(c["clip"]["velChanged"] == 0 for c in C.values()) else "FAIL"
                G["4_info_stored_float32_changes (not gated, rulings 2/3)"] = {k: [l["tangExt"], l["outExt"], l["maxRelTang"]] for k, l in L.items()}
                for k, c in C.items():
                    w = c["gate4b"]; sig = math.sqrt(w["sumVt2"]) * 2 ** -24 / math.sqrt(3); sx = math.sqrt(w["sig2exact"]); n = c["count"]
                    G[f"4b_{k}"] = {"verdict": "PASS" if abs(w["dKEt"]) <= 5 * sig else "FAIL", "signed_sum_dKEt_per_kg": w["dKEt"] / n, "sigma_per_kg": sig / n, "sigma_exact_per_kg": sx / n,
                                    "ratio_to_sigma": w["dKEt"] / sig if sig else None, "ratio_to_sigma_exact": w["dKEt"] / sx if sx else None, "writes": w["writes"], "also_within_5_sigma_exact": abs(w["dKEt"]) <= 5 * sx}
                    k2 = c["clip"]; sg = math.sqrt(k2["sig2X"]); sgy = math.sqrt(k2["sig2Y"])
                    G[f"clip_{k}"] = {"in_code_double": {"max_rel_magnitude_mismatch_of_pair_pushes": k2["maxMagDiffRel"], "max_pair_COM_shift_over_h": k2["maxComOverH"], "max_clipped_over_unclipped": k2["maxGrow"], "verdict": "PASS" if k2["maxMagDiffRel"] <= 1e-12 and k2["maxComOverH"] <= 1e-12 else "FAIL"},
                                      "stored_float32_per_pair": {"pairs": k2["pairs"], "clipped": k2["clipped"], "violations_of_2^-24(|xA|+|xB|)": k2["storedPairViol"], "max_ratio_to_allowance": k2["storedPairMaxRel"], "verdict": "PASS" if k2["storedPairViol"] == 0 else "FAIL"},
                                      "stored_COM_signed_sum_guard": {"sum_x": k2["sumX"], "sum_y": k2["sumY"], "sigma_x": sg, "sigma_y": sgy, "ratio_x": k2["sumX"] / sg if sg else None, "ratio_y": k2["sumY"] / sgy if sgy else None, "verdict": "PASS" if abs(k2["sumX"]) <= 5 * sg and abs(k2["sumY"]) <= 5 * sgy else "FAIL"}}
                    sr = c["separation_row"]
                    G[f"separation_row_{k}"] = {"verdict": "PASS" if sr["cumulative"] <= TOL and sr["max_sliding_0p5s"] <= TOL and sr["max_fixed_0p5s"] <= TOL else "FAIL", **{x: y for x, y in sr.items() if x != "fixed_25_30s"}, "fixed_25_30s_max": max(sr["fixed_25_30s"])}
        b["gates"] = G
        summ["builds"][tag] = b
    # plots
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    fig2, ax2 = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    sty = {"pre": ("tab:red", "pre (item-1 build)"), "final": ("tab:green", "item 2 final (strip + sep. clip)")}
    for j, kind in enumerate(["box", "ramp"]):
        for tag, _ in BUILDS:
            b = summ["builds"].get(tag)
            if not b or kind not in b["cases"]: continue
            t, gb, pr_, rem, slb, slp, s = b["cases"][kind]["_series"]
            col, lab = sty[tag]
            ax[0, j].plot(t, np.cumsum(pr_), color=col, label=lab + ": push-out row")
            ax[0, j].plot(t, np.cumsum(s[:, COL["pPE"]]), color=col, ls=":", label=lab + ": push-out PE")
            ax[1, j].plot(t, np.cumsum(gb), color=col, label=lab + ": push-out + removal")
            ax[1, j].plot(t, np.cumsum(rem), color=col, ls="--", label=lab + ": removal")
            t0 = np.concatenate([[0], t[:-1]])
            ax2[0, j].plot(t0, slp, color=col, label=lab)
            ax2[1, j].plot(t0, slb, color=col, label=lab)
        for a in (ax[0, j], ax[1, j]): a.axhline(TOL, color="k", lw=0.7, ls="--"); a.grid(alpha=.3)
        for a in (ax2[0, j], ax2[1, j]): a.axhline(TOL, color="k", lw=0.7, ls="--", label="gate +1e-4"); a.grid(alpha=.3)
        ax[0, j].set_title(f"{kind}: cumulative push-out row (gate 2)"); ax[1, j].set_title(f"{kind}: cumulative push-out + removal (gate b)")
        ax2[0, j].set_title(f"{kind}: push-out row, sliding 0.5 s window"); ax2[1, j].set_title(f"{kind}: push-out + removal, sliding 0.5 s window")
        ax[1, j].set_xlabel("t [s]"); ax2[1, j].set_xlabel("window start [s]")
    for a in list(ax.flat) + list(ax2.flat): a.set_ylabel("J/kg")
    ax[0, 0].legend(fontsize=7); ax[1, 0].legend(fontsize=7); ax2[0, 0].legend(fontsize=7)
    for a in ax2.flat: a.set_yscale("symlog", linthresh=1e-6)
    fig.suptitle("Item 2, Run 2 (CFL 1, 30 s): energy rows"); fig.tight_layout(); fig.savefig(os.path.join(V, "item2_energy.png"), dpi=110)
    fig2.suptitle("Item 2, Run 2: 0.5 s window rows (symlog; dashed = +1e-4 gate)"); fig2.tight_layout(); fig2.savefig(os.path.join(V, "item2_window_rows.png"), dpi=110)
    for b in summ["builds"].values():
        for c in b["cases"].values(): c.pop("_series", None)
    json.dump(summ, open(os.path.join(V, "summary_item2.json"), "w"), indent=1)
    print(json.dumps({k: v["gates"] for k, v in summ["builds"].items()}, indent=1))

R_PARTICLE = float(os.environ.get("R_PARTICLE", "0"))
if __name__ == "__main__": main()
