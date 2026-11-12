"""Figures F1-F8 (script 11), 300 dpi, English labels throughout."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import OUT, PROC, ROOT

# internal module ids -> manuscript display names (capitalisation matches the text:
# DosR markers / iniBAC / TCA cycle ...)
MODULE_DISPLAY = {
    "respiration_oxphos": "Respiration (OXPHOS)",
    "dosr_markers": "DosR markers",
    "inibac": "iniBAC",
    "ribosome": "Ribosome",
    "fatty_acid_biosyn": "Fatty-acid biosynthesis",
    "fatty_acid_degrad": "Fatty-acid degradation",
    "tca_cycle": "TCA cycle",
    "glycolysis": "Glycolysis",
    "pentose_phosphate": "Pentose phosphate",
    "peptidoglycan": "Peptidoglycan",
    "trehalose_glycerolipid": "Trehalose-glycerolipid",
    "two_component_system": "Two-component system",
    "abc_transporters": "ABC transporters",
    "sulfur_metabolism": "Sulfur metabolism",
    "nitrogen_metabolism": "Nitrogen metabolism",
}


def mdisp(mid):
    return MODULE_DISPLAY.get(str(mid), str(mid).replace("_", " "))


FIG = os.path.join(ROOT, "work", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "figure.dpi": 300})


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.tight_layout()
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[FIG] {name}")


def f1_data_map():
    inv = pd.read_csv(os.path.join(OUT, "GEO_INVENTORY.tsv"), sep="\t")
    n = inv.groupby("gse").agg(samples=("n_samples", "first"), platform=("platform", "first")).reset_index()
    # only the six analysed series (consistent with Table 1 and the caption; excluded/deferred series omitted)
    order = ["GSE165673", "GSE292332", "GSE292331", "GSE123294", "GSE288598", "GSE287913"]
    n = n[n["gse"].isin(order)].copy()
    n["o"] = n["gse"].apply(lambda g: order.index(g))
    n = n.sort_values("o")
    fig, ax = plt.subplots(figsize=(5.8, 2.6))
    bars = ax.barh(n["gse"], pd.to_numeric(n["samples"], errors="coerce").fillna(0),
                   color="#2c7fb8", height=0.55)
    ax.set_xlabel("Reported samples (GEO)")
    for b, s in zip(bars, n["samples"]):
        ax.text(b.get_width() + 3, b.get_y() + b.get_height() / 2, str(s), va="center", fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, float(pd.to_numeric(n["samples"], errors="coerce").fillna(0).max()) * 1.15)
    save(fig, "F1_data_map.png")


def f2_positive_control():
    r = pd.read_csv(os.path.join(OUT, "positive_control_report.csv"))
    r = r[r["verdict"].notna()]
    labels = [f"{mdisp(row['module'])}\n{row['contrast']}" for _, row in r.iterrows()]
    colors = ["#d73027" if v in ("PASS",) else ("#43a2ca" if v == "DIR_OK" else "#bdbdbd") for v in r["verdict"]]
    fig, ax = plt.subplots(figsize=(6.6, 4))
    ax.barh(range(len(r)), r["hedges_g"], color=colors)
    ax.set_yticks(range(len(r)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Hedges g (module score, drug/stress vs control)")
    ax.set_title("")
    # pad both sides so verdict labels next to negative bars do not overlap the y-axis labels
    lo = float(min(r["hedges_g"].min(), 0)); hi = float(max(r["hedges_g"].max(), 0))
    pad = 0.16 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    for i, (g, v) in enumerate(zip(r["hedges_g"], r["verdict"])):
        ax.text(g + (0.02 * (hi - lo) if g > 0 else -0.02 * (hi - lo)), i, str(v),
                va="center", ha="left" if g > 0 else "right", fontsize=6.5)
    ax.invert_yaxis()
    save(fig, "F2_positive_control.png")


def f3_meta_forest():
    m = pd.read_csv(os.path.join(OUT, "meta_module_effects.csv"))
    m = m[m["k"] >= 2].copy()
    m["abs_g"] = m["g_dl"].abs()
    m = m.sort_values("abs_g", ascending=False).head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    y = range(len(m))
    ax.errorbar(m["g_dl"], y, xerr=1.96 * m["se"], fmt="o", ms=3, lw=0.8, color="#252525", capsize=2)
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{mdisp(r.module)} | {r.drug} {int(r.time_h)}h (k={r.k})" for r in m.itertuples()], fontsize=6.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Random-effects Hedges g (95% CI), cross-dataset merges (k \u2265 2)")
    ax.set_title("")
    save(fig, "F3_meta_forest.png")


def f4_module_drug_heatmap():
    p = pd.read_csv(os.path.join(OUT, "module_by_drug_primary_effects.csv"))
    piv = p.pivot_table(index="module", columns="drug", values="g_a2")
    q = p.pivot_table(index="module", columns="drug", values="q_a2")
    kcd = {d: int(p.loc[p.drug.eq(d), "kill_class"].iloc[0]) for d in piv.columns}
    order_c = sorted(piv.columns, key=lambda d: (-kcd[d], d))
    piv, q = piv[order_c], q[order_c]
    fig, ax = plt.subplots(figsize=(6.7, 4.5))
    im = ax.imshow(piv.values, cmap="coolwarm", vmin=-10, vmax=10, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(list(piv.columns), rotation=0, fontsize=6.5)
    kc = [kcd[c] for c in piv.columns]
    # kill-class groups: white separators + group labels (kept out of tick labels to avoid crowding)
    for j in range(1, len(kc)):
        if kc[j] != kc[j - 1]:
            ax.axvline(j - 0.5, color="white", lw=2.2)
    starts = [0] + [j for j in range(1, len(kc)) if kc[j] != kc[j - 1]]
    ends = [j - 1 for j in starts[1:]] + [len(kc) - 1]
    names = {1: "kill=1", 0: "kill=0", -1: "kill=n.d."}
    for s, e in zip(starts, ends):
        ax.text((s + e) / 2, -0.75, names.get(kc[s], f"kill={kc[s]}"), ha="center", va="center", fontsize=6.5)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([mdisp(i) for i in piv.index], fontsize=6.5)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                star = "*" if (q.values[i, j] < 0.05) else ""
                ax.text(j, i, f"{v:.1f}{star}", ha="center", va="center", fontsize=5.0,
                        color="white" if abs(v) > 6 else "black")
    ax.set_title("")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Hedges g (A2)")
    save(fig, "F4_module_drug_heatmap.png")


def f5_timecourse():
    e = pd.read_csv(os.path.join(OUT, "effects_per_dataset.csv"))
    e = e[(e.gse == "GSE165673") & (e.drug != "UNKNOWN")]
    mods = ["respiration_oxphos", "dosr_markers", "inibac", "ribosome", "fatty_acid_biosyn"]
    drugs = ["BDQ", "INH", "RIF", "PTN", "LZD", "CFZ"]
    # two-row layout (3 + 2), larger fonts; legend in the 6th cell, no suptitle
    fig, axes = plt.subplots(2, 3, figsize=(6.7, 4.1), sharey=False)
    axes = axes.ravel()
    for i, mod in enumerate(mods):
        ax = axes[i]
        for d in drugs:
            s = e[(e.module == mod) & (e.drug == d)].sort_values("time_h")
            if len(s):
                ax.plot(s["time_h"], s["hedges_g"], marker="o", ms=3.5, lw=1.1, label=d)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(mdisp(mod), fontsize=9)
        ax.set_xticks([4, 24, 72])
        ax.tick_params(labelsize=7.5)
        ax.set_xlabel("Time (h)", fontsize=8)
        if i % 3 == 0:
            ax.set_ylabel("Hedges g vs vehicle", fontsize=8)
    leg_ax = axes[len(mods)]
    leg_ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    leg_ax.legend(handles, labels, fontsize=8, ncol=2, frameon=False, loc="center")
    save(fig, "F5_timecourse.png")


def f6_hypoxia_validation():
    sc = pd.read_csv(os.path.join(PROC, "GSE292331_modscores.csv"), index_col=0)
    sdf = pd.read_csv(os.path.join(PROC, "GSE292331_samples.tsv"), sep="\t").set_index("gsm")
    sdf["day"] = sdf["time_h"].str.replace("d", "", regex=False).astype(float)
    sc = sc.reindex(sdf.index)
    sc["day"] = sdf["day"]
    agg = sc.groupby("day").mean()
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    for mod in ["respiration_oxphos", "dosr_markers", "inibac", "tca_cycle"]:
        ax.plot(agg.index, agg[mod], marker="o", ms=3, lw=1.1, label=mdisp(mod))
    ax.axvspan(0, 0.5, alpha=0.08, color="grey")
    ax.axvline(7, color="grey", ls="--", lw=0.7)
    ax.text(3.5, ax.get_ylim()[1] * 0.92, "Hypoxia", ha="center", fontsize=8)
    ax.text(9.8, ax.get_ylim()[1] * 0.92, "Reaeration", ha="center", fontsize=8)
    ax.set_xlabel("Day (Wayne-type hypoxia\u2013reaeration, GSE292331)")
    ax.set_ylabel("Module score (dataset-level z mean)")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    save(fig, "F6_hypoxia_validation.png")


def f7_sensitivity():
    vals, labels = [], []
    loo = pd.read_csv(os.path.join(OUT, "sens_loo.csv"))
    vals.append(loo[loo.layer == "A1"]["direction_stable"].mean()); labels.append("LOO (A1)")
    vals.append(loo[loo.layer == "A2"]["direction_stable"].mean()); labels.append("LOO (A2)")
    e = pd.read_csv(os.path.join(OUT, "sens_dose.csv"))
    ag = []
    for _, grp in e.groupby(["module", "drug", "time_h"]):
        if len(grp) == 2:
            ag.append(np.sign(grp["hedges_g"].iloc[0]) == np.sign(grp["hedges_g"].iloc[1]))
    vals.append(float(np.mean(ag))); labels.append("dose strata")
    r = pd.read_csv(os.path.join(OUT, "sens_rankscore.csv"))
    a1 = pd.read_csv(os.path.join(OUT, "effects_per_dataset.csv")).rename(columns={"hedges_g": "g_z"})
    mg = a1.merge(r, on=["gse", "module", "drug", "time_h"]).dropna(subset=["g_z", "hedges_g_rank"])
    vals.append(float((np.sign(mg.g_z) == np.sign(mg.hedges_g_rank)).mean())); labels.append("rank scores")
    cc = pd.read_csv(os.path.join(OUT, "a1_a2_concordance.csv"))
    vals.append(float(cc[cc["q_bh"] < 0.05]["direction_match"].mean())); labels.append("A1 vs A2 (sig.)")
    fig, ax = plt.subplots(figsize=(5.6, 3))
    bars = ax.bar(labels, [v * 100 for v in vals], color="#2c7fb8")
    ax.axhline(50, color="grey", ls=":", lw=0.8)
    ax.set_ylabel("direction agreement (%)")
    ax.set_ylim(0, 105)
    ax.set_title("")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v * 100 + 2, f"{v*100:.1f}", ha="center", fontsize=7.5)
    ax.tick_params(axis="x", labelrotation=20)
    save(fig, "F7_sensitivity.png")


def f8_a1a2_scatter():
    c = pd.read_csv(os.path.join(OUT, "a1_a2_concordance.csv")).dropna(subset=["g_dl", "g_a2"])
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    sig = c["q_bh"] < 0.05
    ax.scatter(c.loc[~sig, "g_dl"], c.loc[~sig, "g_a2"], s=6, alpha=0.35, color="#969696", label="q \u2265 0.05")
    ax.scatter(c.loc[sig, "g_dl"], c.loc[sig, "g_a2"], s=9, alpha=0.75, color="#d73027", label="q < 0.05")
    lim = 30
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_ylabel("A2 g (gene-level mean effect)")
    ax.set_title("")
    ax.legend(fontsize=7)
    ax.set_xlabel("A1 Hedges g (sample-level module SMD)")
    save(fig, "F8_a1a2_scatter.png")


def main():
    f1_data_map()
    f2_positive_control()
    f3_meta_forest()
    f4_module_drug_heatmap()
    f5_timecourse()
    f6_hypoxia_validation()
    f7_sensitivity()
    f8_a1a2_scatter()
    print(f"\nfigures → {FIG}")


if __name__ == "__main__":
    main()
