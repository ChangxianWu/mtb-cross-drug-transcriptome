"""Cross-dataset effect sizes and random-effects meta-analysis (script 06)."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import (OUT, PROC, load_kegg_modules, dataset_expr_log, module_scores,
                    hedges_g, cliffs_delta, dl_meta, bh_fdr,
                    plan_165673, plan_292332, plan_123294, plan_xlsx)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def se_of_g(g, na, nb):
    if any(v is None or np.isnan(v) for v in (g, na, nb)) or na < 2 or nb < 2:
        return np.nan
    return float(np.sqrt((na + nb) / (na * nb) + g**2 / (2 * (na + nb))))


def load_scores(gse):
    p = os.path.join(PROC, f"{gse}_modscores.csv")
    if not os.path.exists(p):
        return None
    return pd.read_csv(p, index_col=0)


def effect_rows(gse, sdf, scores, modules, plan):
    """plan: (drug_list, time_col_values, control_rule) provided by per-dataset functions."""
    rows = []
    for drug, t, ctrl_mask, drug_mask in plan:
        for mod in modules:
            a = scores.loc[drug_mask.reindex(scores.index).fillna(False), mod].dropna()
            b = scores.loc[ctrl_mask.reindex(scores.index).fillna(False), mod].dropna()
            g, p, nn = hedges_g(a, b)
            rows.append(dict(gse=gse, module=mod, drug=drug, time_h=t,
                             hedges_g=g, p_dataset=p, n_drug=nn[0], n_ctrl=nn[1],
                             cliffs_delta=cliffs_delta(a, b),
                             se=se_of_g(g, nn[0], nn[1])))
    return rows


# contrast-plan functions live in mtblib (shared by 06/07; a priori)


def mode_b_199263(modules):
    p = os.path.join(PROC, "GSE199263_lfc.csv")
    df = pd.read_csv(p)
    df["gene"] = df["gene"].astype(str).str.replace("gene:", "", regex=False)
    rows = []
    for stress in df["stress"].unique():
        sub = df[df["stress"].eq(stress)].dropna(subset=["log2FoldChange", "lfcSE"])
        sub = sub[sub["lfcSE"] > 0]
        for mod, genes in modules.items():
            # DESeq table gene column is Rv#### verbatim; case-insensitive double matching
            hit = sub[sub["gene"].isin(genes) | sub["gene"].str.upper().isin([g.upper() for g in genes])]
            if len(hit) < 3:
                continue
            zz = (hit["log2FoldChange"] / hit["lfcSE"]).values
            z_mod = float(zz.sum() / np.sqrt(len(zz)))
            pval = 2 * __import__("scipy").stats.norm.sf(abs(z_mod))
            rows.append(dict(stress=stress, module=mod, k_genes=len(hit),
                             mean_lfc=float(hit["log2FoldChange"].mean()),
                             stouffer_z=round(z_mod, 2), p=f"{pval:.2e}"))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "modeB_module_effects.csv"), index=False, encoding="utf-8-sig")
    return out


def main():
    modules, sha = load_kegg_modules()
    mod_names = list(modules.keys())
    all_rows = []
    jobs = [
        ("GSE165673", plan_165673),
        ("GSE292332", plan_292332),
        ("GSE123294_H37Rv", plan_123294),
        ("GSE123294_041195", plan_123294),
        ("GSE288598", lambda s: plan_xlsx(s, "LZD")),
        ("GSE287913", lambda s: plan_xlsx(s, "DLM")),
    ]
    for gse, planner in jobs:
        sdf = pd.read_csv(os.path.join(PROC, f"{gse}_samples.tsv"), sep="\t").set_index("gsm")
        scores = load_scores(gse)
        if scores is None:
            print(f"[SKIP] {gse}: no module scores")
            continue
        scores = scores.reindex(sdf.index.intersection(scores.index))
        plan = planner(sdf)
        rows = effect_rows(gse, sdf, scores, mod_names, plan)
        all_rows.extend(rows)
        print(f"[OK] {gse}: {len(rows)} effect rows ({len(plan)} contrasts)")
    eff = pd.DataFrame(all_rows)
    eff.to_csv(os.path.join(OUT, "effects_per_dataset.csv"), index=False, encoding="utf-8-sig")

    # random-effects meta: pooled by (module, drug, time)
    metas = []
    for (mod, drug, t), grp in eff.dropna(subset=["hedges_g"]).groupby(["module", "drug", "time_h"]):
        m = dl_meta(grp["hedges_g"].values, grp["se"].values)
        m["cliff_median"] = float(grp["cliffs_delta"].median())
        metas.append(dict(module=mod, drug=drug, time_h=t, datasets=", ".join(grp["gse"]),
                          **{k: (round(v, 3) if isinstance(v, float) else v) for k, v in m.items()}))
    meta = pd.DataFrame(metas)
    if len(meta):
        meta["q_bh"] = bh_fdr(meta["p"].values)
        meta = meta.sort_values(["q_bh", "g_dl"], ascending=[True, False])
    meta.to_csv(os.path.join(OUT, "meta_module_effects.csv"), index=False, encoding="utf-8-sig")

    modeb = mode_b_199263(modules)
    print(f"\neffect rows {len(eff)}; meta groups {len(meta)}; Mode B rows {len(modeb)}")
    print("Top significant meta groups (q_bh<0.05):")
    if len(meta):
        print(meta[meta["q_bh"] < 0.05].head(25).to_string(index=False))
    print(f"\noutput -> {OUT}")


if __name__ == "__main__":
    main()
