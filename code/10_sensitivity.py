"""Sensitivity suite (script 10)."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import (OUT, PROC, dataset_expr_log, load_kegg_modules, module_scores,
                    dl_meta, hedges_g, plan_165673, plan_292332)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def se_of_g(g, na, nb):
    if any((isinstance(v, float) and np.isnan(v)) for v in (g,)) or na < 2 or nb < 2:
        return np.nan
    return float(np.sqrt((na + nb) / (na * nb) + g**2 / (2 * (na + nb))))


def s1_loo():
    rows = []
    for tag, f, gcol in [("A1", "effects_per_dataset.csv", "hedges_g"),
                         ("A2", "gene_effects_per_dataset.csv", "g_gene_mean")]:
        e = pd.read_csv(os.path.join(OUT, f)).dropna(subset=[gcol, "se_between" if tag == "A2" else "se"])
        secol = "se_between" if tag == "A2" else "se"
        for (mod, drug, t), grp in e.groupby(["module", "drug", "time_h"]):
            if len(grp) < 2:
                continue
            ds = list(grp["gse"])
            ys = grp[gcol].values.astype(float)
            ses = grp[secol].values.astype(float)
            ok = np.isfinite(ys) & np.isfinite(ses) & (ses > 0)
            if ok.sum() < 2:
                continue
            ds, ys, ses = [d for d, o in zip(ds, ok) if o], ys[ok], ses[ok]
            base_sign = np.sign(float(dl_meta(ys, ses)["g_dl"]))
            stable, loo_gs = True, []
            for i in range(len(ds)):
                keep = [j for j in range(len(ds)) if j != i]
                mm = dl_meta(ys[keep], ses[keep])
                loo_gs.append(mm["g_dl"])
                if np.sign(mm["g_dl"]) != base_sign:
                    stable = False
            rows.append(dict(layer=tag, module=mod, drug=drug, time_h=t, k=len(ds),
                             direction_stable=stable,
                             loo_g_min=round(min(loo_gs), 3), loo_g_max=round(max(loo_gs), 3)))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "sens_loo.csv"), index=False, encoding="utf-8-sig")
    rate = out.groupby("layer")["direction_stable"].mean()
    print(f"S1 LOO: overall stability {out['direction_stable'].mean():.3f};"
          f"A1={rate.get('A1', np.nan):.3f}, A2={rate.get('A2', np.nan):.3f}({len(out)} groups)")
    return out


def s2_dose():
    expr, sdf = dataset_expr_log("GSE165673")
    modules, _ = load_kegg_modules()
    scores, hits = module_scores(expr, modules)
    sdf = sdf.set_index("gsm")
    sdf["dose_num"] = pd.to_numeric(sdf["dose_code"], errors="coerce")
    sdf["time_h"] = pd.to_numeric(sdf["time_h"], errors="coerce")
    scores = scores.reindex(sdf.index)
    rows = []
    for drug in sorted(set(sdf.loc[sdf["dose_num"] > 0, "drug"])):
        dvals = sorted(sdf.loc[sdf["drug"].eq(drug) & sdf["dose_num"].gt(0), "dose_num"].unique())
        for t in [4, 24, 72]:
            veh = scores.index[sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].eq(0)]
            for dlevel, dval in [("low", dvals[0]), ("high", dvals[-1])]:
                tr = scores.index[sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].eq(dval)]
                if len(tr) < 2 or len(veh) < 2:
                    continue
                for mod in modules:
                    g, p, nn = hedges_g(scores.loc[tr, mod], scores.loc[veh, mod])
                    rows.append(dict(module=mod, drug=drug, time_h=t, dose_level=dlevel,
                                     hedges_g=g, n_drug=nn[0]))
    df = pd.DataFrame(rows)
    # direction agreement: low/high dose sign consistency within (module, drug, time)
    agree = []
    for (mod, drug, t), grp in df.groupby(["module", "drug", "time_h"]):
        if len(grp) == 2:
            agree.append(np.sign(grp["hedges_g"].iloc[0]) == np.sign(grp["hedges_g"].iloc[1]))
    rate = float(np.mean(agree)) if agree else np.nan
    df.to_csv(os.path.join(OUT, "sens_dose.csv"), index=False, encoding="utf-8-sig")
    print(f"S2 dose strata: low/high dose agreement {rate:.3f} ({len(agree)} pairs)")
    return df


def s3_rankscore():
    rows = []
    for gse, planner in [("GSE165673", plan_165673), ("GSE292332", plan_292332)]:
        expr, sdf = dataset_expr_log(gse)
        sdf = sdf.set_index("gsm")
        expr = expr[[c for c in expr.columns if c in sdf.index]]
        modules, _ = load_kegg_modules()
        # rank scoring: genes ranked within each sample (column-wise); module score = mean module-gene rank - genome mean rank (per sample)
        rk = expr.rank(axis=0)
        gmean_rank = rk.mean(axis=0)
        scores = pd.DataFrame(index=rk.columns)
        for mod, genes in modules.items():
            gset = set(g.upper() for g in genes)
            hit = [g for g in rk.index if str(g).upper() in gset]
            scores[mod] = rk.loc[hit].mean(axis=0) - gmean_rank if hit else np.nan
        scores = scores.reindex(sdf.index)
        for drug, t, veh_m, drug_m in planner(sdf):
            ia = scores.index[drug_m.reindex(scores.index).fillna(False)]
            ib = scores.index[veh_m.reindex(scores.index).fillna(False)]
            if len(ia) < 2 or len(ib) < 2:
                continue
            for mod in modules:
                g, p, nn = hedges_g(scores.loc[ia, mod], scores.loc[ib, mod])
                rows.append(dict(gse=gse, module=mod, drug=drug, time_h=t,
                                 hedges_g_rank=g, n_drug=nn[0]))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "sens_rankscore.csv"), index=False, encoding="utf-8-sig")
    a1 = pd.read_csv(os.path.join(OUT, "effects_per_dataset.csv")).rename(
        columns={"hedges_g": "hedges_g_z"})
    mg = a1.merge(df, on=["gse", "module", "drug", "time_h"]).dropna(
        subset=["hedges_g_z", "hedges_g_rank"])
    rate = float((np.sign(mg["hedges_g_z"]) == np.sign(mg["hedges_g_rank"])).mean())
    print(f"S3 rank-score alternative: agreement with z scores {rate:.3f} ({len(mg)} pairs)")
    return df


def s4_timecons():
    e = pd.read_csv(os.path.join(OUT, "effects_per_dataset.csv"))
    e = e[e["gse"].eq("GSE165673") & e["time_h"].isin([4, 24, 72])].dropna(subset=["hedges_g"])
    rows = []
    for (mod, drug), grp in e.groupby(["module", "drug"]):
        if len(grp) < 2:
            continue
        signs = np.sign(grp["hedges_g"]).astype(int)
        dom = signs.value_counts().iloc[0] / len(signs)
        rows.append(dict(module=mod, drug=drug, n_times=len(grp), dominant_frac=float(dom)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "sens_timecons.csv"), index=False, encoding="utf-8-sig")
    print(f"S4 cross-time: three-timepoint direction consistency {(df['dominant_frac'] == 1.0).mean():.3f} ({len(df)} groups)")
    return df


def main():
    s1_loo()
    s2_dose()
    s3_rankscore()
    s4_timecons()
    print(f"\noutput -> {OUT}")


if __name__ == "__main__":
    main()
