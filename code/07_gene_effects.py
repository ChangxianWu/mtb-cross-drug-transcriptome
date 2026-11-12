"""Gene-level effect layer (script 07, Mode A2 robustness arm)."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import (OUT, PROC, load_kegg_modules, dataset_expr_log, dl_meta, bh_fdr,
                    plan_165673, plan_292332, plan_123294, plan_xlsx)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

JOBS = [
    ("GSE165673", plan_165673),
    ("GSE292332", plan_292332),
    ("GSE123294_H37Rv", plan_123294),
    ("GSE123294_041195", plan_123294),
    ("GSE288598", lambda s: plan_xlsx(s, "LZD")),
    ("GSE287913", lambda s: plan_xlsx(s, "DLM")),
]


def gene_effects(expr, cols_a, cols_b):
    A = expr[cols_a]
    B = expr[cols_b]
    na, nb = A.shape[1], B.shape[1]
    ma, mb = A.mean(axis=1), B.mean(axis=1)
    sa, sb = A.std(axis=1, ddof=1), B.std(axis=1, ddof=1)
    sp2 = ((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = (ma - mb) / np.sqrt(sp2)
    J = 1 - 3 / (4 * (na + nb) - 9)
    g = d * J
    v = sa**2 / na + sb**2 / nb
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ma - mb) / np.sqrt(v)
        df = v**2 / ((sa**2 / na)**2 / (na - 1) + (sb**2 / nb)**2 / (nb - 1))
    p = 2 * stats.t.sf(np.abs(t), df)
    return pd.DataFrame({"g": g, "p": p})


def main():
    modules, sha = load_kegg_modules()
    rows = []
    for gse, planner in JOBS:
        expr, sdf = dataset_expr_log(gse)
        if expr is None:
            print(f"[SKIP] {gse}")
            continue
        sdf = sdf.set_index("gsm")
        plan = planner(sdf)
        for drug, t, veh_m, drug_m in plan:
            cols_a = [c for c in expr.columns if drug_m.get(c, False)]
            cols_b = [c for c in expr.columns if veh_m.get(c, False)]
            if len(cols_a) < 2 or len(cols_b) < 2:
                continue
            ge = gene_effects(expr, cols_a, cols_b)
            ge = ge.dropna(subset=["g"])
            for mod, genes in modules.items():
                gset = set(genes)
                hit = ge[ge.index.astype(str).isin(gset) | ge.index.str.upper().isin({g.upper() for g in genes})]
                k = len(hit)
                if k < 2:
                    rows.append(dict(gse=gse, module=mod, drug=drug, time_h=t, g_gene_mean=np.nan,
                                     se_between=np.nan, k_genes=k, frac_sign=np.nan, k_zero_var=0))
                    continue
                gm = hit["g"].mean()
                se_b = hit["g"].std(ddof=1) / np.sqrt(k) if k > 1 else np.nan
                frac = float((np.sign(hit["g"]) == np.sign(gm)).mean())
                rows.append(dict(gse=gse, module=mod, drug=drug, time_h=t,
                                 g_gene_mean=gm, se_between=se_b, k_genes=k,
                                 frac_sign=frac, k_zero_var=0))
        print(f"[OK] {gse}: {len(plan)} gene-level contrasts done")
    eff = pd.DataFrame(rows)
    eff.to_csv(os.path.join(OUT, "gene_effects_per_dataset.csv"), index=False, encoding="utf-8-sig")

    # A2 meta
    metas = []
    for (mod, drug, t), grp in eff.dropna(subset=["g_gene_mean"]).groupby(["module", "drug", "time_h"]):
        if (grp["se_between"] <= 0).any() or grp["se_between"].isna().any():
            k1 = len(grp)
            g = float(grp["g_gene_mean"].mean())
            se = float(grp["g_gene_mean"].std(ddof=1) / np.sqrt(k1)) if k1 > 1 else np.nan
        else:
            m = dl_meta(grp["g_gene_mean"].values, grp["se_between"].values)
            g, se = m["g_dl"], m["se"]
        z = g / se if (se and se > 0) else np.nan
        p = 2 * stats.norm.sf(abs(z)) if not np.isnan(z) else np.nan
        metas.append(dict(module=mod, drug=drug, time_h=t, k=len(grp),
                          datasets=", ".join(grp["gse"]), g_a2=round(g, 3), se_a2=round(se, 3) if se == se else np.nan,
                          p_a2=p))
    meta2 = pd.DataFrame(metas)
    meta2["q_a2"] = bh_fdr(meta2["p_a2"].fillna(1.0).values)
    meta2 = meta2.sort_values("q_a2")
    meta2.to_csv(os.path.join(OUT, "gene_meta.csv"), index=False, encoding="utf-8-sig")

    # A1 vs A2 direction agreement
    a1 = pd.read_csv(os.path.join(OUT, "meta_module_effects.csv"))
    key = ["module", "drug", "time_h"]
    mg = a1.merge(meta2, on=key, how="inner", suffixes=("_a1", "_a2"))
    mg = mg.dropna(subset=["g_dl", "g_a2"])
    mg["direction_match"] = np.sign(mg["g_dl"]) == np.sign(mg["g_a2"])
    conc = mg.groupby("direction_match").size()
    match_rate = float(mg["direction_match"].mean()) if len(mg) else np.nan
    sig = mg[mg["q_bh"] < 0.05]
    sig_rate = float(sig["direction_match"].mean()) if len(sig) else np.nan
    mg.to_csv(os.path.join(OUT, "a1_a2_concordance.csv"), index=False, encoding="utf-8-sig")
    print(f"A1xA2 comparable groups {len(mg)}: direction agreement {match_rate:.3f} (A1 q_bh<0.05 subset: {sig_rate})")
    print(conc.to_string())
    print(f"output -> {OUT}")


if __name__ == "__main__":
    main()
