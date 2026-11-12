"""Negative controls (script 08): size-matched random gene sets form the empirical null."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import (OUT, PROC, dataset_expr_log, plan_165673, plan_292332, plan_123294,
                    plan_xlsx)

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

R = 100
SIZES = [3, 30, 60]
SEED = 20260918
JOBS = [
    ("GSE165673", plan_165673),
    ("GSE292332", plan_292332),
    ("GSE123294_H37Rv", plan_123294),
    ("GSE123294_041195", plan_123294),
    ("GSE288598", lambda s: plan_xlsx(s, "LZD")),
    ("GSE287913", lambda s: plan_xlsx(s, "DLM")),
]


def smd_matrix(S, ia, ib):
    """S: samples x R random-set scores; returns the per-permutation Hedges g vector."""
    A, B = S[ia, :], S[ib, :]
    na, nb = A.shape[0], B.shape[0]
    ma, mb = A.mean(axis=0), B.mean(axis=0)
    sa, sb = A.std(axis=0, ddof=1), B.std(axis=0, ddof=1)
    sp = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    J = 1 - 3 / (4 * (na + nb) - 9)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (ma - mb) / sp * J


def build_random_scores(Zv, genes, rng):
    S_by_size = {}
    for size in SIZES:
        M = np.zeros((len(genes), R))
        for r in range(R):
            M[rng.choice(len(genes), size=size, replace=False), r] = 1.0
        S_by_size[size] = Zv @ M / size
    return S_by_size


def main():
    rng = np.random.default_rng(SEED)
    null_rows, merged = [], []
    eff = pd.read_csv(os.path.join(OUT, "effects_per_dataset.csv"))
    hits = pd.read_csv(os.path.join(OUT, "module_hits.csv"), index_col=0).T  # index=module, cols=gse
    for gse, planner in JOBS:
        expr, sdf = dataset_expr_log(gse)
        if expr is None:
            continue
        sdf = sdf.set_index("gsm")
        expr = expr[[c for c in expr.columns if c in sdf.index]]
        genes = list(expr.index)
        Z = expr.T
        Z = (Z - Z.mean()) / Z.std(ddof=0).replace(0, np.nan)
        Zv = Z.values
        S_by_size = build_random_scores(Zv, genes, rng)
        plan = planner(sdf)
        plan_lookup = {(d, float(t)): (vm, dm) for d, t, vm, dm in plan}
        sub = eff[eff["gse"].eq(gse)].copy()
        out_q = []
        for drug, t, veh_m, drug_m in plan:
            ia = np.where(drug_m.reindex(Z.index).fillna(False).values)[0]
            ib = np.where(veh_m.reindex(Z.index).fillna(False).values)[0]
            if len(ia) < 2 or len(ib) < 2:
                continue
            for size in SIZES:
                gp = smd_matrix(S_by_size[size], ia, ib)
                gp = gp[np.isfinite(gp)]
                if len(gp) == 0:
                    continue
                null_rows.append(dict(gse=gse, drug=drug, time_h=t, size=size, R=R,
                                      q95=float(np.quantile(np.abs(gp), 0.95)),
                                      q99=float(np.quantile(np.abs(gp), 0.99)),
                                      median_abs=float(np.median(np.abs(gp)))))
        for _, row in sub.iterrows():
            key = (row["drug"], float(row["time_h"]))
            if key not in plan_lookup:
                out_q.append(np.nan)
                continue
            vm, dm = plan_lookup[key]
            ia = np.where(dm.reindex(Z.index).fillna(False).values)[0]
            ib = np.where(vm.reindex(Z.index).fillna(False).values)[0]
            size = min(SIZES, key=lambda s: abs(s - float(hits.loc[row["module"], gse]))
                       if (gse in hits.columns and row["module"] in hits.index) else 30)
            gp = smd_matrix(S_by_size[size], ia, ib)
            gp = gp[np.isfinite(gp)]
            gobs = abs(float(row["hedges_g"])) if pd.notna(row["hedges_g"]) else np.nan
            q = (1 + np.sum(np.abs(gp) >= gobs)) / (1 + len(gp)) if np.isfinite(gobs) else np.nan
            out_q.append(round(float(q), 4))
        sub["q_emp"] = out_q
        merged.append(sub)
        print(f"[OK] {gse}: null {len(null_rows)} rows accumulated, observed q_emp {sum(1 for q in out_q if q == q)}")
    pd.DataFrame(null_rows).to_csv(os.path.join(OUT, "negative_control_summary.csv"),
                                   index=False, encoding="utf-8-sig")
    effq = pd.concat(merged, ignore_index=True)
    effq.to_csv(os.path.join(OUT, "effects_with_qemp.csv"), index=False, encoding="utf-8-sig")
    valid = effq["q_emp"].dropna()
    frac = float((valid <= 0.05).mean())
    print(f"\nobserved module effects with q_emp<=0.05: {frac:.3f} (random expectation ~0.05; markedly higher = compositional-inflation warning)")
    print(f"output -> {OUT}")


if __name__ == "__main__":
    main()
