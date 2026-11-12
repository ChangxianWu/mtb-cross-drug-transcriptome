"""Module-by-drug kill-class association (script 09)."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import OUT

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# kill classification (anchored to probe records; simplification stated in docstring)
CLASSIFICATION = [
    ("INH", 1, ""), ("RIF", 1, ""), ("RPE", 1, ""), ("BDQ", 1, ""),
    ("PTN", 1, ""), ("MXF", 1, ""), ("LVF", 1, ""), ("CPF", 1, ""),
    ("AMK", 1, ""), ("POA", 1, "active form of PZA; bactericidal at acidic pH"),
    ("LZD", 1, "EBA positive but dose-dependent; caveat"),
    ("CFZ", 0, "mainly bacteriostatic, low EBA"), ("DLM", 0, "bacteriostatic/bactericidal mixed in literature; uncertainty caveat; sensitivity flip"),
]


def main():
    cls = pd.DataFrame(CLASSIFICATION, columns=["drug", "kill_class", "note"])
    cls.to_csv(os.path.join(OUT, "drug_classification.csv"), index=False, encoding="utf-8-sig")

    a2 = pd.read_csv(os.path.join(OUT, "gene_meta.csv"))
    effq = pd.read_csv(os.path.join(OUT, "effects_with_qemp.csv"))
    # per (module, drug): A2 |g|-maximal time point
    a2_nn = a2.dropna(subset=["g_a2"])
    idx = a2_nn.groupby(["module", "drug"])["g_a2"].apply(lambda s: s.abs().idxmax())
    prim = a2_nn.loc[idx][["module", "drug", "time_h", "g_a2", "q_a2", "k", "datasets"]].copy()
    prim = prim.merge(cls, on="drug", how="left")
    prim["kill_class"] = prim["kill_class"].fillna(-1).astype(int)  # -1 = unclassified (should not occur)
    # q_emp corroboration (minimum across time points per module x drug)
    qmin = effq.groupby(["module", "drug"])["q_emp"].min().rename("q_emp_min")
    prim = prim.merge(qmin, on=["module", "drug"], how="left")

    # module x kill-class association: two-group Mann-Whitney + point-biserial (kill-0 has only CFZ/DLM; limited power - descriptive)
    tests = []
    for mod, grp in prim.groupby("module"):
        g1 = grp.loc[grp["kill_class"].eq(1), "g_a2"]
        g0 = grp.loc[grp["kill_class"].eq(0), "g_a2"]
        if len(g1) >= 3 and len(g0) >= 1:
            u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
        else:
            u, p = np.nan, np.nan
        r, pr = stats.pointbiserialr(grp["kill_class"], grp["g_a2"]) if grp["kill_class"].nunique() > 1 else (np.nan, np.nan)
        tests.append(dict(module=mod, n_kill1=len(g1), n_kill0=len(g0),
                          median_g_kill1=round(float(g1.median()), 3) if len(g1) else np.nan,
                          median_g_kill0=round(float(g0.median()), 3) if len(g0) else np.nan,
                          u=u, p_mw=p, r_pb=round(float(r), 3), p_pb=pr))
    tests_df = pd.DataFrame(tests).sort_values("p_mw")
    tests_df.to_csv(os.path.join(OUT, "killclass_tests.csv"), index=False, encoding="utf-8-sig")
    prim.to_csv(os.path.join(OUT, "module_by_drug_primary_effects.csv"), index=False, encoding="utf-8-sig")

    print("Module x kill class (top Mann-Whitney results per module):")
    print(tests_df.head(10).to_string(index=False))
    print("\nDLM flip caveat: kill-0 has only CFZ/DLM, limited Mann-Whitney power - descriptive only, no causal or strong statistical claims (pre-specified).")
    print(f"output -> {OUT}")


if __name__ == "__main__":
    main()
