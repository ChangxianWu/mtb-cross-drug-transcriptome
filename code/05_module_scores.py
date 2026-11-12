"""Module scoring (script 05)."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mtblib import PROC, OUT, load_kegg_modules, dataset_expr_log, module_scores

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

DATASETS = ["GSE165673", "GSE292332", "GSE292331",
            "GSE123294_H37Rv", "GSE123294_041195", "GSE288598", "GSE287913"]


def main():
    modules, sha = load_kegg_modules()
    print(f"[KEGG] module panel: {len(modules)} modules, kegg_sha256_12={sha}")
    hit_rows = []
    for gse in DATASETS:
        expr, sdf = dataset_expr_log(gse)
        if expr is None:
            print(f"[SKIP] {gse}: no processed expression matrix")
            continue
        scores, hits = module_scores(expr, modules)
        scores.to_csv(os.path.join(PROC, f"{gse}_modscores.csv"), encoding="utf-8-sig")
        hit_rows.append({"gse": gse, "n_samples": expr.shape[1], **hits})
        flag = "" if all(v > 0 for v in hits.values()) else "  <- note: zero-hit modules present"
        print(f"[OK] {gse}: {expr.shape[1]} samples x {scores.shape[1]} modules{flag}")
    hits_df = pd.DataFrame(hit_rows).set_index("gse").T
    hits_df.to_csv(os.path.join(OUT, "module_hits.csv"), encoding="utf-8-sig")
    print("\nModule gene hit counts per dataset:")
    print(hits_df.to_string())
    print(f"\nhit table -> {os.path.join(OUT, 'module_hits.csv')}")


if __name__ == "__main__":
    main()
