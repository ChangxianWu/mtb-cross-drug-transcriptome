"""Shared statistics library: module scoring (z-mean), Hedges g (small-sample
correction), DerSimonian-Laird random effects, BH-FDR, KEGG module panel loading
and dataset loaders."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import hashlib
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = BASE / "."
PROC = os.path.join(ROOT, "work", "data", "processed")
OUT = os.path.join(ROOT, "work", "output")
CACHE = os.path.join(ROOT, "work", "cache")
KEGG_CACHE = os.path.join(CACHE, "kegg_mtu_pathway_link_2026-09-18.tsv")

# literature-fixed marker sets (Tier-1)
DOSR_MARKERS = ["Rv2031c", "Rv1738", "Rv1737c", "Rv2007c", "Rv2623", "Rv2626c"]
INIBAC = ["Rv0341", "Rv0342", "Rv0343"]

# module panel (KEGG mtu pathways -> module families; one pathway per module for reproducibility)
MODULE_PANEL = {
    "respiration_oxphos": "mtu00190",
    "glycolysis": "mtu00010",
    "tca_cycle": "mtu00020",
    "pentose_phosphate": "mtu00030",
    "fatty_acid_biosyn": "mtu00061",
    "fatty_acid_degrad": "mtu00071",
    "peptidoglycan": "mtu00550",
    "trehalose_glycerolipid": "mtu00561",
    "ribosome": "mtu03010",
    "two_component_system": "mtu02020",
    "abc_transporters": "mtu02010",
    "sulfur_metabolism": "mtu00920",
    "nitrogen_metabolism": "mtu00910",
    "dosr_markers": None,      # fixed set
    "inibac": None,            # fixed set
}


def load_kegg_modules():
    if os.path.exists(KEGG_CACHE):
        body = open(KEGG_CACHE, encoding="utf-8").read()
    else:
        import urllib.request
        req = urllib.request.Request("https://rest.kegg.jp/link/pathway/mtu",
                                     headers={"User-Agent": "mtb-pipeline/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read().decode()
        os.makedirs(CACHE, exist_ok=True)
        with open(KEGG_CACHE, "w", encoding="utf-8") as f:
            f.write(body)
    sha = hashlib.sha256(body.encode()).hexdigest()[:12]
    pw2genes = {}
    for line in body.strip().split("\n"):
        gene, pw = line.split("\t")
        pw2genes.setdefault(pw.replace("path:", ""), []).append(gene.replace("mtu:", ""))
    modules = {}
    for mod, pw in MODULE_PANEL.items():
        if pw is None:
            modules[mod] = list({"dosr_markers": DOSR_MARKERS, "inibac": INIBAC}[mod])
        else:
            modules[mod] = sorted(set(pw2genes.get(pw, [])))
    return modules, sha


def load_dataset(gse):
    """Return (expr DataFrame genes x samples [log scale], samples DataFrame, log_scale: bool).
    expr is the raw parse output; TPM/normalized counts are unified to log2(x+1) here."""
    import glob
    expr_files = glob.glob(os.path.join(PROC, f"{gse}_*_expr.csv.gz")) or \
        glob.glob(os.path.join(PROC, f"{gse}_expr.csv.gz"))
    if not expr_files:
        return None, None, None
    expr = pd.read_csv(expr_files[0], index_col=0)
    samp_file = os.path.join(PROC, f"{gse}_samples.tsv")
    sdf = pd.read_csv(samp_file, sep="\t") if os.path.exists(samp_file) else None
    return expr, sdf, expr_files[0]


def log_transform(expr, kind):
    """kind: log2cpm (already log) or tpm/norm (needs log2(x+1))."""
    if kind == "log2cpm":
        return expr
    return np.log2(expr.clip(lower=0) + 1)


def dataset_expr_log(gse):
    """Loader with log-transform detection: re-log decided by the parse source_type."""
    expr, sdf, src_file = load_dataset(gse)
    if expr is None:
        return None, None
    if "log2cpm" in os.path.basename(src_file):
        pass  # script 03 already outputs log2CPM
    else:
        if sdf is not None and "source_type" in sdf.columns:
            st = str(sdf["source_type"].iloc[0])
        else:
            st = ""
        if "tpm" in st or "xlsx" in st or "tsv_counts_norm" in st or gse == "GSE72459":
            expr = log_transform(expr, "norm")
    return expr, sdf


def module_scores(expr, modules):
    """Return (samples x modules score matrix, module -> hit count dict); per-gene z within each dataset."""
    idx_upper = {str(g).upper(): str(g) for g in expr.index}
    z = expr.T
    z = (z - z.mean()) / z.std(ddof=0).replace(0, np.nan)
    scores, hits = {}, {}
    for mod, genes in modules.items():
        present = [idx_upper[g.upper()] for g in genes if g.upper() in idx_upper]
        hits[mod] = len(present)
        scores[mod] = z[present].mean(axis=1) if present else np.nan
    return pd.DataFrame(scores), hits


def hedges_g(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan, np.nan, (na, nb)
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    if sp == 0 or np.isnan(sp):
        return np.nan, np.nan, (na, nb)
    d = (a.mean() - b.mean()) / sp
    J = 1 - 3 / (4 * (na + nb) - 9)
    g = d * J
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return g, p, (na, nb)


def cliffs_delta(a, b):
    """Rank effect size [-1, 1]: bounded alternative when tiny variances blow up SMD."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (len(a) * len(b)) - 1.0


def dl_meta(ys, ses):
    """DerSimonian-Laird random effects; k=1 degenerates to the single-study effect."""
    ys = np.asarray(ys, float)
    ses = np.asarray(ses, float)
    k = len(ys)
    w = 1.0 / ses**2
    mu_fe = float((w * ys).sum() / w.sum())
    Q = float((w * (ys - mu_fe) ** 2).sum())
    df = k - 1
    C = w.sum() - (w**2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / C) if (df > 0 and C > 0) else 0.0
    w2 = 1.0 / (ses**2 + tau2)
    mu = float((w2 * ys).sum() / w2.sum())
    se_mu = float(np.sqrt(1.0 / w2.sum()))
    zval = mu / se_mu if se_mu > 0 else np.nan
    pval = 2 * stats.norm.sf(abs(zval)) if not np.isnan(zval) else np.nan
    i2 = max(0.0, (Q - df) / Q) if (Q > 0 and df > 0) else 0.0
    return dict(k=k, g_dl=mu, se=se_mu, Q=Q, I2=i2, p=pval)


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ---------------- contrast plans (shared by 06/07; a priori) ----------------

def plan_165673(sdf):
    plans = []
    sdf = sdf.copy()
    sdf["dose_num"] = pd.to_numeric(sdf["dose_code"], errors="coerce")
    sdf["time_h"] = pd.to_numeric(sdf["time_h"], errors="coerce")
    for drug in sorted(sdf["drug"].unique()):
        if drug in ("UNKNOWN",):
            continue
        for t in sorted(sdf["time_h"].dropna().unique()):
            drug_m = sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].gt(0)
            veh_m = sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].eq(0)
            if drug_m.sum() >= 2 and veh_m.sum() >= 2:
                plans.append((drug, t, veh_m, drug_m))
    return plans


def plan_292332(sdf):
    sdf = sdf.copy()
    sdf["time_h"] = pd.to_numeric(sdf["time_h"], errors="coerce")
    ctrl = sdf["drug"].eq("CTRL") & sdf["time_h"].eq(16)
    plans = []
    for drug in sorted(set(sdf.loc[sdf["role"].eq("drug"), "drug"])):
        drug_m = sdf["drug"].eq(drug) & sdf["time_h"].eq(16)
        if drug_m.sum() >= 2:
            plans.append((drug, 16, ctrl, drug_m))
    return plans


def plan_123294(sdf):
    sdf = sdf.copy()
    sdf["time_h"] = pd.to_numeric(sdf["time_h"], errors="coerce")
    plans = []
    for t in sorted(sdf["time_h"].dropna().unique()):
        drug_m = sdf["drug"].eq("DLM") & sdf["time_h"].eq(t)
        veh_m = sdf["drug"].eq("DMSO") & sdf["time_h"].eq(t)
        if drug_m.sum() >= 2 and veh_m.sum() >= 2:
            plans.append(("DLM", t, veh_m, drug_m))
    return plans


def plan_xlsx(sdf, drug):
    sdf = sdf.copy()
    drug_m = sdf["drug"].eq(drug)
    veh_m = sdf["drug"].eq("CTRL")
    if drug_m.sum() >= 2 and veh_m.sum() >= 2:
        t = pd.to_numeric(sdf["time_h"], errors="coerce").iloc[0]
        return [(drug, float(t), veh_m, drug_m)]
    return []
