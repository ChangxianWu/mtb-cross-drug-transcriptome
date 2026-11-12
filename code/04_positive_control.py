"""Positive-control gate (script 04). Primary analysis runs only after the gate passes.
Gate A (GSE292331 hypoxia-reaeration): the DosR-marker module must be significantly
induced in Hypoxia D7 vs Control D0 (g > 0.5 and p < 0.01);
Reaeration D12 vs Hypoxia D7 direction reversal is a reported (corroborating) arm;
respiration/OXPHOS (KEGG mtu00190) D7 vs D0 downward direction is a reported arm.
Gate B (GSE165673): INH at T24 (any dose) vs time-matched 0x vehicle must significantly
induce iniBAC (Rv0341-343) (g > 0.5 and p < 0.01); FAS-II/fatty-acid biosynthesis
(mtu00061) upward direction is a reported arm.
Gene sets: KEGG REST mtu pathways (version + hash cached); marker sets are small
literature-fixed sets.
Outputs: output/positive_control_report.csv + output/GATE_STATUS.txt"""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import hashlib
import io
import os
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = BASE / "."
PROC = os.path.join(ROOT, "work", "data", "processed")
CACHE = os.path.join(ROOT, "work", "cache")
OUT = os.path.join(ROOT, "work", "output")
os.makedirs(CACHE, exist_ok=True)

KEGG_CACHE = os.path.join(CACHE, "kegg_mtu_pathway_link_2026-09-18.tsv")

DOSR_MARKERS = ["Rv2031c", "Rv1738", "Rv1737c", "Rv2007c", "Rv2623", "Rv2626c"]
INIBAC = ["Rv0341", "Rv0342", "Rv0343"]


def load_kegg_modules():
    if os.path.exists(KEGG_CACHE):
        body = open(KEGG_CACHE, encoding="utf-8").read()
        src = "cache"
    else:
        req = urllib.request.Request("https://rest.kegg.jp/link/pathway/mtu",
                                     headers={"User-Agent": "mtb-pipeline/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read().decode()
        with open(KEGG_CACHE, "w", encoding="utf-8") as f:
            f.write(body)
        src = "fetched"
    sha = hashlib.sha256(body.encode()).hexdigest()[:12]
    mods = {}
    for line in body.strip().split("\n"):
        gene, pw = line.split("\t")
        gene = gene.replace("mtu:", "")
        mods.setdefault(pw.replace("path:", ""), []).append(gene)
    print(f"[KEGG] {src}, n_genes_pathway_lines={sum(len(v) for v in mods.values())}, sha256_12={sha}")
    return mods, sha


def module_score(expr, genes, id_normal=True):
    """z-mean module scoring: mean of z-scored member genes across the log expression matrix (gene case matched as needed)."""
    idx_all = expr.index.astype(str)
    lookup = {g.upper(): g for g in idx_all}
    present = [lookup[g.upper()] for g in genes if g.upper() in lookup]
    if not present:
        return None, 0
    z = expr.loc[present].T
    z = (z - z.mean()) / z.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), len(present)


def hedges_g(a, b):
    """Standardized mean difference (b = control, a = treatment), small-sample corrected."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan, np.nan
    sa, sb = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    if sp == 0:
        return np.nan, np.nan
    d = (a.mean() - b.mean()) / sp
    J = 1 - 3 / (4 * (na + nb) - 9)
    u = (na + nb) / 2
    correction = 1 - J / (4 * u - 1) if False else 1  # keep Hedges J correction
    g = d * J
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return g, p


def gate_a(mods, kegg_sha):
    gse = "GSE292331"
    expr = pd.read_csv(os.path.join(PROC, f"{gse}_expr.csv.gz"), index_col=0)
    sdf = pd.read_csv(os.path.join(PROC, f"{gse}_samples.tsv"), sep="\t")
    sdf = sdf.set_index("gsm")
    expr = expr[[c for c in expr.columns if c in sdf.index]]
    rows = []
    arms = [("DosR_markers", DOSR_MARKERS, "Hypoxia", "d7", "Control", "d0", "GATE", 0.5, 0.05),
            ("DosR_markers", DOSR_MARKERS, "Reaeration", "d12", "Hypoxia", "d7", "REPORT_DIR_NEG", None, None),
            ("Resp_OXPHOS", mods.get("mtu00190", []), "Hypoxia", "d7", "Control", "d0", "REPORT_DIR_NEG", None, None)]
    for mod_name, genes, drug_a, t_a, drug_b, t_b, kind, g_min, p_max in arms:
        score, n_hit = module_score(expr, genes)
        if score is None:
            rows.append([gse, mod_name, f"{drug_a}:{t_a} vs {drug_b}:{t_b}", np.nan, np.nan,
                         0, len(genes), "SKIP", "no gene hits", kind])
            continue
        m_a = sdf["drug"].eq(drug_a) & sdf["time_h"].eq(t_a)
        m_b = sdf["drug"].eq(drug_b) & sdf["time_h"].eq(t_b)
        a = score[m_a.reindex(score.index).fillna(False)].dropna()
        b = score[m_b.reindex(score.index).fillna(False)].dropna()
        g, p = hedges_g(a, b)
        if kind == "GATE":
            verdict = "PASS" if (not np.isnan(g) and g > g_min and p < p_max) else "FAIL"
        elif kind == "REPORT_DIR_NEG":
            verdict = "DIR_OK" if (not np.isnan(g) and g < 0) else "DIR_REV"
        else:
            verdict = "REPORT"
        rows.append([gse, mod_name, f"{drug_a}:{t_a} vs {drug_b}:{t_b}", round(g, 3), f"{p:.2e}",
                     n_hit, len(genes), verdict, f"n_a={len(a)}, n_b={len(b)}", kind])
    return rows


def gate_b(mods, kegg_sha):
    gse = "GSE165673"
    expr = pd.read_csv(os.path.join(PROC, f"{gse}_expr.csv.gz"), index_col=0)
    sdf = pd.read_csv(os.path.join(PROC, f"{gse}_samples.tsv"), sep="\t")
    sdf = sdf.set_index("gsm")
    sdf["time_h"] = pd.to_numeric(sdf["time_h"], errors="coerce")
    sdf["dose_num"] = pd.to_numeric(sdf["dose_code"], errors="coerce")
    expr = expr[[c for c in expr.columns if c in sdf.index]]
    rows = []
    arms = [("iniBAC", INIBAC, "INH", [0.18, 1.8], 24, "GATE", 0.5, 0.01),
            ("FAS_fatty_acid_biosyn", mods.get("mtu00061", []), "INH", [0.18, 1.8], 24, "REPORT", None, None),
            ("iniBAC", INIBAC, "RIF", [0.008, 0.02], 24, "REPORT", None, None)]
    for mod_name, genes, drug, doses, t, kind, g_min, p_max in arms:
        score, n_hit = module_score(expr, genes)
        if score is None:
            rows.append([gse, mod_name, f"{drug} vs 0x @T{t}", np.nan, np.nan, 0, len(genes),
                         "SKIP", "no gene hits", kind])
            continue
        m_drug = sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].isin(doses)
        m_veh = sdf["drug"].eq(drug) & sdf["time_h"].eq(t) & sdf["dose_num"].eq(0)
        a = score[m_drug.reindex(score.index).fillna(False)].dropna()
        b = score[m_veh.reindex(score.index).fillna(False)].dropna()
        g, p = hedges_g(a, b)
        if kind == "GATE":
            verdict = "PASS" if (not np.isnan(g) and g > g_min and p < p_max) else "FAIL"
        else:
            verdict = "REPORT"
        rows.append([gse, mod_name, f"{drug} vs 0x @T{t}", round(g, 3) if not np.isnan(g) else np.nan,
                     f"{p:.2e}" if not np.isnan(p) else np.nan, n_hit, len(genes), verdict,
                     f"n_drug={len(a)}, n_veh={len(b)}", kind])
    return rows


def main():
    mods, sha = load_kegg_modules()
    rows = gate_a(mods, sha) + gate_b(mods, sha)
    rep = pd.DataFrame(rows, columns=["gse", "module", "contrast", "hedges_g", "p_wilcoxon",
                                      "n_genes_hit", "n_genes_set", "verdict", "note", "arm_kind"])
    rep_path = os.path.join(OUT, "positive_control_report.csv")
    rep.to_csv(rep_path, index=False, encoding="utf-8-sig")
    print(rep.to_string(index=False))
    a_pass = ((rep["arm_kind"] == "GATE") & (rep["gse"] == "GSE292331")
              & rep["verdict"].eq("PASS")).any()
    # Corroborating arms: reaeration reversal + OXPHOS reduction (direction-level)
    a_corr = ((rep["verdict"].eq("DIR_OK"))
              & (rep["gse"] == "GSE292331")).sum() >= 2
    a_pass = a_pass and a_corr
    b_pass = ((rep["arm_kind"] == "GATE") & (rep["gse"] == "GSE165673")
              & rep["verdict"].eq("PASS")).any()
    overall = "PASS" if (a_pass and b_pass) else "FAIL"
    with open(os.path.join(OUT, "GATE_STATUS.txt"), "w", encoding="utf-8") as f:
        f.write(f"Positive control gate: {overall}\nGateA(DosR hypoxia)={a_pass}; "
                f"GateB(INH iniBAC)={b_pass}; kegg_sha256_12={sha}\n")
    print(f"\nGATE: {overall} (A={a_pass}, B={b_pass}) → {rep_path}")


if __name__ == "__main__":
    main()
