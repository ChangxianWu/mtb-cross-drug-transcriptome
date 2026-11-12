"""Unified parsing of downloaded GEO series (script 03).
Outputs: data/processed/{GSE}_expr.csv.gz (genes x samples) or {GSE}_lfc.csv (Mode B);
         data/processed/{GSE}_samples.tsv (standardized sample table);
         output/PARSE_LOG.tsv (per-step status; deferred/excluded recorded explicitly).
Sample-table columns: gse, gsm, drug, dose_raw, dose_code, time_h, batch, replicate,
role, source_type, notes.
Rules: no merging outside defined modules; every mapping assumption goes into notes;
a parse failure does not abort the remaining series."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import gzip
import os
import re
import sys
import tarfile

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = BASE / "."
RAW = os.path.join(ROOT, "work", "data", "raw")
PROC = os.path.join(ROOT, "work", "data", "processed")
OUT = os.path.join(ROOT, "work", "output")
LOG = os.path.join(OUT, "PARSE_LOG.tsv")
os.makedirs(PROC, exist_ok=True)

log_rows = []


def log(gse, step, status, note=""):
    log_rows.append([gse, step, status, note])
    print(f"[{status}] {gse} {step} :: {note}")


def save(gse, tag, df, kind="expr"):
    p = os.path.join(PROC, f"{gse}_{tag}_{kind}.csv.gz" if kind == "expr" else f"{gse}_{tag}_{kind}.csv")
    df.to_csv(p, index=False, compression="gzip" if kind == "expr" else None)
    return p


def save_samples(gse, df):
    p = os.path.join(PROC, f"{gse}_samples.tsv")
    df.to_csv(p, sep="\t", index=False, encoding="utf-8-sig")
    return p


def read_series_matrix(path):
    """GSEMatrix parsing: returns meta(rows), expr DataFrame (when a value table exists), platform_id."""
    meta_rows, expr_lines, in_expr, plat_id = [], [], False, ""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                plat_id = line.rstrip("\n").split("\t")[1]
            if line.startswith("!series_matrix_table_begin"):
                in_expr = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_expr = False
                continue
            if in_expr:
                expr_lines.append(line.rstrip("\n"))
            elif line.startswith("!"):
                meta_rows.append(line.rstrip("\n").split("\t"))
    meta = {}
    for row in meta_rows:
        key = row[0]
        meta.setdefault(key, []).append(row[1:])
    expr = None
    if expr_lines:
        header = expr_lines[0].lstrip("!").split("\t")
        data = [l.split("\t") for l in expr_lines[1:]]
        expr = pd.DataFrame(data[0:], columns=header).set_index(header[0])
        expr = expr.apply(pd.to_numeric, errors="coerce")
    return meta, expr, plat_id


def unquote(x):
    return x.strip().strip('"') if isinstance(x, str) else x


# ---------------- GSE72459: series matrix carries expression values ----------------

def parse_72459():
    gse = "GSE72459"
    p = os.path.join(RAW, gse, f"{gse}_series_matrix.txt.gz")
    meta, expr, plat = read_series_matrix(p)
    if expr is None:
        log(gse, "matrix", "deferred", "no expression table")
        return
    accs = [unquote(x) for x in meta["!Sample_geo_accession"][0]]
    n_char = len(meta.get("!Sample_characteristics_ch1", []))
    rows = []
    for i, gsm in enumerate(accs):
        # this series uses one characteristics row per sample (three segments)
        line_vals = meta["!Sample_characteristics_ch1"][i] if i < n_char else []
        dose = time_h = None
        for v in line_vals:
            v = unquote(v)
            m = re.search(r"dose:\s*([\d.]+)\s*μ?M", v)
            if m:
                dose = float(m.group(1))
            m2 = re.search(r"for\s+(\d+)\s*hr", v)
            if m2:
                time_h = int(m2.group(1))
        rows.append(dict(gse=gse, gsm=gsm, drug="BDQ" if (dose or 0) > 0 else "CTRL",
                         dose_raw=f"{dose}uM" if dose is not None else "0uM",
                         dose_code="drug" if (dose or 0) > 0 else "0x",
                         time_h=time_h, batch="", replicate="", role="drug" if (dose or 0) > 0 else "vehicle",
                         source_type="matrix_values",
                         notes=f"plat={plat}; characteristics per-sample-row"))
    sdf = pd.DataFrame(rows)
    expr.columns = accs
    save(gse, "bdq", expr)
    save_samples(gse, sdf)
    log(gse, "matrix", "ok", f"{expr.shape[0]} genes x {expr.shape[1]} samples; plat={plat}")


# ---------------- GSE165673 SuperSeries (207 per-sample TPM files) ----------------

DRUG_ALIAS = {"PA824": "PTN", "BDQ": "BDQ", "INH": "INH", "RIF": "RIF", "EMB": "EMB",
              "PZA": "PZA", "LZD": "LZD", "DLM": "DLM", "POA": "POA", "MXF": "MXF",
              "CFZ": "CFZ", "AMK": "AMK", "LVF": "LVF", "CPF": "CPF", "RPE": "RPE"}


def parse_165673():
    gse = "GSE165673"
    tar_p = os.path.join(RAW, gse, f"{gse}_RAW.tar")
    # metadata: series matrices of two platforms
    gsm_meta = {}
    for plat_file in [f"{gse}-GPL22688_series_matrix.txt.gz", f"{gse}-GPL23592_series_matrix.txt.gz"]:
        meta, _, _ = read_series_matrix(os.path.join(RAW, gse, plat_file))
        accs = [unquote(x) for x in meta["!Sample_geo_accession"][0]]
        titles = [unquote(x) for x in meta["!Sample_title"][0]]
        char_rows = [[unquote(v) for v in r] for r in meta.get("!Sample_characteristics_ch1", [])]
        for i, gsm in enumerate(accs):
            drug_txt = char_rows[0][i] if char_rows and i < len(char_rows[0]) else ""
            time_txt = char_rows[2][i] if len(char_rows) > 2 and i < len(char_rows[2]) else ""
            t = titles[i] if i < len(titles) else ""
            m = re.match(r"^(?:[A-Za-z0-9]+_)?T(\d+)_([A-Za-z0-9]+)_([\d.]+)x?_([A-Za-z])$", t)
            if m:
                time_h, drug_tok, dose, rep = int(m.group(1)), m.group(2), m.group(3), m.group(4)
            else:
                time_h, drug_tok, dose, rep = None, "", "", ""
                log(gse, "title_regex", "needs_review", f"unmatched title: {t}")
            gsm_meta[gsm] = dict(title=t, drug_tok=drug_tok, dose_code=dose, time_h=time_h,
                                 replicate=rep, drug_txt=drug_txt, time_txt=time_txt)
    log(gse, "meta", "ok", f"{len(gsm_meta)} GSM metadata entries (two platform matrices)")
    # count matrices: read counts/TPM per tar member (READS_M raw counts preferred; unified to log2CPM; missing columns skipped and logged)
    cols, name_map, kinds = {}, {}, {}
    with tarfile.open(tar_p) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            gsm = m.name.split("_")[0]
            df = pd.read_csv(tf.extractfile(m), sep="\t", compression="gzip", dtype=str)
            gene_col = "GENE_ID" if "GENE_ID" in df.columns else df.columns[0]
            name_col = "GENE_NAME" if "GENE_NAME" in df.columns else None
            val_col, kind = next(
                ((c, k) for c, k in [("READS_M", "counts"), ("TPM_M", "tpm"), ("TPM", "tpm")]
                 if c in df.columns), (None, None))
            if val_col is None:
                log(gse, "member_cols", "needs_review", f"{m.name}: columns={list(df.columns)[:6]}")
                continue
            s = pd.Series(pd.to_numeric(df[val_col], errors="coerce").values,
                          index=df[gene_col].astype(str), name=gsm)
            cols[gsm] = s
            kinds[gsm] = kind
            if name_col:
                nm = df.set_index(df[gene_col].astype(str))[name_col]
                for k, v in nm.items():
                    if isinstance(v, str) and v:
                        name_map.setdefault(k, v)
    expr = pd.DataFrame(cols)
    dup = expr.columns[expr.columns.duplicated()].tolist()
    if dup:
        log(gse, "dup_gsm", "needs_review", f"duplicate GSM columns: {dup}")
        expr = expr.loc[:, ~expr.columns.duplicated()]
    kind_counts = {"counts": 0, "tpm": 0}
    for c in expr.columns:
        kind_counts[kinds.get(c, "?")] = kind_counts.get(kinds.get(c, "?"), 0) + 1
    if kind_counts.get("counts", 0) and kind_counts.get("tpm", 0):
        log(gse, "mixed_units", "needs_review",
            f"mixed units counts={kind_counts['counts']} tpm={kind_counts['tpm']}; unified to log2CPM via counts")
    counts = expr.loc[:, [c for c in expr.columns if kinds.get(c) == "counts"]]
    if counts.shape[1] > 0:
        cpm = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
        expr_final = np.log2(cpm + 1)
        src_type = "counts_log2cpm"
    else:
        expr_final = expr
        src_type = "tpm"
    expr_final.to_csv(os.path.join(PROC, f"{gse}_expr.csv.gz"), compression="gzip")
    pd.Series(name_map, name="gene_name").to_csv(
        os.path.join(PROC, f"{gse}_gene_names.csv"), encoding="utf-8-sig")
    # 0x vehicle samples without a drug token: inherit from nearest GSM (within-block, explicit assumption)
    gsm_sorted = sorted(gsm_meta.keys(), key=lambda g: int(g[3:]))
    for gsm in gsm_sorted:
        d = gsm_meta[gsm]
        if d["drug_tok"] == "" and d["dose_code"] == "0":
            i = gsm_sorted.index(gsm)
            inherited = None
            for j in list(range(i - 1, -1, -1)) + list(range(i + 1, len(gsm_sorted))):
                dd = gsm_meta[gsm_sorted[j]]
                if dd["drug_tok"] and dd["dose_code"] == "0":
                    inherited = dd["drug_tok"]
                    break
            if inherited:
                d["drug_tok"] = inherited
                d["notes_inherit"] = f"vehicle samples lack a drug token; inherit from nearest GSM {inherited}"
                log(gse, "vehicle_inherit", "ok", f"{gsm} ← {inherited}（nearest-neighbour inheritance assumption）")
    rows = []
    for gsm, d in gsm_meta.items():
        if gsm not in expr_final.columns:
            log(gse, "missing_member", "needs_review", f"matrix missing {gsm}")
            continue
        drug = DRUG_ALIAS.get(d["drug_tok"], d["drug_tok"])
        try:
            is_vehicle = float(d["dose_code"]) == 0
        except ValueError:
            is_vehicle = False
        role = "vehicle" if is_vehicle else "drug"
        if not d["drug_tok"]:
            drug, role = "UNKNOWN", "needs_review"
        rows.append(dict(gse=gse, gsm=gsm, drug=drug, dose_raw=d["drug_txt"],
                         dose_code=d["dose_code"], time_h=d["time_h"], batch="",
                         replicate=d["replicate"], role=role, source_type=src_type,
                         notes=f"title={d['title']}; {d.get('notes_inherit','')}"))
    sdf = pd.DataFrame(rows)
    save_samples(gse, sdf)
    log(gse, "matrix", "ok", f"{expr_final.shape[0]} genes x {expr_final.shape[1]} samples ({src_type})")


# ---------------- GSE292332 / GSE292331: processed TPM + series-matrix metadata ----------------

def parse_tpm_series(gse, fname, mode="drug"):
    """mode="drug": 292332 (Treatment/Time/Batch metadata inside files; labels are sample IDs).
    mode='hypoxia': 292331 (Strain column only; treatment derived from _D{day}_ labels, cross-checked against matrix features)."""
    p = os.path.join(RAW, gse, fname)
    raw = pd.read_csv(p, sep="\t", compression="gzip", index_col=0)
    meta, _, _ = read_series_matrix(os.path.join(RAW, gse, f"{gse}_series_matrix.txt.gz"))
    accs = [unquote(x) for x in meta["!Sample_geo_accession"][0]]
    titles = [unquote(x) for x in meta["!Sample_title"][0]]
    chars = [[unquote(v) for v in r] for r in meta.get("!Sample_characteristics_ch1", [])]

    def char_at(i, key):
        for row in chars:
            if row and row[0].lower().startswith(key):
                return row[i] if i < len(row) else ""
        return ""

    rv_like_cols = [c for c in raw.columns if re.match(r"^(Rv|Rvnc|MTS)", str(c), re.I)]
    meta_cols = [c for c in raw.columns if not re.match(r"^(Rv|Rvnc|MTS)", str(c), re.I)]
    if len(rv_like_cols) > 100:  # transpose: rows are samples, labels are sample IDs (no GSM mapping; recorded)
        meta_frame = raw[meta_cols].copy()
        expr = raw[rv_like_cols].T
        labels = list(raw.index)
        log(gse, "orientation", "ok",
            f"transposed layout: {len(labels)} samples x {len(rv_like_cols)} genes; labels are IDs; metadata columns={meta_cols}")
    else:
        raise ValueError(f"non-transposed layout not supported: columns={list(raw.columns)[:5]}")

    if mode == "drug":
        rows = []
        for lab in labels:
            treat = str(meta_frame.loc[lab].get("Treatment", ""))
            tme = str(meta_frame.loc[lab].get("Time", ""))
            batch = str(meta_frame.loc[lab].get("Batch", ""))
            drug, dose_raw, role, veh = "CTRL", "", "control", ""
            m = re.match(r"^([A-Za-z0-9]+)\+([A-Za-z0-9]+),\s*(.+)$", treat)
            m2 = re.match(r"^([A-Za-z0-9]+),\s*([\d.]+)\s*(\S*)$", treat)
            if m and m.group(1).upper() != "UNTREATED":
                drug = DRUG_ALIAS.get(m.group(1).upper(), m.group(1).upper())
                veh = m.group(2)
                dose_raw = m.group(3)
                role = "vehicle" if dose_raw.strip("0 ") == "" else "drug"
            elif m2:
                drug = DRUG_ALIAS.get(m2.group(1).upper(), m2.group(1).upper())
                dose_raw = m2.group(2) + (m2.group(3) or "")
                role = "vehicle" if float(m2.group(2)) == 0 else "drug"
            elif treat.lower().startswith("untreated"):
                drug, role = "CTRL", "control"
            rows.append(dict(gse=gse, gsm=lab, drug=drug, dose_raw=dose_raw, dose_code="",
                             time_h=str(tme).replace("hours", "").strip(), batch=batch, replicate="",
                             role=role, source_type="processed_tpm",
                             notes=f"treatment={treat}; vehicle={veh}; labels are IDs (no GSM mapping)"))
        sdf = pd.DataFrame(rows)
    else:  # hypoxia: treatment derived from _D{day}_ labels
        rows = []
        mismatch = 0
        for i, lab in enumerate(labels):
            m = re.search(r"_D(\d+)_", str(lab))
            day = int(m.group(1)) if m else None
            if day is None:
                log(gse, "label_day", "needs_review", f"label lacks D{n}: {lab}")
                continue
            treat = "Control" if day == 0 else ("Hypoxia" if day <= 7 else "Reaeration")
            # positional cross-check against series matrix (count consistency only)
            c_treat = char_at(i, "treatment")
            if c_treat and c_treat.split(":")[-1].strip().lower() != treat.lower():
                mismatch += 1
            rows.append(dict(gse=gse, gsm=lab, drug=treat, dose_raw="", dose_code="",
                             time_h=f"d{day}", batch=str(meta_frame.loc[lab].get("Strain", "")),
                             replicate="", role="stress" if treat != "Control" else "control",
                             source_type="processed_tpm",
                             notes=f"label-derived day={day}; cross-checked against matrix features"))
        sdf = pd.DataFrame(rows)
        log(gse, "crosscheck", "ok" if mismatch == 0 else "needs_review",
            f"label-derived vs matrix-feature mismatches={mismatch}/{len(labels)}")
    expr.to_csv(os.path.join(PROC, f"{gse}_expr.csv.gz"), compression="gzip")
    save_samples(gse, sdf)
    log(gse, "matrix", "ok", f"{expr.shape[0]} genes x {expr.shape[1]} samples")


# ---------------- GSE199263: DESeq2 LFC table (Mode B) ----------------

def parse_199263():
    gse = "GSE199263"
    stress_map = {"stationary": "stationary", "oxy": "microaerophilic", "PH5": "pH5", "PBS": "PBS"}
    frames = []
    for f, stress in stress_map.items():
        p = os.path.join(RAW, gse, f"GSE199263_{f}_exp_LFC.deseq.gene.txt.gz")
        df = pd.read_csv(p, sep="\t", compression="gzip")
        df["stress"] = stress
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(os.path.join(PROC, f"{gse}_lfc.csv"), index=False)
    sdf = pd.DataFrame([dict(gse=gse, gsm="", drug=s, dose_raw="vs standard medium?", dose_code="",
                             time_h="", batch="", replicate="", role="lfc_table", source_type="deseq_lfc",
                         notes="Mode B: gene-level LFC aggregation, no per-sample values") for s in stress_map.values()])
    save_samples(gse, sdf)
    log(gse, "lfc", "ok", f"{out.shape[0]} rows x 4 stresses")


# ---------------- GSE288598 / GSE287913: xlsx counts (positional mapping assumption) ----------------

def parse_xlsx_counts(gse, fname, drug, veh_name):
    p = os.path.join(RAW, gse, fname)
    df = pd.read_excel(p, sheet_name=0)
    gene_col = df.columns[0]
    expr = df.set_index(gene_col).apply(pd.to_numeric, errors="coerce")
    meta, _, _ = read_series_matrix(os.path.join(RAW, gse, f"{gse}_series_matrix.txt.gz"))
    accs = [unquote(x) for x in meta["!Sample_geo_accession"][0]]
    chars = [[unquote(v) for v in r] for r in meta.get("!Sample_characteristics_ch1", [])]

    def char_for(i, key):
        for row in chars:
            if row and row[0].lower().startswith(key):
                return row[i] if i < len(row) else ""
        return ""

    rows = []
    if expr.shape[1] == len(accs):
        expr.columns = accs  # positional mapping assumption (recorded)
        for i, gsm in enumerate(accs):
            treat = char_for(i, "treatment")
            t_low = treat.lower()
            is_drug = ("exposed" in t_low) and ("unexposed" not in t_low)
            rows.append(dict(gse=gse, gsm=gsm,
                             drug=drug if is_drug else "CTRL",
                             dose_raw=char_for(i, "time"), dose_code="MIC/4" if is_drug else "",
                             time_h="24", batch="", replicate="",
                             role="drug" if is_drug else "vehicle",
                             source_type="xlsx_counts",
                             notes="xlsx columns <-> GSM aligned by position (assumption, recorded)"))
        expr.to_csv(os.path.join(PROC, f"{gse}_expr.csv.gz"), compression="gzip")
        save_samples(gse, pd.DataFrame(rows))
        log(gse, "matrix", "ok", f"{expr.shape[0]} genes x {expr.shape[1]} samples(positional mapping assumption)")
    else:
        log(gse, "colmap", "deferred", f"column count mismatch expr={expr.shape[1]} vs GSM={len(accs)}")


# ---------------- GSE123294: TSV text 'complete' table ----------------

def parse_123294():
    gse = "GSE123294"
    # gene-name -> Rv bridge table from the GSE165673 GENE_ID/GENE_NAME mapping
    bridge_p = os.path.join(PROC, "GSE165673_gene_names.csv")
    bridge = {}
    if os.path.exists(bridge_p):
        b = pd.read_csv(bridge_p, index_col=0).dropna()
        bridge = dict(zip(b["gene_name"].astype(str), b.index.astype(str)))
    for f, strain in [("GSE123294_H37RV-del.complete.xls.gz", "H37Rv"),
                      ("GSE123294_041195-del.complete.xls.gz", "041195")]:
        df = pd.read_csv(os.path.join(RAW, gse, f), sep="\t", compression="gzip")
        name_col = "Name" if "Name" in df.columns else df.columns[1]
        id_col = df.columns[0]
        # sample-column pattern: {norm.}{04-}{0|6|24}h{Dela|DMSO}-{rep}; prefer norm. columns; exclude DESeq stat columns
        sample_re = re.compile(r"^(norm\.)?H?(?:04-)?(\d+)h(Dela|DMSO)-([A-Za-z])$")
        cands = [c for c in df.columns if sample_re.match(str(c))]
        norm_cands = [c for c in cands if str(c).startswith("norm.")]
        sample_cols = norm_cands if norm_cands else cands
        if not sample_cols:
            log(gse, f"cols_{strain}", "deferred", "no sample columns matched")
            continue
        log(gse, f"cols_{strain}", "ok",
            f"sample columns {len(sample_cols)} (norm-first={bool(norm_cands)}), DESeq stat columns excluded")
        expr = df.set_index(name_col)[sample_cols].apply(pd.to_numeric, errors="coerce")
        rows, keep = [], []
        for c in sample_cols:
            m = sample_re.match(str(c))
            time_h, arm, rep = int(m.group(2)), m.group(3), m.group(4)
            keep.append(c)
            rows.append(dict(gse=gse, gsm=str(c),
                             drug="DLM" if arm == "Dela" else "DMSO",
                             dose_raw="", dose_code="",
                             time_h=time_h, batch=strain, replicate=rep,
                             role="baseline" if time_h == 0 else ("drug" if arm == "Dela" else "vehicle"),
                             source_type="tsv_counts_norm" if norm_cands else "tsv_counts",
                             notes=f"arm={arm}; strain={strain}(gsm = file column name; strain recorded in batch)"))
        # name -> Rv mapping (unmapped genes keep original names; mapping counts logged)
        new_idx = [bridge.get(str(n), str(n)) for n in expr.index]
        n_mapped = sum(1 for a, b in zip(expr.index, new_idx) if a != b)
        log(gse, "name_bridge", "ok", f"{strain}: {n_mapped}/{len(new_idx)} gene names mapped to Rv (unmapped keep original names)")
        expr.index = new_idx
        expr = expr[keep]
        expr = expr[~expr.index.duplicated(keep="first")]
        expr.to_csv(os.path.join(PROC, f"{gse}_{strain}_expr.csv.gz"), compression="gzip")
        save_samples(f"{gse}_{strain}", pd.DataFrame(rows))
        log(gse, f"matrix_{strain}", "ok", f"{expr.shape[0]} genes x {expr.shape[1]} samples")


def main():
    for gse, fn in [("GSE72459", parse_72459),
                    ("GSE165673", parse_165673),
                    ("GSE292332", lambda: parse_tpm_series("GSE292332",
                        "GSE292332_FrontiersInTb_INDIGO1Sherman_GEO_processed_tpm.txt.gz", mode="drug")),
                    ("GSE292331", lambda: parse_tpm_series("GSE292331",
                        "GSE292331_FrontiersInTb_Hypoxia_GEO_processed_tpm.txt.gz", mode="hypoxia")),
                    ("GSE199263", parse_199263),
                    ("GSE288598", lambda: parse_xlsx_counts("GSE288598",
                        "GSE288598_Expression_Estimation_Linezolid.xlsx", "LZD", "unexposed")),
                    ("GSE287913", lambda: parse_xlsx_counts("GSE287913",
                        "GSE287913_Expression_Estimation_Delamanid.xlsx", "DLM", "control")),
                    ("GSE123294", parse_123294)]:
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            log(gse, "exception", "FAIL", f"{type(e).__name__}: {str(e)[:200]}")
    # explicit exclusions/deferrals (deviation log)
    log("GSE154627", "modality", "excluded", "TnSeq fitness data (not expression) -> logged as R2 candidate")
    log("GSE319987", "id_mapping", "deferred", "Erdman-strain ERDMAN_RS locus tags; needs Rv mapping before use")
    log("GSE119585", "no_suppl", "deferred", "SuperSeries no supplementary files (HTTP 404); use SubSeries/SRA route")
    import csv as _csv
    with open(LOG, "w", newline="", encoding="utf-8-sig") as f:
        w = _csv.writer(f, delimiter="\t")
        w.writerow(["gse", "step", "status", "note"])
        w.writerows(log_rows)
    print(f"\nPARSE_LOG → {LOG}")


if __name__ == "__main__":
    main()
