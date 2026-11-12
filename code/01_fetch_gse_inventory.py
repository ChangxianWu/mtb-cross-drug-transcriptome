"""GEO inventory for anti-tuberculosis drug-perturbation transcriptome series:
metadata plus FTP file listings for candidate GSE accessions."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

UA = "horzontal-pipeline mtb-pipeline/1.0 (research use)"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
OUT_TSV = BASE / "." / "work" / "output" / "GEO_INVENTORY.tsv"
META_DIR = BASE / "." / "work" / "data" / "gse_meta"

GSES = [
    "GSE72459", "GSE166546", "GSE165671", "GSE166271", "GSE288598",
    "GSE287913", "GSE165673", "GSE154627", "GSE319987", "GSE292332",
    "GSE119585", "GSE123294", "GSE292331", "GSE199263",
]

# supplementary-file classification: counts/expr/tpm/fpkm/normalized -> download; fasta/fastq/raw/sra -> skip
KEEP_RE = re.compile(r"counts|expr|tpm|fpkm|norm|quant|feature|read_?gene", re.I)
SKIP_RE = re.compile(r"fasta|fastq|\.sra|raw_?read|_md5\.txt$", re.I)
MAX_BYTES = 300 * 1024 * 1024


def get(url, tries=4, timeout=45, head=False):
    last = None
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA},
                                         method="HEAD" if head else "GET")
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if head:
                    return r.headers.get("Content-Length"), r.status
                return r.read().decode("utf-8", errors="replace"), r.status
        except Exception as e:  # noqa: BLE001
            last = e
            if t < tries - 1:
                time.sleep(2.0 * (t + 1))
    print(f"  !! GET failed ({tries} tries): {url} :: {last}")
    return None, None


def gse_prefix(gse):
    n = int(gse[3:])
    return f"GSE{str(n)[:-3]}nnn"


def fetch_uid(gse):
    u = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "term": f"{gse}[ACCN]", "retmode": "json"})
    body, _ = get(u)
    if body is None:
        return None
    ids = json.loads(body)["esearchresult"].get("idlist", [])
    return ids[0] if ids else None


def fetch_summary(uid):
    u = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "id": uid, "retmode": "json"})
    body, _ = get(u)
    if body is None:
        return None
    d = json.loads(body).get("result", {})
    r = d.get(uid, {}) if uid in d else (d.get("uids", [None])[0] and d[d["uids"][0]])
    return r or {}


def ftp_list(url):
    body, _ = get(url)
    if body is None:
        return []
    # keep relative names only: exclude directories (trailing /), absolute URLs, footer links
    names = []
    for f in re.findall(r'href="([^"?]+)"', body):
        f = f.rstrip("/")
        if not f or f.startswith("/") or "://" in f:
            continue
        names.append(f)
    return names


def head_size(url):
    cl, _ = get(url, head=True)
    try:
        return int(cl) if cl else None
    except Exception:
        return None


def classify(fname, size):
    if fname.endswith("series_matrix.txt.gz"):
        return "matrix"
    if SKIP_RE.search(fname):
        return "skip_raw"
    if KEEP_RE.search(fname):
        return "counts_expr" if (size is not None and size <= MAX_BYTES) else "skip_large_or_unknown"
    # other files (xlsx/csv/txt/gz, unknown content): recorded for manual review
    return "review" if (size is not None and size <= MAX_BYTES) else "skip_large_or_unknown"


def main():
    os.makedirs(META_DIR, exist_ok=True)
    rows = []
    for gse in GSES:
        print(f"--- {gse}")
        uid = fetch_uid(gse)
        meta = fetch_summary(uid) if uid else {}
        title = (meta.get("title") or "")[:120]
        platform = meta.get("gpl") or ""
        n_samples = meta.get("n_samples") or ""
        pmid = ""
        pl = meta.get("pubmedids") or meta.get("pubmedid") or []
        if isinstance(pl, list) and pl:
            pmid = pl[0]
        summary = (meta.get("summary") or "")[:300]

        base = f"{FTP}/{gse_prefix(gse)}/{gse}"
        matrix_url = f"{base}/matrix/"
        suppl_url = f"{base}/suppl/"
        matrix_files = [f for f in ftp_list(matrix_url) if f.endswith("series_matrix.txt.gz")]
        suppl_files = ftp_list(suppl_url)

        for f in matrix_files:
            u = f"{matrix_url}/{f}"
            size = head_size(u)
            rows.append([gse, uid, title, platform, n_samples, pmid, "matrix", f, u,
                         size, classify(f, size), summary])
        for f in suppl_files:
            u = f"{suppl_url}/{f}"
            size = head_size(u)
            rows.append([gse, uid, title, platform, n_samples, pmid, "suppl", f, u,
                         size, classify(f, size), summary])
        if not matrix_files and not suppl_files:
            rows.append([gse, uid, title, platform, n_samples, pmid, "none", "", "", None,
                         "listing_failed_or_empty", summary])
        with open(f"{META_DIR}\\{gse}.json", "w", encoding="utf-8") as fj:
            json.dump({"uid": uid, "meta": meta, "matrix_files": matrix_files,
                       "suppl_files": suppl_files}, fj, ensure_ascii=False, indent=1)
        time.sleep(0.5)

    with open(OUT_TSV, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gse", "uid", "title", "platform", "n_samples", "pmid", "kind",
                    "file", "url", "bytes", "fetch_class", "summary_head"])
        w.writerows(rows)
    print(f"\ninventory done: {len(rows)} file rows -> {OUT_TSV}")


if __name__ == "__main__":
    main()
