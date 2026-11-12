"""Selective GEO download with byte-length checks, per-file manifest and SHA-256."""
import os as _os
BASE = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import csv
import hashlib
import os
import sys
import time
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

UA = "horzontal-pipeline mtb-pipeline/1.0 (research use)"
INV = BASE / "." / "work" / "output" / "GEO_INVENTORY.tsv"
RAW = BASE / "." / "work" / "data" / "raw"
MANIFEST = BASE / "." / "work" / "data" / "MANIFEST.tsv"
VALID = BASE / "." / "work" / "data" / "VALIDATION.tsv"

# deterministic download plan: GSE -> filename list (matrix files always included)
PLAN = {
    "GSE72459": ["GSE72459_series_matrix.txt.gz"],
    "GSE166546": ["GSE166546_series_matrix.txt.gz", "GSE166546_RAW.tar"],
    "GSE165671": ["GSE165671_series_matrix.txt.gz", "GSE165671_RAW.tar"],
    "GSE166271": ["GSE166271_series_matrix.txt.gz", "GSE166271_RAW.tar"],
    "GSE288598": ["GSE288598_series_matrix.txt.gz", "GSE288598_Expression_Estimation_Linezolid.xlsx"],
    "GSE287913": ["GSE287913_series_matrix.txt.gz", "GSE287913_Expression_Estimation_Delamanid.xlsx"],
    "GSE165673": ["GSE165673-GPL22688_series_matrix.txt.gz", "GSE165673-GPL23592_series_matrix.txt.gz",
                  "GSE165673_RAW.tar"],
    "GSE154627": ["GSE154627_series_matrix.txt.gz", "GSE154627_RAW.tar"],
    "GSE319987": ["GSE319987_series_matrix.txt.gz", "GSE319987_RAW.tar"],
    "GSE292332": ["GSE292332_series_matrix.txt.gz",
                  "GSE292332_FrontiersInTb_INDIGO1Sherman_GEO_processed_tpm.txt.gz",
                  "GSE292332_FrontiersInTb_INDIGO1Sherman_GEO_processed_reads.txt.gz"],
    "GSE119585": ["GSE119585-GPL22688_series_matrix.txt.gz", "GSE119585-GPL23592_series_matrix.txt.gz"],
    "GSE123294": ["GSE123294_series_matrix.txt.gz", "GSE123294_041195-del.complete.xls.gz",
                  "GSE123294_H37RV-del.complete.xls.gz"],
    "GSE292331": ["GSE292331_series_matrix.txt.gz",
                  "GSE292331_FrontiersInTb_Hypoxia_GEO_processed_tpm.txt.gz",
                  "GSE292331_FrontiersInTb_Hypoxia_GEO_processed_reads.txt.gz"],
    "GSE199263": ["GSE199263_series_matrix.txt.gz", "GSE199263_PBS_exp_LFC.deseq.gene.txt.gz",
                  "GSE199263_PH5_exp_LFC.deseq.gene.txt.gz", "GSE199263_oxy_exp_LFC.deseq.gene.txt.gz",
                  "GSE199263_stationary_exp_LFC.deseq.gene.txt.gz"],
}


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_inventory():
    idx = {}
    with open(INV, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            idx[(row["gse"], row["file"])] = row
    return idx


def download(url, dest, expected_bytes):
    part = dest + ".part"
    last = None
    for t in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                cl = r.headers.get("Content-Length")
                if cl and expected_bytes and int(cl) != int(expected_bytes):
                    print(f"  !! Content-Length mismatch: {cl} vs {expected_bytes}")
                with open(part, "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
            actual = os.path.getsize(part)
            if expected_bytes and actual != int(expected_bytes):
                raise IOError(f"byte-count mismatch {actual} != {expected_bytes}")
            os.replace(part, dest)
            return actual
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2.0 * (t + 1))
    raise IOError(f"download failed (4 tries): {url} :: {last}")


def main():
    idx = load_inventory()
    man_rows, val_rows = [], []
    n_ok = n_fail = 0
    for gse, files in PLAN.items():
        os.makedirs(os.path.join(RAW, gse), exist_ok=True)
        for fname in files:
            row = idx.get((gse, fname))
            if row is None:
                print(f"[MISS] {gse} {fname} not in inventory")
                n_fail += 1
                continue
            url = row["url"]
            dest = os.path.join(RAW, gse, fname)
            if os.path.exists(dest) and os.path.getsize(dest) == int(row["bytes"] or 0):
                actual = os.path.getsize(dest)
                print(f"[SKIP-EXISTS] {gse}/{fname} ({actual} B)")
            else:
                try:
                    actual = download(url, dest, row["bytes"])
                    print(f"[OK] {gse}/{fname} ({actual} B)")
                    n_ok += 1
                except Exception as e:  # noqa: BLE001
                    print(f"[FAIL] {gse}/{fname} :: {e}")
                    val_rows.append([gse, fname, row["bytes"], "", "FAIL", str(e)[:150]])
                    n_fail += 1
                    continue
            digest = sha256_file(dest)
            man_rows.append([gse, fname, url, actual, digest])
            val_rows.append([gse, fname, row["bytes"], actual, "OK", digest])
            time.sleep(0.5)

    with open(MANIFEST, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gse", "file", "url", "bytes", "sha256"])
        w.writerows(man_rows)
    with open(VALID, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gse", "file", "expected_bytes", "actual_bytes", "status", "sha256_or_error"])
        w.writerows(val_rows)
    print(f"\ndownload done: OK={n_ok} FAIL={n_fail} -> {MANIFEST}")


if __name__ == "__main__":
    main()
