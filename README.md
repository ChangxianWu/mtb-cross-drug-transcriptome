# Cross-drug integration of antibiotic-perturbed transcriptomes in *Mycobacterium tuberculosis*

Analysis pipeline for the manuscript **"Cross-drug integration of antibiotic-perturbed
transcriptomes maps metabolic response programs in *Mycobacterium tuberculosis*"
(submitted to *Frontiers in Microbiology*, 2026-11-12).**

This repository contains the **analysis code only** — the scripts that download the
public GEO data, build the module scores, estimate the effect layers, calibrate
inference and draw the figures. Manuscript-assembly and reference-auditing scripts
are not part of this repository; all result tables quoted in the paper are produced
by the scripts below and are archived by the authors.

## Layout

```
code/   Python analysis pipeline, scripts 01-11 + mtblib.py (execution order = filename number)
```

| Script | Role |
|---|---|
| 01 | GEO inventory of candidate series (metadata + FTP listings) |
| 02 | Selective download with byte-length checks, manifest, SHA-256 |
| 03 | Unified parsing to per-series expression matrices + standardized sample tables |
| 04 | Positive-control gate (runs before any primary analysis) |
| 05 | Module scoring (mean z of member genes, per dataset) |
| 06 | Cross-dataset effect sizes + DerSimonian-Laird random-effects meta-analysis |
| 07 | Gene-level effect layer (Mode A2 robustness arm) |
| 08 | Empirical null from size-matched random gene sets (negative controls) |
| 09 | Module-by-drug kill-class association |
| 10 | Sensitivity suite (leave-one-dataset-out, dose strata, rank scores, cross-timepoint) |
| 11 | Figures F1-F8 (300 dpi) |
| mtblib.py | Shared statistics library (Hedges g, DL meta, BH-FDR, KEGG panel, loaders) |

## Data flow

```
GEO series ──01──► inventory ──02──► raw downloads ──03──► parsed matrices
                ──04──► gate check ──05──► module scores
                ──06/07──► effect layers (A1 / A2) ──08──► empirical null
                ──09/10──► kill-class + sensitivity ──11──► figures
```

## Environment

Python 3.13 with `pandas`, `numpy`, `scipy`, `matplotlib`. All inputs are public GEO
series; accessions are listed in Table 1 of the manuscript.

## Running

1. Scripts resolve their workspace from `BASE` (defined in `code/mtblib.py` as the
   parent directory of `code/`); outputs are written beneath `BASE` as
   `data/`, `output/`, `cache/` and `figures/`.
2. Run `code/01` … `code/11` in filename order; each step writes its logs and tables
   under `output/`.
3. The pipeline is deterministic: re-running scripts 03-11 reproduces the archived
   result tables byte-identically (no random seeds in scripts 01-11).

## Manuscript

- Submitted to *Frontiers in Microbiology* (Original Research), November 2026.
- Correspondence: addresses in the manuscript.
