# CUT&RUN Pipeline (HOMER + deepTools + bedtools + GNU Parallel)

Automated CUT&RUN processing + peak calling + blacklisting + BED/annotation +
replicate merging + square plots — with multi-antibody support and optional
condition contrasts (e.g., `Wt` vs `Het`).

- **Robust inputs**: handles TF *and* histone modes per antibody.
- **Idempotent** reruns: BigWigs + matrices + heatmaps can be skipped.
- **Explicit heatmap order** via `.env`.
- **Always-labeled highlights** on plots (`HIGHLIGHT_GENES`).
- **Two-rep overlap** (histone-friendly): pure R, preserves per-rep coords.

---

## Repo layout

```
.
├─ run_cutrun_all.sh              # main pipeline (bash; GNU parallel)
├─ .env.example                   # copy to `.env` and customize
├─ HelperScripts/
│  ├─ homer_clean_and_bed.R       # blacklist + standardize + BED
│  ├─ center_bed.R                # center intervals to 1bp
│  ├─ merge_stats.R               # join stats↔annotation; keeps all columns
│  ├─ merge_replicates.R          # two-rep, pure R overlap; keeps per-rep coords
│  └─ Plots.R                     # square plots; volcano x=-log10(p), y=FC
├─ scripts/
│  └─ install_r_deps.R            # installs required R packages
├─ .gitignore
└─ LICENSE
```

> **Note:** File names are **case-sensitive** on Linux. Ensure your `.env` paths match.

---

## Requirements

System tools
- **bash** (≥4), **GNU parallel**, **bedtools**
- **deepTools**: `bamCoverage`, `computeMatrix`, `plotHeatmap`
- **HOMER**: `findPeaks`, `annotatePeaks.pl` (ensure `HOMER` is on PATH or set `$HOMER`)
- **R** (≥4.0) with packages: `dplyr`, `readr`, `stringr`, `ggplot2`, `ggrepel`

Install R deps:
```bash
Rscript scripts/install_r_deps.R
```

---

## Quickstart

1. **Copy and edit** the env template:
   ```bash
   cp .env.example .env
   # open .env and fill in RAW, ANALYSIS, HOMER_ROOT, BLACKLIST, GENOME, etc.
   ```

2. **Make the pipeline executable**:
   ```bash
   chmod +x run_cutrun_all.sh
   ```

3. **Run**:
   ```bash
   ./run_cutrun_all.sh
   ```

4. **Re-run safely**: BigWigs are skipped if present; matrices and heatmaps are
   skipped if the files already exist (see `.env` and comments in the script).

---

## Configuration (in `.env`)

Key knobs (see `.env.example` for full list):

- **Roots & inputs**
  - `RAW`, `ANALYSIS`, `HOMER_ROOT`, `BLACKLIST`, `GENOME`, `BAM_SUFFIX`
- **Antibodies + control**
  - `ANTIBODIES="RelA_monoclonal MCM2"` (space-separated)
  - `CONTROL_DEFAULT=IgG`, `CONTROL_MAP="MCM2=IgG;RelA_monoclonal=IgG"`
  - `AB_MODE_MAP="RelA_monoclonal=tf;MCM2=histone"`
- **Replicates & contrasts**
  - `REPS="1 2"`
  - `CONTRASTS="pIC:DMSO pIC:Combo"` (treatment contrasts)
- **Condition contrasts (optional)**
  - `COND_PREFIXES="Wt Het"`
  - `RUN_CONDITION_CONTRASTS=1`
  - `COND_CONTRASTS="Het:Wt"`
- **Heatmap order**
  - `HEATMAP_EXPLICIT_ORDER="DMSO:1,2|TNFa:1,2|TPCA1:1,2|Combo:1,2"`
  - or `HEATMAP_EXPLICIT_BASENAMES="Sample_DMSO_${ANTIBODY}_1.bigWig;…"`
  - `SORT_IDX="3 4"` (deepTools row sorting by these track indices)
- **Plot highlights & square size**
  - `HIGHLIGHT_GENES="IFNAR2,STAT1"`
  - `PLOT_SIZE_IN=6`
- **Overlap controls (merge_replicates)**
  - `OVERLAP_PAD_BP=200` (histone-friendly)
  - `OVERLAP_REQUIRE_SAME_GENE=0`
  - `MERGE_VERBOSE=1`

---

## Outputs (under `$ANALYSIS`)

```
Data/
  Beds/                 # per-treatment BEDs + merged/centered
  BigWigs/
  PeakCallBeds/         # per-contrast BEDs (custom/conditions)
Results/
  AnnotatedPeaks/       # annotatePeaks outputs
  Matrices/             # deepTools matrices (per antibody)
  MergedRepsWithStats/  # per-rep merged stats + Overlap.tsv
  PeakCalls_*_Control/  # generic peaks vs control (by antibody)
Plots/
  Heatmaps/
  FoldChange/
  pval/
  Volcanos/
```

---

## Troubleshooting

- **No gene labels / `gene_name` NA**: the annotation file must include `"Gene Name"`; `merge_stats.R` renames it to `gene_name` and propagates it through. Ensure your `.env` points to this helper and file casing matches.
- **Empty overlap** on histone: set `OVERLAP_PAD_BP=200` and rerun the overlap step (the overlap script preserves per-rep coords and FC/P).
- **Skip expensive steps**: matrices/heatmaps are skipped if files already exist.
- **Case-sensitive paths**: Linux is strict — double-check `.env` and file names.

---

## Citation

Please cite the tools this wraps: **HOMER**, **deepTools**, **bedtools**, **GNU parallel**, and the R packages **dplyr**, **readr**, **stringr**, **ggplot2/ggrepel**.

---

## License

MIT — see [LICENSE](LICENSE).
