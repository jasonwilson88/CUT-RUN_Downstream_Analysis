#!/usr/bin/env bash
set -Eeuo pipefail

############################################
# 0) Load config
############################################
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/.env"
[[ -f "$ENV_FILE" ]] || { echo "Missing ${ENV_FILE}"; exit 1; }
# shellcheck disable=SC1090
source "$ENV_FILE"

# Optional debug
: "${DEBUG:=0}"
if [[ "${DEBUG}" == "1" ]]; then set -x; fi

############################################
# 1) Safe defaults for optional settings
############################################
: "${HEATMAP_CONTROL_TOKENS:=}"
: "${HEATMAP_TREATMENT:=}"
: "${HEATMAP_EXPLICIT_BASENAMES:=}"
: "${HEATMAP_EXPLICIT_ORDER:=}"
: "${HEATMAP_CASE_INSENSITIVE:=0}"
: "${TF_FINDPEAKS_OPTS:=-style factor -F 3 -L 3 -center}"
: "${HISTONE_FINDPEAKS_OPTS:=-style histone -size 300 -minDist 600}"
: "${ANTIBODIES:=}"        # space-separated antibodies (overrides ANTIBODY if set)
: "${CONTROL_MAP:=}"       # antibody=control;antibody=control
: "${AB_MODE_MAP:=}"       # antibody=tf|histone;antibody=tf...
: "${HIGHLIGHT_GENES:=}"   # optional; comma-separated genes to label on plots
: "${COND_PREFIXES:=}"     # e.g., "Wt Het" (optional)
: "${RUN_CONDITION_CONTRASTS:=0}"  # 1 to enable condition contrasts
: "${COND_CONTRASTS:=}"    # e.g., "Het:Wt"
: "${SORT_IDX:=}"          # deepTools --sortUsingSamples list (e.g., "3 4")

############################################
# 2) Validate required vars
############################################
require_vars() {
  local missing=0
  for v in "$@"; do
    if [[ -z "${!v+x}" ]]; then
      echo "Config error: '$v' is not set in .env" >&2; missing=1
    fi
  done
  (( missing == 0 )) || exit 2
}

# core requirements (antibody/controls checked below)
require_vars RAW ANALYSIS HOMER_ROOT BLACKLIST GENOME \
             BAM_SUFFIX THREADS PARALLEL_JOBS \
             REPS CONTRASTS CONTROL_LABELS TAGDIR_SUFFIX \
             PEAKS_UNTREATED_DIR_NAME PEAKCALL_BEDS_DIR_REL MERGED_REPS_DIR_NAME ANNOTATED_DIR_NAME PLOTS_DIR_REL \
             HELPERS HOMER_CLEAN_BED_R CENTER_R MERGE_STATS_R MERGE_REPS_R PLOTS_R

# Antibody & control requirements
if [[ -z "${ANTIBODIES:-}" && -z "${ANTIBODY:-}" ]]; then
  echo "Config error: set ANTIBODIES (space-separated) or ANTIBODY in .env" >&2; exit 2;
fi
if [[ -z "${CONTROL_DEFAULT:-}" && -z "${CONTROL:-}" && -z "${CONTROL_MAP:-}" ]]; then
  echo "Config error: set CONTROL_DEFAULT or CONTROL or CONTROL_MAP in .env" >&2; exit 2;
fi

############################################
# 3) Tools check
############################################
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 3; }; }
need bamCoverage
need bedtools
need computeMatrix
need plotHeatmap
need annotatePeaks.pl
need Rscript
need parallel

if [[ -n "${HOMER:-}" ]]; then FINDPEAKS_BIN="${HOMER}/findPeaks"; else FINDPEAKS_BIN="findPeaks"; fi
[[ -x "$FINDPEAKS_BIN" || "$(command -v "$FINDPEAKS_BIN" || true)" ]] \
  || { echo "Cannot find HOMER findPeaks (set \$HOMER or add to PATH)"; exit 3; }

log(){ printf "[%(%F %T)T] %s\n" -1 "$*" >&2; }

# Helpers exist?
for f in "$HOMER_CLEAN_BED_R" "$MERGE_STATS_R" "$MERGE_REPS_R" "$PLOTS_R"; do
  [[ -f "$f" ]] || { echo "Missing helper: $f" >&2; exit 5; }
  log "[helpers] using: $f"
done

############################################
# 4) Directories
############################################
mkdir -p \
  "${ANALYSIS}/Data/Beds" \
  "${ANALYSIS}/Data/BigWigs" \
  "${ANALYSIS}/${PEAKCALL_BEDS_DIR_REL}" \
  "${ANALYSIS}/${PLOTS_DIR_REL}/Heatmaps" \
  "${ANALYSIS}/${PLOTS_DIR_REL}/Volcanos" \
  "${ANALYSIS}/${PLOTS_DIR_REL}/pval" \
  "${ANALYSIS}/${PLOTS_DIR_REL}/FoldChange" \
  "${ANALYSIS}/Results/AnnotatedPeaks" \
  "${ANALYSIS}/Results/Matrices" \
  "${ANALYSIS}/Results/${MERGED_REPS_DIR_NAME}" \
  "${ANALYSIS}/Results/${PEAKS_UNTREATED_DIR_NAME}" \
  "${ANALYSIS}/Results"

BIGWIG_DIR="${ANALYSIS}/Data/BigWigs"
BED_DIR="${ANALYSIS}/Data/Beds"
PEAKCALL_BEDS_DIR="${ANALYSIS}/${PEAKCALL_BEDS_DIR_REL}"
MATRICES_DIR="${ANALYSIS}/Results/Matrices"
ANNOTATED_DIR="${ANALYSIS}/Results/${ANOTATED_DIR_NAME:-${ANNOTATED_DIR_NAME:-AnnotatedPeaks}}"
ANNOTATED_DIR="${ANALYSIS}/Results/${ANNOTATED_DIR_NAME}"
MERGED_REPS_DIR="${ANALYSIS}/Results/${MERGED_REPS_DIR_NAME}"
PEAKS_UNTREATED_DIR="${ANALYSIS}/Results/${PEAKS_UNTREATED_DIR_NAME}"
PLOTS_DIR="${ANALYSIS}/${PLOTS_DIR_REL}"
HEAT_DIR="${PLOTS_DIR}/Heatmaps"

############################################
# 5) Resolve antibody list, control map, mode map
############################################
declare -A CONTROL_BY_AB
if [[ -n "${CONTROL_MAP:-}" ]]; then
  IFS=';' read -r -a _pairs <<< "$CONTROL_MAP"
  for kv in "${_pairs[@]}"; do
    [[ -z "$kv" ]] && continue
    key="${kv%%=*}"; val="${kv#*=}"
    CONTROL_BY_AB["$key"]="$val"
  done
fi

declare -A AB_MODE_BY_AB
if [[ -n "${AB_MODE_MAP:-}" ]]; then
  IFS=';' read -r -a _mpairs <<< "$AB_MODE_MAP"
  for kv in "${_mpairs[@]}"; do
    [[ -z "$kv" ]] && continue
    key="${kv%%=*}"; val="${kv#*=}"
    AB_MODE_BY_AB["$key"]="$val"
  done
fi

if [[ -n "${ANTIBODIES:-}" ]]; then
  read -r -a AB_LIST <<< "$ANTIBODIES"
else
  AB_LIST=("$ANTIBODY")
fi
CONTROL_DEFAULT_EFFECTIVE="${CONTROL_DEFAULT:-${CONTROL:-IgG}}"

############################################
# 6) BAM -> BigWig (idempotent; atomic; antibody-agnostic)
############################################
log "Discovering BAMs ending with '${BAM_SUFFIX}' under ${RAW} ..."
export BLACKLIST BIGWIG_DIR BAM_SUFFIX

find "$RAW" -type f -name "*${BAM_SUFFIX}" -print0 \
| parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
  base="{/}"
  case "$base" in
    *"'"$BAM_SUFFIX"'") sample="${base%'"$BAM_SUFFIX"'}" ;;
    *)                  sample="${base%.bam}" ;;
  esac
  out="'"$BIGWIG_DIR"'/${sample}.bigWig"
  tmp="${out}.$$.tmp"
  if [ -f "$out" ] && [ ! "{}" -nt "$out" ]; then
    echo "[up-to-date BigWig] $sample"
    exit 0
  fi
  echo "[bamCoverage] $sample"
  bamCoverage -b "{}" -bl "'"$BLACKLIST"'" --normalizeUsing CPM -p "max/2" -o "$tmp"
  mv -f "$tmp" "$out"
'

############################################
# 7) Per-antibody loop
############################################
for ANTIBODY in "${AB_LIST[@]}"; do
  CONTROL_THIS="${CONTROL_BY_AB[$ANTIBODY]:-$CONTROL_DEFAULT_EFFECTIVE}"
  AB_MODE="${AB_MODE_BY_AB[$ANTIBODY]:-tf}"   # tf | histone
  export ANTIBODY CONTROL_THIS AB_MODE
  CONTROL="$CONTROL_THIS"

  log "=== Processing ANTIBODY='${ANTIBODY}' (CONTROL='${CONTROL}', MODE='${AB_MODE}') ==="

  ##########################################
  # 7.1) Generic peaks: ANTIBODY vs CONTROL
  ##########################################
  GENERIC_PEAKS_DIR="${ANALYSIS}/Results/PeakCalls_${CONTROL}_Control"
  mkdir -p "$GENERIC_PEAKS_DIR"

  export BED_DIR FINDPEAKS_BIN ANTIBODY CONTROL HOMER_CLEAN_BED_R BLACKLIST GENERIC_PEAKS_DIR TAGDIR_SUFFIX HOMER_ROOT AB_MODE TF_FINDPEAKS_OPTS HISTONE_FINDPEAKS_OPTS

  find "$HOMER_ROOT" -type d -name "*_${ANTIBODY}_*${TAGDIR_SUFFIX}" -print0 \
  | parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
    tagdir="{}"
    dirbase="{/}"
    name="${dirbase%'"$TAGDIR_SUFFIX"'}"
    rep="${name##*_}"
    if [ -z "$rep" ] || ! [[ "$rep" =~ ^[0-9]+$ ]]; then
      echo "[skip generic] bad rep in $name" >&2; exit 0
    fi
    suffix="_'"$ANTIBODY"'_$rep"
    case "$name" in
      *"$suffix") head="${name%"$suffix"}" ;;
      *) echo "[skip generic] pattern mismatch: $name" >&2; exit 0 ;;
    esac

    # head is <treatment> or <cond>_<treatment>; generic peaks treat all of head as 'treatment'
    treatment="$head"

    ctrl_base="${treatment}_'"$CONTROL"'_${rep}"
    ctrl_dir="'"$HOMER_ROOT"'/${ctrl_base}'"$TAGDIR_SUFFIX"'"
    echo "[map] TAGDIR=$dirbase -> CONTROL=${ctrl_base}'"$TAGDIR_SUFFIX"'" >&2
    [ -d "$ctrl_dir" ] || { echo "[skip generic] control dir missing: $ctrl_dir" >&2; exit 0; }

    if [ "'"$AB_MODE"'" = "histone" ]; then
      opts="'"$HISTONE_FINDPEAKS_OPTS"'"
    else
      opts="'"$TF_FINDPEAKS_OPTS"'"
    fi

    out_txt="'"$GENERIC_PEAKS_DIR"'/${name}.txt"
    out_clean="'"$GENERIC_PEAKS_DIR"'/${name}_clean.tsv"
    out_bed="'"$BED_DIR"'/${name}.bed"
    if [ -f "$out_bed" ] && [ -f "$out_clean" ]; then
      echo "[skip generic] $name"
      exit 0
    fi

    echo "[findPeaks] $name ($opts)"
    # shellcheck disable=SC2086
    "'"$FINDPEAKS_BIN"'" "$tagdir" $opts -o "$out_txt" -i "$ctrl_dir"

    echo "[clean+bed] $name"
    Rscript "'"$HOMER_CLEAN_BED_R"'" "$out_txt" "$out_clean" "$out_bed" "'"$BLACKLIST"'" --schema auto
  '

  #####################################################
  # 7.2) Merge + center generic beds -> matrix + heatmap
  #####################################################
  log "Merging generic BEDs for antibody '${ANTIBODY}' ..."
  tmp_all="${BED_DIR}/All_${ANTIBODY}.bed"
  tmp_sorted="${BED_DIR}/AllSorted_${ANTIBODY}.bed"
  merged="${BED_DIR}/AllMerged_${ANTIBODY}.bed"
  centered="${BED_DIR}/AllMergedCentered_${ANTIBODY}.bed"

  : > "$tmp_all"
  while IFS= read -r -d '' f; do
    cat "$f" >> "$tmp_all"
  done < <(find "$BED_DIR" -maxdepth 1 -type f -name "*_${ANTIBODY}_*.bed" -size +0c -print0)

  if [[ ! -s "$tmp_all" ]]; then
    log "[warn] No BEDs found for antibody '${ANTIBODY}'. Skipping heatmap."
  else
    bedtools sort  -i "$tmp_all"    > "$tmp_sorted"
    bedtools merge -i "$tmp_sorted" > "$merged"
    rm -f "$tmp_all" "$tmp_sorted"

    log "Centering merged peaks..."
    Rscript "$CENTER_R" "$merged" "$centered"
    rm -f "$merged"

    # Collect BigWigs for this antibody
    log "Collecting BigWigs for antibody '${ANTIBODY}' ..."
    mapfile -d '' BW_ARR < <(find "$BIGWIG_DIR" -maxdepth 1 -type f -name "*${ANTIBODY}*.bigWig" -print0)
    if (( ${#BW_ARR[@]} == 0 )); then
      log "[warn] No BigWigs matched *${ANTIBODY}*.bigWig; skipping heatmap."
    else
      # -------- Heatmap ordering (explicit-first; robust token parsing) --------
      declare -A _used=()
      declare -A _tokrep_to_file=()
      ordered=()

      # Map token|rep -> file via tail parse: *_<TOKEN>_<ANTIBODY>_<rep>.bigWig
      for f in "${BW_ARR[@]}"; do
        fb="$(basename "$f")"
        if [[ "$fb" =~ _${ANTIBODY}_([0-9]+)\.bigWig$ ]]; then
          rep="${BASH_REMATCH[1]}"
          tail="_${ANTIBODY}_${rep}.bigWig"
          prefix="${fb%$tail}"
          token="${prefix##*_}"
          key="${token,,}|${rep}"
          _tokrep_to_file["$key"]="$f"
        fi
      done

      if [[ -n "${HEATMAP_EXPLICIT_BASENAMES}" ]]; then
        IFS=',;' read -r -a _bases <<< "${HEATMAP_EXPLICIT_BASENAMES}"
        for b in "${_bases[@]}"; do
          b="${b%.bigWig}.bigWig"
          for f in "${BW_ARR[@]}"; do
            fb="$(basename "$f")"
            if [[ -z "${_used[$f]+x}" && "$fb" == "$b" ]]; then
              ordered+=("$f"); _used["$f"]=1; break
            fi
          done
        done
      elif [[ -n "${HEATMAP_EXPLICIT_ORDER}" ]]; then
        IFS='|' read -r -a _groups <<< "${HEATMAP_EXPLICIT_ORDER}"
        for g in "${_groups[@]}"; do
          token="${g%%:*}"; reps="${g#*:}"
          [[ "$reps" == "$g" ]] && reps="1,2"
          IFS=',' read -r -a _reps <<< "$reps"
          ktok="${token,,}"
          for r in "${_reps[@]}"; do
            key="${ktok}|${r}"
            f="${_tokrep_to_file[$key]-}"
            if [[ -n "${f:-}" && -z "${_used[$f]+x}" ]]; then
              ordered+=("$f"); _used["$f"]=1
            fi
          done
        done
      elif [[ -n "${HEATMAP_CONTROL_TOKENS}" || -n "${HEATMAP_TREATMENT}" ]]; then
        declare -A _seenlegacy=()
        _append_by_rx() {
          local rx="$1" ff fb
          for ff in "${BW_ARR[@]}"; do
            [[ -n "${_seenlegacy[$ff]+x}" ]] && continue
            fb="$(basename "$ff")"
            if [[ "$fb" =~ $rx ]]; then
              ordered+=("$ff"); _seenlegacy["$ff"]=1
            fi
          done
        }
        [[ -n "${HEATMAP_CONTROL_TOKENS}" ]] && _append_by_rx "(${HEATMAP_CONTROL_TOKENS}).*_[0]*1\.bigWig$"
        [[ -n "${HEATMAP_CONTROL_TOKENS}" ]] && _append_by_rx "(${HEATMAP_CONTROL_TOKENS}).*_[0]*2\.bigWig$"
        [[ -n "${HEATMAP_TREATMENT}"     ]] && _append_by_rx "(${HEATMAP_TREATMENT}).*_[0]*1\.bigWig$"
        [[ -n "${HEATMAP_TREATMENT}"     ]] && _append_by_rx "(${HEATMAP_TREATMENT}).*_[0]*2\.bigWig$"
        _append_by_rx "(Combo|COMBO).*_[0]*1\.bigWig$"
        _append_by_rx "(Combo|COMBO).*_[0]*2\.bigWig$"
      fi

      # Append remaining in discovery order
      for f in "${BW_ARR[@]}"; do
        [[ -z "${_used[$f]+x}" ]] && ordered+=("$f")
      done

      if ((${#ordered[@]})); then
        printf '[heatmap order %s] %s\n' "$ANTIBODY" "${ordered[@]##*/}"
        # Define outputs BEFORE skip checks to avoid unbound vars
	# Define outputs first (so they're never unbound)
	MATRIX_OUT="${MATRICES_DIR}/${ANTIBODY}.matrix.gz"
	HEATMAP_OUT="${HEAT_DIR}/$(echo "${ANTIBODY}" | sed 's/ /_/g')_heatmap.pdf"

	# --- computeMatrix: only if matrix doesn't already exist ---
	if [[ -f "$MATRIX_OUT" ]]; then
	  log "[skip computeMatrix] exists: ${MATRIX_OUT}"
	else
	  log "computeMatrix -> ${MATRIX_OUT}"
	  if [[ -n "${SORT_IDX}" ]]; then
	    # shellcheck disable=SC2086
	    computeMatrix reference-point -R "$centered" -S "${ordered[@]}" \
	      -b 2500 -a 2500 -bs 50 --numberOfProcessors "$THREADS" \
	      --missingDataAsZero --sortUsingSamples $SORT_IDX \
	      -o "$MATRIX_OUT"
	  else
	    computeMatrix reference-point -R "$centered" -S "${ordered[@]}" \
	      -b 2500 -a 2500 -bs 50 --numberOfProcessors "$THREADS" \
	      --missingDataAsZero \
	      -o "$MATRIX_OUT"
	  fi
	fi

	# --- plotHeatmap: only if heatmap doesn't already exist ---
	if [[ -f "$HEATMAP_OUT" ]]; then
	  log "[skip plotHeatmap] exists: ${HEATMAP_OUT}"
	else
	  log "plotHeatmap -> ${HEATMAP_OUT}"
	  plotHeatmap -m "$MATRIX_OUT" -out "$HEATMAP_OUT" \
	    --colorList white,steelblue --sortUsing mean --refPointLabel ""
	fi
      else
        printf '[heatmap order %s] <none>\n' "$ANTIBODY"
      fi
    fi
  fi

  ###################################################
  # 7.3) Custom treatment contrasts (e.g., pIC:DMSO)
  ###################################################
  log "Running treatment contrasts for antibody '${ANTIBODY}' ..."
  mkdir -p "$PEAKS_UNTREATED_DIR" "$PEAKCALL_BEDS_DIR" "$ANNOTATED_DIR" "$MERGED_REPS_DIR"

  declare -A CTRL_ALIAS=()
  IFS=';' read -r -a __pairs <<< "${CONTROL_LABELS}"
  for kv in "${__pairs[@]}"; do key="${kv%%=*}"; val="${kv#*=}"; CTRL_ALIAS["$key"]="$val"; done

  TMP_TSV="$(mktemp)"
  {
    for pair in ${CONTRASTS}; do
      q="${pair%%:*}"; c="${pair#*:}"; ctrl_label="${CTRL_ALIAS[$c]:-$c}"
      for rep in ${REPS}; do
        tag_q="${HOMER_ROOT}/${q}_${ANTIBODY}_${rep}${TAGDIR_SUFFIX}"
        tag_c="${HOMER_ROOT}/${c}_${ANTIBODY}_${rep}${TAGDIR_SUFFIX}"
        out_txt="${PEAKS_UNTREATED_DIR}/${ANTIBODY}_${q}_vs_${ctrl_label}_${rep}.txt"
        [[ -d "$tag_q" && -d "$tag_c" ]] || { echo "WARN: missing $tag_q or $tag_c" >&2; continue; }
        [[ -f "$out_txt" ]] && { echo "[skip custom findPeaks] $(basename "$out_txt")" >&2; continue; }
        printf '%s\t%s\t%s\n' "$tag_q" "$tag_c" "$out_txt"
      done
    done
  } > "$TMP_TSV"

  if [[ -s "$TMP_TSV" ]]; then
    export AB_MODE TF_FINDPEAKS_OPTS HISTONE_FINDPEAKS_OPTS FINDPEAKS_BIN
    parallel -j "${PARALLEL_JOBS}" --colsep '\t' --bar --halt now,fail=1 '
      q="{1}"; c="{2}"; out="{3}";
      if [ "'"$AB_MODE"'" = "histone" ]; then
        opts="'"$HISTONE_FINDPEAKS_OPTS"'"
      else
        opts="'"$TF_FINDPEAKS_OPTS"'"
      fi
      echo "[custom findPeaks] $(basename "$q") vs $(basename "$c") -> $out ($opts)"
      # shellcheck disable=SC2086
      '"$FINDPEAKS_BIN"' "$q" $opts -o "$out" -i "$c"
    ' :::: "$TMP_TSV"
  fi
  rm -f "$TMP_TSV"

  # Clean + BED (all custom)
  export HOMER_CLEAN_BED_R BLACKLIST PEAKCALL_BEDS_DIR
  find "$PEAKS_UNTREATED_DIR" -type f -name "${ANTIBODY}_*.txt" -print0 \
  | parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
      in="{}"
      base="{/.}"
      out_clean="{.}_clean.tsv"
      out_bed="'"$PEAKCALL_BEDS_DIR"'/${base}.bed"
      if [ -f "$out_bed" ] && [ -f "$out_clean" ]; then
        echo "[skip clean+bed] $base"
        exit 0
      fi
      echo "[clean+bed] $base"
      Rscript "'"$HOMER_CLEAN_BED_R"'" "$in" "$out_clean" "$out_bed" "'"$BLACKLIST"'" --schema auto
  '

  # Annotate (all custom)
  export GENOME ANNOTATED_DIR
  find "$PEAKCALL_BEDS_DIR" -type f -name "${ANTIBODY}_*.bed" -print0 \
  | parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
      bed="{}"
      base="{/.}"
      out="'"$ANNOTATED_DIR"'/${base}.txt"
      if [ -f "$out" ]; then
        echo "[skip annotate] $base"
        exit 0
      fi
      echo "[annotatePeaks] $base"
      annotatePeaks.pl "$bed" "'"$GENOME"'" > "$out"
  '

  ###################################################
  # 7.4) Merge & plot (treatment contrasts)
  ###################################################
  for pair in ${CONTRASTS}; do
    q="${pair%%:*}"
    c="${pair#*:}"
    declare -A CTRL_ALIAS=()
    IFS=';' read -r -a __pairs <<< "${CONTROL_LABELS}"
    for kv in "${__pairs[@]}"; do key="${kv%%=*}"; val="${kv#*=}"; CTRL_ALIAS["$key"]="$val"; done
    ctrl_label="${CTRL_ALIAS[$c]:-$c}"
    base="${ANTIBODY}_${q}_vs_${ctrl_label}"

    cat_tmp="${PEAKCALL_BEDS_DIR}/${base}.cat.bed"
    sort_tmp="${PEAKCALL_BEDS_DIR}/${base}.sort.bed"
    merged_bed="${PEAKCALL_BEDS_DIR}/${base}.merged.bed"
    merged_anno="${ANNOTATED_DIR}/${base}_Merged.txt"

    : > "$cat_tmp"
    have_any=0
    rep_stats=()
    for rep in ${REPS}; do
      bed="${PEAKCALL_BEDS_DIR}/${base}_${rep}.bed"
      clean_tsv="${PEAKS_UNTREATED_DIR}/${base}_${rep}_clean.tsv"
      anno_txt="${ANNOTATED_DIR}/${base}_${rep}.txt"
      out_stats="${MERGED_REPS_DIR}/${base}_${rep}.tsv"

      if [[ -s "$bed" ]]; then
        cat "$bed" >> "$cat_tmp"
        have_any=1
      fi

      if [[ -f "$clean_tsv" && -f "$anno_txt" ]]; then
        Rscript "$MERGE_STATS_R" "$anno_txt" "$clean_tsv" "$out_stats"
        rep_stats+=("$out_stats")
      fi
    done

    (( have_any )) || { rm -f "$cat_tmp"; continue; }

    bedtools sort  -i "$cat_tmp" > "$sort_tmp"
    bedtools merge -i "$sort_tmp" > "$merged_bed"
    rm -f "$sort_tmp" "$cat_tmp"

    annotatePeaks.pl "$merged_bed" "$GENOME" > "$merged_anno"

    if (( ${#rep_stats[@]} >= 2 )); then
      overlap="${MERGED_REPS_DIR}/${base}_Overlap.tsv"
      Rscript "$MERGE_REPS_R" "$overlap" "${rep_stats[@]}"

      pval_pdf="${PLOTS_DIR}/pval/${base}_pval.pdf"
      fc_pdf="${PLOTS_DIR}/FoldChange/${base}_FC.pdf"
      volcano1_pdf="${PLOTS_DIR}/Volcanos/${base}_Rep1Volcano.pdf"
      volcano2_pdf="${PLOTS_DIR}/Volcanos/${base}_Rep2Volcano.pdf"

      if [[ -n "${HIGHLIGHT_GENES}" ]]; then
        Rscript "$PLOTS_R" "$overlap" "$pval_pdf" "$fc_pdf" "$volcano1_pdf" "$volcano2_pdf" --highlight="${HIGHLIGHT_GENES}"
      else
        Rscript "$PLOTS_R" "$overlap" "$pval_pdf" "$fc_pdf" "$volcano1_pdf" "$volcano2_pdf"
      fi
    fi
  done

  ###################################################
  # 7.5) Condition contrasts (e.g., Het vs Wt)
  ###################################################
  if [[ "${RUN_CONDITION_CONTRASTS:-0}" == "1" && -n "${COND_CONTRASTS:-}" ]]; then
    log "Running condition contrasts for antibody '${ANTIBODY}' ..."
    mkdir -p "$PEAKS_UNTREATED_DIR" "$PEAKCALL_BEDS_DIR" "$ANNOTATED_DIR" "$MERGED_REPS_DIR"

    # Build map of available tagdirs: key = "<cond>|<treatment>|<rep>"
    declare -A TD_MAP=()
    declare -A TREATMENT_SET=()
    read -r -a _COND_LIST <<< "${COND_PREFIXES:-}"

    while IFS= read -r -d '' d; do
      b="$(basename "$d")"
      name="${b%$TAGDIR_SUFFIX}"
      rep="${name##*_}"
      [[ -n "$rep" && "$rep" =~ ^[0-9]+$ ]] || continue
      suffix="_${ANTIBODY}_${rep}"
      [[ "$name" == *"$suffix" ]] || continue
      head="${name%$suffix}"

      cond="" ; treatment="$head"
      for cp in "${_COND_LIST[@]:-}"; do
        if [[ -n "$cp" && "$head" == "${cp}_"* ]]; then
          cond="$cp"
          treatment="${head#${cp}_}"
          break
        fi
      done

      key="${cond}|${treatment}|${rep}"
      TD_MAP["$key"]="$d"
      TREATMENT_SET["$treatment"]=1
    done < <(find "$HOMER_ROOT" -type d -name "*_${ANTIBODY}_*${TAGDIR_SUFFIX}" -print0)

    TMP_COND_TSV="$(mktemp)"
    {
      for cpair in ${COND_CONTRASTS}; do
        qcond="${cpair%%:*}"; ccond="${cpair#*:}"
        for treatment in "${!TREATMENT_SET[@]}"; do
          for rep in ${REPS}; do
            q_key="${qcond}|${treatment}|${rep}"
            c_key="${ccond}|${treatment}|${rep}"
            q_dir="${TD_MAP[$q_key]-}"; c_dir="${TD_MAP[$c_key]-}"
            [[ -n "$q_dir" && -n "$c_dir" ]] || continue
            base="${ANTIBODY}_${treatment}_${qcond}_vs_${ccond}_${rep}"
            out_txt="${PEAKS_UNTREATED_DIR}/${base}.txt"
            [[ -f "$out_txt" ]] && { echo "[skip cond findPeaks] $(basename "$out_txt")" >&2; continue; }
            printf '%s\t%s\t%s\n' "$q_dir" "$c_dir" "$out_txt"
          done
        done
      done
    } > "$TMP_COND_TSV"

    if [[ -s "$TMP_COND_TSV" ]]; then
      export AB_MODE TF_FINDPEAKS_OPTS HISTONE_FINDPEAKS_OPTS FINDPEAKS_BIN
      parallel -j "${PARALLEL_JOBS}" --colsep '\t' --bar --halt now,fail=1 '
        q="{1}"; c="{2}"; out="{3}";
        if [ "'"$AB_MODE"'" = "histone" ]; then
          opts="'"$HISTONE_FINDPEAKS_OPTS"'"
        else
          opts="'"$TF_FINDPEAKS_OPTS"'"
        fi
        echo "[cond findPeaks] $(basename "$q") vs $(basename "$c") -> $out ($opts)"
        # shellcheck disable=SC2086
        '"$FINDPEAKS_BIN"' "$q" $opts -o "$out" -i "$c"
      ' :::: "$TMP_COND_TSV"
    fi
    rm -f "$TMP_COND_TSV"

    # Clean + BED (condition)
    export HOMER_CLEAN_BED_R BLACKLIST PEAKCALL_BEDS_DIR
    find "$PEAKS_UNTREATED_DIR" -type f -name "${ANTIBODY}_*.txt" -print0 \
    | parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
        in="{}"
        base="{/.}"
        out_clean="{.}_clean.tsv"
        out_bed="'"$PEAKCALL_BEDS_DIR"'/${base}.bed"
        if [ -f "$out_bed" ] && [ -f "$out_clean" ]; then
          echo "[skip clean+bed] $base"
          exit 0
        fi
        echo "[clean+bed] $base"
        Rscript "'"$HOMER_CLEAN_BED_R"'" "$in" "$out_clean" "$out_bed" "'"$BLACKLIST"'" --schema auto
    '

    # Annotate (condition)
    export GENOME ANNOTATED_DIR
    find "$PEAKCALL_BEDS_DIR" -type f -name "${ANTIBODY}_*.bed" -print0 \
    | parallel -0 -j "${PARALLEL_JOBS}" --bar --halt now,fail=1 '
        bed="{}"
        base="{/.}"
        out="'"$ANNOTATED_DIR"'/${base}.txt"
        if [ -f "$out" ]; then
          echo "[skip annotate] $base"
          exit 0
        fi
        echo "[annotatePeaks] $base"
        annotatePeaks.pl "$bed" "'"$GENOME"'" > "$out"
    '

    # Merge & plot (condition)
    for cpair in ${COND_CONTRASTS}; do
      qcond="${cpair%%:*}"; ccond="${cpair#*:}"
      for treatment in "${!TREATMENT_SET[@]}"; do
        base="${ANTIBODY}_${treatment}_${qcond}_vs_${ccond}"

        cat_tmp="${PEAKCALL_BEDS_DIR}/${base}.cat.bed"
        sort_tmp="${PEAKCALL_BEDS_DIR}/${base}.sort.bed"
        merged_bed="${PEAKCALL_BEDS_DIR}/${base}.merged.bed"
        merged_anno="${ANNOTATED_DIR}/${base}_Merged.txt"

        : > "$cat_tmp"
        have_any=0
        rep_stats=()
        for rep in ${REPS}; do
          bed="${PEAKCALL_BEDS_DIR}/${base}_${rep}.bed"
          clean_tsv="${PEAKS_UNTREATED_DIR}/${base}_${rep}_clean.tsv"
          anno_txt="${ANNOTATED_DIR}/${base}_${rep}.txt"
          out_stats="${MERGED_REPS_DIR}/${base}_${rep}.tsv"

          if [[ -s "$bed" ]]; then
            cat "$bed" >> "$cat_tmp"
            have_any=1
          fi
          if [[ -f "$clean_tsv" && -f "$anno_txt" ]]; then
            Rscript "$MERGE_STATS_R" "$anno_txt" "$clean_tsv" "$out_stats"
            rep_stats+=("$out_stats")
          fi
        done

        (( have_any )) || { rm -f "$cat_tmp"; continue; }

        bedtools sort  -i "$cat_tmp" > "$sort_tmp"
        bedtools merge -i "$sort_tmp" > "$merged_bed"
        rm -f "$sort_tmp" "$cat_tmp"

        annotatePeaks.pl "$merged_bed" "$GENOME" > "$merged_anno"

        if (( ${#rep_stats[@]} >= 2 )); then
          overlap="${MERGED_REPS_DIR}/${base}_Overlap.tsv"
          Rscript "$MERGE_REPS_R" "$overlap" "${rep_stats[@]}"

          pval_pdf="${PLOTS_DIR}/pval/${base}_pval.pdf"
          fc_pdf="${PLOTS_DIR}/FoldChange/${base}_FC.pdf"
          volcano1_pdf="${PLOTS_DIR}/Volcanos/${base}_Rep1Volcano.pdf"
          volcano2_pdf="${PLOTS_DIR}/Volcanos/${base}_Rep2Volcano.pdf"

          if [[ -n "${HIGHLIGHT_GENES}" ]]; then
            Rscript "$PLOTS_R" "$overlap" "$pval_pdf" "$fc_pdf" "$volcano1_pdf" "$volcano2_pdf" --highlight="${HIGHLIGHT_GENES}"
          else
            Rscript "$PLOTS_R" "$overlap" "$pval_pdf" "$fc_pdf" "$volcano1_pdf" "$volcano2_pdf"
          fi
        fi
      done
    done
  else
    log "Condition contrasts disabled (RUN_CONDITION_CONTRASTS=${RUN_CONDITION_CONTRASTS:-0})."
  fi

done # end per-antibody loop

log "All steps complete."

