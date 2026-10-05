#!/usr/bin/env bash
# Expand the shell command of every rule for one published project, without running anything.
#
# Used by audit/export_parameter_table.py to decide which argparse flags the pipeline really passes.
# Snakemake's dry-run needs its inputs to exist, so a scratch tree is built first:
#   - every entry of the real results directory is symlinked, nothing is copied;
#   - the sample directory and the trimmed-FASTQ directory get placeholder (empty) files, which
#     exist only in the scratch tree, so that dry-run can plan the rules that need them;
#   - de/loo points at the LOO output directory, which lives outside the results directory.
# The real results directory is never written to.
#
# usage: dryrun_shell_dump.sh <project_config> <real_results_dir> <loo_dir> <scratch_root> <srr> <repo_dir>
#   scratch_root must be spelled exactly as the results_dir in the config (e.g. /home3/... not ~/...),
#   because Snakemake compares target paths as strings.
set -eu

CFG_SRC=$1; REAL=$2; LOO=$3; SCRATCH=$4; SRR=$5; REPO=$6
# SNAKEMAKE overrides the default; the default is the ciriquant conda env used for all published runs.
SM=${SNAKEMAKE:-$HOME/miniconda3/envs/ciriquant/bin/snakemake}
NAME=$(basename "$REAL")
SCR="$SCRATCH/$NAME"
TRIM="$SCRATCH/trimmed"
RAW="$SCRATCH/raw"

mirror_dir() {  # mirror_dir <real_dir> <scratch_dir> <skip_name>...
    local src=$1 dst=$2; shift 2
    mkdir -p "$dst"
    for f in "$src"/* "$src"/.[!.]*; do
        [ -e "$f" ] || [ -L "$f" ] || continue
        local b; b=$(basename "$f")
        local skip=0
        for s in "$@"; do
            if [ "$b" = "$s" ]; then skip=1; fi
        done
        if [ "$skip" = 1 ]; then continue; fi
        if [ ! -e "$dst/$b" ] && [ ! -L "$dst/$b" ]; then ln -s "$f" "$dst/$b"; fi
    done
    return 0
}

rm -rf "$SCRATCH"; mkdir -p "$SCR" "$TRIM" "$RAW"
mirror_dir "$REAL" "$SCR" de circRNA
mirror_dir "$REAL/de" "$SCR/de" loo
ln -s "$LOO" "$SCR/de/loo"
mirror_dir "$REAL/circRNA" "$SCR/circRNA" "$SRR"
mirror_dir "$REAL/circRNA/$SRR" "$SCR/circRNA/$SRR" mate1 mate2 Aligned.sortedByCoord.out.bam
mkdir -p "$SCR/circRNA/$SRR/mate1" "$SCR/circRNA/$SRR/mate2"
for p in "$SCR/circRNA/$SRR/Aligned.sortedByCoord.out.bam" \
         "$SCR/circRNA/$SRR/Chimeric.out.junction" \
         "$SCR/circRNA/$SRR/mate1/Chimeric.out.junction" \
         "$SCR/circRNA/$SRR/mate2/Chimeric.out.junction" \
         "$TRIM/${SRR}_1.fastq.gz" "$TRIM/${SRR}_2.fastq.gz" \
         "$RAW/${SRR}_1.fastq.gz" "$RAW/${SRR}_2.fastq.gz"; do
    if [ ! -e "$p" ] && [ ! -L "$p" ]; then touch "$p"; fi
done

# The scratch config differs from the project config only in the three directories it points at.
sed -e "s#^results_dir:.*#results_dir: $SCR#" \
    -e "s#^trimmed_dir:.*#trimmed_dir: $TRIM#" \
    -e "s#^raw_dir:.*#raw_dir: $RAW#" "$CFG_SRC" > "$SCRATCH/cfg.yaml"

cd "$REPO"
run_rule() {  # run_rule <rule> <target>
    echo "### RULE $1"
    # --cores: without it a dry run caps every rule at 1 thread, so the expanded commands would show
    # -t 1 / --runThreadN 1 instead of the configured threads. Dry run does not use these cores.
    timeout 300 "$SM" --snakefile workflow/Snakefile --configfile "$SCRATCH/cfg.yaml" \
        --allowed-rules "$1" --forcerun "$1" -n -p --nolock --cores 96 "$2" 2>&1 || true
    echo "### EXIT $1"
}

run_rule annotate_circbase "$SCR/circRNA/circbase_annotated.tsv"
run_rule assign_isoforms "$SCR/circRNA/isoform_groups.tsv"
run_rule check_ciriquant_config "config/.ciriquant_ready"
run_rule ciriquant "$SCR/circRNA/$SRR/$SRR.gtf"
run_rule consensus_filter "$SCR/circRNA/$SRR/high_confidence.bed"
run_rule dcc "$SCR/circRNA/$SRR/DCC/CircCoordinates"
run_rule de_analysis "$SCR/de/de_results.tsv"
run_rule download_fastq "$RAW/${SRR}_1.fastq.gz"
run_rule fastp_trim "$TRIM/${SRR}_1.fastq.gz"
run_rule fastqc_raw "$SCR/qc/fastqc_raw.done"
run_rule generate_report "$SCR/report.html"
run_rule isoform_switching "$SCR/de/isoform_switching.tsv"
run_rule loo_stability "$SCR/de/loo/run_info.tsv"
run_rule merge_counts "$SCR/circRNA/count_matrix.tsv"
run_rule multiqc "$SCR/qc/multiqc_report.html"
run_rule predict_interactions "$SCR/de/interactions.json"
run_rule rank_biomarkers "$SCR/de/biomarker_candidates.tsv"
run_rule rank_biomarkers_v2 "$SCR/de/biomarker_candidates_v2_edgeR.tsv"
run_rule star_align "$SCR/circRNA/$SRR/Chimeric.out.junction"
run_rule star_align_mate1 "$SCR/circRNA/$SRR/mate1/Chimeric.out.junction"
run_rule star_align_mate2 "$SCR/circRNA/$SRR/mate2/Chimeric.out.junction"
