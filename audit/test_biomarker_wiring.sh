#!/usr/bin/env bash
# Dry-run assertion test for rank_biomarkers_v2 wiring (read-only; no output is written).
#
# Every de.biomarker key that deviates from its default must appear in the rank_biomarkers_v2
# command line; the control (no overlay) must contain none of the optional flags.
#
# usage: test_biomarker_wiring.sh <scratch_config> <scratch_target> <repo_dir> <work_dir>
#   scratch_config  a config whose results_dir is a scratch tree (see audit/dryrun_shell_dump.sh)
#   scratch_target  the rank_biomarkers_v2 output path inside that tree
set -u
CFG=$1; TARGET=$2; REPO=$3; WORK=$4
SM=${SNAKEMAKE:-$HOME/miniconda3/envs/ciriquant/bin/snakemake}
mkdir -p "$WORK"
cd "$REPO" || exit 2
fail=0

check() {  # check <name> <config file to use> <must-have flag>... -- <must-not flag>...
  local name=$1 cfgfile=$2; shift 2
  local must=() mustnot=() seen=0
  for a in "$@"; do
    if [ "$a" = "--" ]; then seen=1; continue; fi
    if [ $seen = 0 ]; then must+=("$a"); else mustnot+=("$a"); fi
  done
  "$SM" --snakefile workflow/Snakefile --configfile "$cfgfile" \
      --allowed-rules rank_biomarkers_v2 --forcerun rank_biomarkers_v2 \
      -n -p --nolock "$TARGET" > "$WORK/$name.txt" 2>&1 || true
  local cmd
  cmd=$(grep "scripts/rank_biomarkers.py" "$WORK/$name.txt" | tail -1)
  if [ -z "$cmd" ]; then echo "FAIL $name: no rank_biomarkers_v2 command in dry-run (see $WORK/$name.txt)"; fail=1; return; fi
  local ok=1
  for f in ${must[@]+"${must[@]}"}; do
    if ! grep -q -- "$f" <<< "$cmd"; then echo "FAIL $name: missing '$f'"; ok=0; fi
  done
  for f in ${mustnot[@]+"${mustnot[@]}"}; do
    if grep -q -- "$f" <<< "$cmd"; then echo "FAIL $name: unexpected '$f'"; ok=0; fi
  done
  if [ $ok = 1 ]; then echo "PASS $name"; else fail=1; fi
}

# A variant is a full copy of the base config with one de.biomarker key changed. A separate overlay
# file is not used: Snakemake replaces (not deep-merges) top-level keys across --configfile files.
PY=${PYTHON:-$HOME/miniconda3/envs/ciriquant/bin/python}
variant() {  # variant <name> <dotted key path> <value as YAML>
  $PY - "$CFG" "$WORK/cfg_$1.yaml" "$2" "$3" <<'PYEOF'
import sys, yaml
src, dst, key, val = sys.argv[1:5]
cfg = yaml.safe_load(open(src))
node = cfg
parts = key.split(".")
for p in parts[:-1]:
    node = node.setdefault(p, {})
node[parts[-1]] = yaml.safe_load(val)
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False)
PYEOF
  echo "$WORK/cfg_$1.yaml"
}

check control "$CFG" "--cross-method-min 3" "--stability-min-ratio 0.5" -- \
      "--sig-cap" "--fc-cap" "--stability-max-candidates" "--components-json"

check sig_cap "$(variant sig_cap de.biomarker.sig_cap 8)" "--sig-cap 8" -- "--fc-cap"
check fc_cap "$(variant fc_cap de.biomarker.fc_cap 4)" "--fc-cap 4" -- "--sig-cap"
check stability_max_candidates "$(variant smc de.biomarker.stability_max_candidates 300)" "--stability-max-candidates 300" --
check cross_method_min "$(variant cmm de.biomarker.cross_method_min 2)" "--cross-method-min 2" --
check stability_min_ratio "$(variant smr de.biomarker.stability_min_ratio 0.6)" "--stability-min-ratio 0.6" --
check components "$(variant comp de.biomarker.components '{sig: {enabled: true, weight: 2.0}, fc: {enabled: true, weight: 1.0}, conf: {enabled: true, weight: 1.0}, known: {enabled: true, weight: 1.0}, mirna: {enabled: true, weight: 1.0}, rbp: {enabled: true, weight: 1.0}}')" "--components-json" --

if [ $fail = 0 ]; then echo "ALL_PASS"; else echo "SOME_FAILED"; exit 1; fi
