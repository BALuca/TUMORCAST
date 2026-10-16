#!/usr/bin/env bash
set -euo pipefail
cd /work
base=almanac_conditionals/state_conditioned_bridge/program_first
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
Rscript "$base/scripts/05_select_programs.R" > "$base/artifacts/program_selection.log" 2>&1
mkdir -p "$base/sensitive_resistant/three_panels"
cp almanac_conditionals/state_conditioned_bridge/sensitivity_selected/sensitive_resistant/three_panels/lapatinib_rows.tsv "$base/sensitive_resistant/three_panels/"
cp almanac_conditionals/state_conditioned_bridge/sensitivity_selected/sensitive_resistant/three_panels/crizotinib_rows.tsv "$base/sensitive_resistant/three_panels/"
Rscript "$base/scripts/06_selected_program_heatmaps.R" > "$base/artifacts/heatmaps.log" 2>&1 &
heatmap_pid=$!
Rscript "$base/scripts/08_program_themes.R" > "$base/artifacts/themes.log" 2>&1 &
theme_pid=$!
Rscript "$base/scripts/12_gene_cross_correlations.R" > "$base/artifacts/cross_correlations.log" 2>&1 &
cross_pid=$!
Rscript "$base/scripts/07_plot_program_selection.R" > "$base/artifacts/selection_figure.log" 2>&1
/work/.venv/bin/python "$base/scripts/10_verify_frozen_bank.py" > "$base/artifacts/final_bank_audit.log" 2>&1
Rscript "$base/scripts/11_verify_source_normalization.R" > "$base/artifacts/normalization_identity.log" 2>&1
Rscript "$base/scripts/09_verify_program_first.R" > "$base/artifacts/verification.log" 2>&1
wait "$heatmap_pid"
wait "$theme_pid"

wait "$cross_pid"
