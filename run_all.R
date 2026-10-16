#!/usr/bin/env Rscript
# Run from the project root. Each stage runs in a fresh R process so one
# analysis cannot accidentally inherit another script's intermediate objects.

# Small matrix operations use one BLAS thread. Large stages explicitly set their
# own worker/thread counts, avoiding accidental oversubscription on this server.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")

scripts <- c(
  "01_prepare_expression.R",
  "02_prepare_tahoe.R",
  "03_build_gene_relationships.R",
  "04_generate_programs.R",
  "05_score_molecular_inputs.R",
  "06_prepare_functional_outcomes.R",
  "07_learn_response_signatures.R",
  "08_fit_single_drug_bridge.R",
  "09_reconstruct_hgf_sources.R",
  "10_prepare_ligand_assays.R",
  "11_transport_hgf_responses.R",
  "12_fit_ligand_rescue.R",
  "13_annotate_programs.R",
  "14_plot_model_diagnostics.R",
  "15_plot_polished_panels.R",
  "16_assemble_figures.R",
  "17_verify_reproduction.R"
)
stopifnot(file.exists("AGENTS.md"), all(file.exists(file.path("scripts", scripts))))
for (script in scripts) {
  message(format(Sys.time(), tz = "UTC"), " | Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", file.path("scripts", script)))
  if (status != 0L) stop("Stage failed: ", script, "; inspect its output before continuing")
}
message("Raw-data reconstruction and figure verification completed.")
