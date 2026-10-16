# Recreate only query/scoring/plot outputs from completed124 fit caches.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
expressions <- parse('/work/src/program_pilot/124_HER2_RNA_coordinate_pair_decoder.R')
assignment_to <- function(x,name)is.call(x)&&identical(x[[1]],as.name('<-'))&&
  is.symbol(x[[2]])&&identical(as.character(x[[2]]),name)
start_fit <- which(vapply(expressions,assignment_to,logical(1),name='original_final'))
start_query <- which(vapply(expressions,assignment_to,logical(1),name='model'))
stopifnot(length(start_fit)==1L,length(start_query)==1L)
for(expression in expressions[seq_len(start_fit-1L)])eval(expression,envir=.GlobalEnv)
final_models <- list(normalized_cell_count=qread(file.path(root,'models/final_normalized_cell_count.qs')),
  raw_cell_count=qread(file.path(root,'models/final_raw_cell_count.qs')))
mema_pred <- fread(file.path(root,'data/mema_predictions_frozen.tsv'))
combo <- fread(file.path(root,'data/combination_predictions_frozen.tsv'))
external_metrics <- fread(file.path(root,'data/external_held_cell_correlations.tsv'))
for(expression in expressions[start_query:length(expressions)])eval(expression,envir=.GlobalEnv)
