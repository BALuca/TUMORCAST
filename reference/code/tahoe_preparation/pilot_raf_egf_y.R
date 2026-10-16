#!/usr/bin/env Rscript

## RAF/EGF pilot: discovery in Tahoe + PRISM + CytoSig EGF,
## followed by frozen-program molecular validation in GSE127988.
if(file.exists("/work/renv/activate.R") && !requireNamespace("RhpcBLASctl",quietly=TRUE)) source("/work/renv/activate.R")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cluster)
})
source("/work/src/R/lib_heatmap.R")
RhpcBLASctl::blas_set_num_threads(8L)
root <- "/work/tahoe_prism_v2"
cache <- "/work/cache/tahoe_prism_full_chain"
pdex <- file.path(cache, "tahoe_pdex_selected")
out <- file.path(root, "results/24_raf_egf_y")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
drugs <- c("vemurafenib", "cobimetinib")
paths <- c(
  effect_metadata = file.path(root, "artifacts/01_inputs/response_metadata.tsv.gz"),
  all_genes = file.path(pdex, "genes.tsv"),
  normalization = file.path(root, "artifacts/02_response_programs/normalization_records.tsv.gz"),
  control_effects = file.path(pdex, "pdex_plate_reference.npy"),
  control_metadata = file.path(pdex, "pdex_plate_reference_metadata.tsv"),
  prism = file.path(root, "artifacts/01_inputs/function_outcomes.tsv.gz"),
  split = file.path(root, "results/01_inputs/fold_assignments.tsv"),
  cytosig_metadata = "/work/data/01_CytoSig_data/meta_info.txt",
  cytosig_effects = "/work/data/01_CytoSig_data/diff.merge.gz")
stopifnot(all(file.exists(paths)))
paths <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
meta <- as_tibble(fread(paths[["effect_metadata"]])) %>%
  mutate(compound = tolower(trimws(drug))) %>% filter(compound %in% drugs, available == 1L)
genes <- fread(paths[["all_genes"]]) %>% as_tibble() %>% arrange(gene_index)
stopifnot(nrow(genes) == 62710L, !anyDuplicated(genes$feature),
  identical(as.integer(genes$gene_index), seq_len(nrow(genes)) - 1L),
  "vemurafenib" %in% meta$compound)
folds <- as_tibble(fread(paths[["split"]])) %>%
  transmute(depmap_id, tahoe_cell_name, split = if_else(outer_fold == 1L, "Holdout", "Training"))
stopifnot(!anyDuplicated(folds$depmap_id), sum(folds$split == "Training") == 26L,
  sum(folds$split == "Holdout") == 7L)
prism <- as_tibble(fread(paths[["prism"]])) %>% filter(chain_observed == 1L) %>%
  mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"))
norm <- as_tibble(fread(paths[["normalization"]]))
stopifnot(!anyDuplicated(prism$observation_id), !anyDuplicated(norm$observation_id))
meta <- meta %>% left_join(folds, by = "depmap_id") %>%
  left_join(select(prism, observation_id, prism_sensitization, prism_dose_uM), by = "observation_id") %>%
  left_join(select(norm, observation_id, treated_total, reference_total,
    treated_correction_factor, contributing_plates, negative_tolerance), by = "observation_id") %>%
  arrange(compound, dose, depmap_id)
stopifnot(!anyNA(meta$split), !anyNA(meta$treated_total), !anyDuplicated(meta$observation_id))
## Adequacy is assessed from overlap counts only, before expression selection.
coverage <- meta %>% group_by(compound, dose, split) %>% summarise(
  n_molecular_lines = n_distinct(depmap_id),
  n_prism_lines = n_distinct(depmap_id[is.finite(prism_sensitization)]), .groups = "drop")
adequate <- coverage %>% filter(split == "Training") %>% group_by(compound) %>%
  summarise(include_in_discovery = any(n_prism_lines >= 8L), .groups = "drop")
stopifnot(adequate$include_in_discovery[adequate$compound == "vemurafenib"])
drugs <- adequate$compound[adequate$include_in_discovery]
meta <- meta %>% filter(compound %in% drugs)
fwrite(left_join(coverage, adequate, by = "compound"), file.path(out, "drug_overlap.tsv"), sep = "\t")
meta$matrix_row <- seq_len(nrow(meta))

## Read rows of saved processed pseudobulk arrays in place, using their original
## gene order. The same saved sample-level normalization factors apply to ALL
## genes; no raw single cells or upstream feature-selection code is rerun.
read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  stopifnot(identical(magic, "\x93NUMPY"), version[1] %in% 1:3)
  len <- readBin(con, "integer", 1L, size = if (version[1] == 1L) 2L else 4L,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", len))
  stopifnot(grepl("'descr': '<f4'", header, fixed = TRUE),
    grepl("'fortran_order': False", header, fixed = TRUE))
  shape <- as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*", "\\1", header), ",")[[1L]]))
  shape <- shape[!is.na(shape)]
  stopifnot(length(shape) == 2L, all(rows_zero_based >= 0 & rows_zero_based < shape[1]))
  offset <- seek(con)
  answer <- matrix(NA_real_, length(rows_zero_based), shape[2])
  for (i in seq_along(rows_zero_based)) {
    seek(con, offset + 4 * shape[2] * rows_zero_based[i], origin = "start")
    v <- readBin(con, "numeric", n = shape[2], size = 4L, endian = "little")
    stopifnot(length(v) == shape[2]); answer[i, ] <- v
  }
  answer
}
control_meta <- as_tibble(fread(paths[["control_metadata"]]))
plate_keys <- unlist(lapply(seq_len(nrow(meta)), function(i)
  paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1]], sep = "__")))
control_meta$key <- paste(control_meta$depmap_id, control_meta$plate, sep = "__")
control_meta <- control_meta %>% filter(key %in% plate_keys)
stopifnot(!anyDuplicated(control_meta$key), setequal(control_meta$key, plate_keys))
controls <- read_npy_rows(paths[["control_effects"]], control_meta$array_row)
rownames(controls) <- control_meta$key
effects <- matrix(NA_real_, nrow(meta), nrow(genes), dimnames = list(meta$observation_id, genes$feature))
abundance <- effects
for (path in unique(meta$shard_path)) {
  ids <- which(meta$shard_path == path)
  raw <- read_npy_rows(file.path(cache, path), meta$shard_row[ids])
  for (j in seq_along(ids)) {
    i <- ids[j]
    keys <- paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1]], sep = "__")
    reference <- colMeans(controls[keys, , drop = FALSE])
    treated <- raw[j, ] + reference
    stopifnot(all(is.finite(treated)), min(treated) >= -meta$negative_tolerance[i])
    treated <- pmax(treated, 0)
    abundance[i, ] <- pmax(treated / meta$treated_total[i] * 10000 / meta$treated_correction_factor[i],
      reference / meta$reference_total[i] * 10000)
    effects[i, ] <- log2(1 + treated / meta$treated_total[i] * 10000 /
      meta$treated_correction_factor[i]) - log2(1 + reference / meta$reference_total[i] * 10000)
  }
}
stopifnot(all(is.finite(effects)))
rm(controls, raw)
message("Reconstructed existing normalized pseudobulk effects for all ", ncol(effects),
  " genes, ", nrow(effects), " observations; no activity or old gene gate.")

## Discovery only: exact EGF profiles; never use Gerosa/GSE127988 here.
cyto_meta <- as_tibble(fread(paths[["cytosig_metadata"]], skip = 1L, header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose"))) %>%
  filter(toupper(trimws(treatment)) == "EGF") %>%
  mutate(study = sub(".*@(GSE[0-9]+|E-MTAB-[0-9]+).*", "\\1", profile_id))
stopifnot(nrow(cyto_meta) >= 3L, !any(grepl("GSE127988", cyto_meta$profile_id)),
  all(grepl("^(GSE|E-MTAB-)", cyto_meta$study)))
con <- gzfile(paths[["cytosig_effects"]], "rt")
profile_header <- strsplit(readLines(con, 1L), "\t", fixed = TRUE)[[1]]; close(con)
stopifnot(all(cyto_meta$profile_id %in% profile_header))
cyto_dt <- suppressWarnings(fread(paths[["cytosig_effects"]], header = TRUE,
  col.names = c("gene", cyto_meta$profile_id), check.names = FALSE,
  select = c(1L, match(cyto_meta$profile_id, profile_header) + 1L), showProgress = FALSE))
stopifnot(!anyDuplicated(cyto_dt$gene))
universe <- intersect(genes$feature, cyto_dt$gene)
cyto <- t(as.matrix(cyto_dt[match(universe, gene), -1L, with = FALSE]))
colnames(cyto) <- universe; rownames(cyto) <- cyto_meta$profile_id
cyto[!is.finite(cyto)] <- NA_real_
studies <- sort(unique(cyto_meta$study))
study_effect <- t(vapply(studies, function(s)
  apply(cyto[cyto_meta$study == s, , drop = FALSE], 2L, median, na.rm = TRUE), numeric(length(universe))))
colnames(study_effect) <- universe; rownames(study_effect) <- studies
egf <- tibble(gene = universe, n_egf_profiles = colSums(is.finite(cyto)),
  n_egf_studies = colSums(is.finite(study_effect)),
  egf_median_effect = apply(study_effect, 2L, median, na.rm = TRUE))

## New technical rule addresses the observed sparse-gene failure directly.
## This is gene detectability, not a drug-activity gate: all doses remain.
## 0.1 per 10,000 = 10 CPM, in at least half the matched training lines.
## A gene's PRISM association must retain its sign and nominal p<.05 after
## deleting ANY one training line. A failed/undefined deletion is a failure.
training <- which(meta$split == "Training")
gene_tests <- list()
for (drug in drugs) for (dose in sort(unique(meta$dose[meta$compound == drug]))) {
  ids <- which(meta$compound == drug & meta$dose == dose & meta$split == "Training" & is.finite(meta$prism_sensitization))
  x <- effects[ids, universe, drop = FALSE]; y <- meta$prism_sensitization[ids]; n <- length(y)
  xc <- sweep(x, 2, colMeans(x)); yc <- y - mean(y)
  ssx <- colSums(xc^2); ssy <- sum(yc^2); sxy <- colSums(xc * yc)
  r <- sxy / sqrt(ssx * ssy)
  p <- 2 * pt(abs(r) * sqrt((n-2)/pmax(0,1-r^2)), n-2, lower.tail = FALSE)
  loo <- matrix(NA_real_, n, ncol(x))
  for (i in seq_len(n)) {
    vx <- ssx - n/(n-1)*xc[i, ]^2; vy <- ssy - n/(n-1)*yc[i]^2
    loo[i, ] <- (sxy - n/(n-1)*xc[i, ]*yc[i]) / sqrt(pmax(vx,0)*pmax(vy,0))
    loo[i, vx <= ssx*1e-10 | !is.finite(loo[i, ])] <- NA_real_
  }
  loo[is.finite(loo)] <- pmax(-1,pmin(1,loo[is.finite(loo)]))
  stable <- colSums(!is.finite(loo) | sweep(sign(loo), 2, sign(r), "!=") , na.rm = TRUE) == 0
  worst_r <- apply(abs(loo),2,min); worst_r[!is.finite(worst_r)] <- 0
  worst_p <- 2 * pt(worst_r * sqrt((n-3)/pmax(0,1-worst_r^2)),n-3,lower.tail=FALSE)
  detected <- colSums(abundance[ids, universe, drop = FALSE] >= .1)
  dominance <- apply(xc^2,2,max)/ssx
  gene_tests[[length(gene_tests)+1L]] <- tibble(compound=drug,dose_uM=dose,gene=universe,
    n_training_lines=n,pearson_r=r,nominal_p=p,minimum_loo_abs_r=worst_r,
    maximum_loo_p=worst_p,loo_sign_stable=stable,n_detected_lines=detected,
    median_abundance_cpm=apply(abundance[ids,universe,drop=FALSE],2,median)*100,
    max_single_line_variance_fraction=dominance,
    eligible=is.finite(p)&p<.05&stable&worst_p<.05&detected>=max(8L,ceiling(n/2)))
}
gene_tests <- bind_rows(gene_tests)
drug_selected <- gene_tests %>% filter(eligible) %>% arrange(compound,desc(minimum_loo_abs_r),nominal_p,gene,dose_uM) %>%
  distinct(compound,gene,.keep_all=TRUE) %>% group_by(compound) %>% slice_head(n=100L) %>%
  mutate(selection_rank=row_number()) %>% ungroup()
## EGF candidates require consistent direction across independent study medians,
## plus measurable Tahoe abundance. Rank by the median effect, cap at 100.
egf$n_positive_studies <- colSums(study_effect>0,na.rm=TRUE)
egf$n_nonzero_studies <- colSums(is.finite(study_effect)&study_effect!=0)
egf$sign_test_p <- vapply(seq_len(nrow(egf)),function(i) {
  n <- egf$n_nonzero_studies[i]
  if(n<3) return(NA_real_)
  binom.test(egf$n_positive_studies[i],n,p=.5)$p.value
},numeric(1))
egf$detected_training_fraction <- colMeans(abundance[training,universe,drop=FALSE]>=.1)
egf_selected <- egf %>% filter(n_egf_studies>=3,sign_test_p<.05,detected_training_fraction>=.5) %>%
  arrange(desc(abs(egf_median_effect)),gene) %>% slice_head(n=100L) %>% mutate(selection_rank=row_number())
selected <- sort(union(drug_selected$gene,egf_selected$gene))
stopifnot(length(selected)>=3L,length(selected)<=300L)
## Balance CytoSig studies before correlation; repeated doses/time points from
## one study do not count as multiple independent source contexts.
r_tahoe <- cor(effects[training,selected,drop=FALSE])
r_cyto <- cor(study_effect[,selected,drop=FALSE],use="pairwise.complete.obs")
pair_n <- crossprod(is.finite(study_effect[,selected,drop=FALSE]))
r_cyto[pair_n<3] <- NA_real_
pair_audit <- tibble(gene=selected,observed_studies=colSums(is.finite(study_effect[,selected,drop=FALSE])),used=TRUE)
invalid <- !is.finite(r_tahoe)|!is.finite(r_cyto)
while(any(invalid)) {
 b <- order(-rowSums(invalid),colSums(is.finite(study_effect[,selected,drop=FALSE])),selected)[1]
 pair_audit$used[pair_audit$gene==selected[b]] <- FALSE
 selected <- selected[-b];r_tahoe <- r_tahoe[-b,-b,drop=FALSE];r_cyto <- r_cyto[-b,-b,drop=FALSE]
 invalid <- !is.finite(r_tahoe)|!is.finite(r_cyto)
}
fwrite(pair_audit,file.path(out,"candidate_correlation_coverage.tsv"),sep="\t")
matrices <- list(Joint=(r_tahoe+r_cyto)/2,Tahoe=r_tahoe,CytoSig_EGF=r_cyto)
tree <- hclust(as.dist(1-matrices$Joint),method="average")
## Cut at a stated correlation strength, not at a requested branch count.
## Average-linkage height <=0.6 corresponds to between-child joint r>=0.4.
## Groups smaller than eight genes stay unassigned instead of being absorbed.
cut_height <- .6
raw_groups <- cutree(tree,h=cut_height)
gene_order <- selected[tree$order]
group_sizes <- table(raw_groups)
supported_groups <- unique(raw_groups[gene_order])[group_sizes[as.character(unique(raw_groups[gene_order]))]>=8L]
program_levels <- if(length(supported_groups)) paste0("P",seq_along(supported_groups)) else character()
programs <- setNames(ifelse(raw_groups[gene_order] %in% supported_groups,
  paste0("P",match(raw_groups[gene_order],supported_groups)),"Unassigned"),gene_order)
membership <- tibble(gene=gene_order,program=unname(programs),gene_order=seq_along(gene_order),
  selected_vemurafenib=gene %in% drug_selected$gene[drug_selected$compound=="vemurafenib"],
  selected_cobimetinib=gene %in% drug_selected$gene[drug_selected$compound=="cobimetinib"],
  selected_egf=gene %in% egf_selected$gene) %>% left_join(egf,by="gene")
fwrite(gene_tests,file.path(out,"robust_gene_prism_training_tests.tsv.gz"),sep="\t",na="NA")
fwrite(drug_selected,file.path(out,"selected_genes_by_drug.tsv"),sep="\t",na="NA")
fwrite(egf_selected,file.path(out,"selected_egf_genes.tsv"),sep="\t",na="NA")
fwrite(membership,file.path(out,"candidate_gene_membership.tsv"),sep="\t",na="NA")
fwrite(cyto_meta,file.path(out,"cytosig_egf_profiles.tsv"),sep="\t")
for(s in names(matrices)) fwrite(data.table(gene=gene_order,matrices[[s]][gene_order,gene_order],check.names=FALSE),
  file.path(out,paste0("candidate_",tolower(s),"_correlations.tsv.gz")),sep="\t")
fwrite(data.table(gene=gene_order,t(effects[training,gene_order]),check.names=FALSE),
  file.path(out,"candidate_tahoe_gene_sample_effects.tsv.gz"),sep="\t")
fwrite(data.table(gene=gene_order,t(study_effect[,gene_order]),check.names=FALSE),
  file.path(out,"candidate_cytosig_study_effects.tsv.gz"),sep="\t",na="NA")
message("Discovery candidates and supported modules:")
print(drug_selected %>% select(compound,gene,dose_uM,minimum_loo_abs_r,median_abundance_cpm))
print(egf_selected %>% select(gene,egf_median_effect))
print(membership %>% dplyr::count(program))

## Annotate every candidate from HGNC, preserving original symbols and the
## exact source used to justify inclusion. Annotation is not outcome selection.
hgnc_path <- file.path(root,"artifacts/23_raf_egf_repair/hgnc_complete_set.tsv")
if(!file.exists(hgnc_path)) download.file(
  "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",hgnc_path,mode="wb",quiet=TRUE)
hgnc <- fread(hgnc_path,na.strings="")
resolve_hgnc <- function(g) {
  exact <- match(g,hgnc$symbol)
  if(!is.na(exact)) return(exact)
  hits <- which(vapply(seq_len(nrow(hgnc)),function(i)
    g %in% c(strsplit(ifelse(is.na(hgnc$prev_symbol[i]),"",hgnc$prev_symbol[i]),"|",fixed=TRUE)[[1]],
      strsplit(ifelse(is.na(hgnc$alias_symbol[i]),"",hgnc$alias_symbol[i]),"|",fixed=TRUE)[[1]]),logical(1)))
  if(length(hits)==1L) hits else NA_integer_
}
hi <- vapply(gene_order,resolve_hgnc,integer(1))
membership$approved_symbol <- hgnc$symbol[hi]
membership$gene_name <- hgnc$name[hi]
membership$locus_type <- hgnc$locus_type[hi]
membership$ensembl_gene_id <- hgnc$ensembl_gene_id[hi]
membership$hgnc_id <- hgnc$hgnc_id[hi]
membership$inclusion_evidence <- vapply(gene_order,function(g) {
  a <- drug_selected %>% filter(gene==g)
  b <- if(g %in% egf_selected$gene) "Consistent EGF response across study medians" else character()
  paste(c(if(nrow(a)) paste0(a$compound,": detectable across training lines; all single-line deletions retain nominal p<.05 at ",a$dose_uM," uM"),b),collapse="; ")
},character(1))
membership$program_decision <- ifelse(membership$program=="Unassigned",
  "Not retained in a >=8-gene module at joint-tree height 0.6", "Retained coherent module")
## Keep all candidate decisions auditable. The requested 30-60 gene scale is
## implemented as >=8 genes per module at the existing 0.6 height cut. This
## discovery-only size choice yields 52 genes; no GSE values enter this decision.
## Prune excluded leaves without reclustering or changing retained merge heights.
fwrite(membership,file.path(out,"candidate_gene_inclusion_evidence.tsv"),sep="\t",na="NA")
candidate_membership <- membership
candidate_tree <- tree
candidate_matrices <- matrices
candidate_selected <- selected
core_genes <- gene_order[programs!="Unassigned"]
stopifnot(length(core_genes)>=30L,length(core_genes)<=60L)
prune_dendrogram <- function(d) {
  if(is.leaf(d)) return(if(attr(d,"label") %in% core_genes) d else NULL)
  children <- Filter(Negate(is.null),lapply(d,prune_dendrogram))
  if(!length(children)) return(NULL)
  if(length(children)==1L) return(children[[1L]])
  a <- attributes(d); d <- children; attributes(d) <- a
  attr(d,"members") <- sum(vapply(d,function(x) attr(x,"members"),numeric(1)))
  d
}
pruned <- prune_dendrogram(as.dendrogram(candidate_tree))
leaf_index <- 0L
pruned <- dendrapply(pruned,function(d) {
  if(is.leaf(d)) {a <- attributes(d);leaf_index <<- leaf_index+1L;d <- leaf_index;attributes(d) <- a};d
})
tree <- as.hclust(pruned)
stopifnot(identical(tree$labels[tree$order],core_genes))
gene_order <- core_genes
selected <- tree$labels
programs <- programs[gene_order]
membership <- membership %>% filter(gene %in% core_genes) %>% mutate(gene_order=row_number())
stopifnot(identical(membership$gene,core_genes))
hi <- hi[candidate_membership$gene %in% gene_order]
matrices <- lapply(matrices,function(m) m[selected,selected])
fwrite(tibble(parameter=c("candidate_cap_per_drug","candidate_cap_egf","tree_cut_height","minimum_program_genes","target_total_min","target_total_max"),
  value=c(100,100,cut_height,8,30,60),basis=c(rep("Earlier candidate cap",2),"Existing coherence cut", "Size choice after user requested 30-60 total; discovery only",rep("User preference",2))),
  file.path(out,"selection_parameters.tsv"),sep="\t")
for(s in names(matrices)) fwrite(data.table(gene=gene_order,matrices[[s]][gene_order,gene_order],check.names=FALSE),
  file.path(out,paste0(tolower(s),"_correlations.tsv.gz")),sep="\t")
fwrite(data.table(gene=gene_order,t(effects[training,gene_order]),check.names=FALSE),file.path(out,"tahoe_gene_sample_effects.tsv.gz"),sep="\t")
fwrite(data.table(gene=gene_order,t(study_effect[,gene_order]),check.names=FALSE),file.path(out,"cytosig_study_effects.tsv.gz"),sep="\t",na="NA")

## Raw expression orientation is fixed at +1, matching the gene heatmaps.
## Association with PRISM is measured, never used to reverse a program's color.
membership$frozen_sign <- 1
membership$frozen_weight <- vapply(seq_len(nrow(membership)),function(i)
  if(membership$program[i]=="Unassigned") 0 else 1/sum(membership$program==membership$program[i]),numeric(1))
palette <- c(setNames(grDevices::hcl.colors(max(1,length(program_levels)),"Dark 3"),
  if(length(program_levels)) program_levels else "No module"),Unassigned="#888888")
weights <- matrix(0,length(gene_order),length(program_levels),dimnames=list(gene_order,program_levels))
for(p in program_levels) weights[programs==p,p] <- membership$frozen_weight[programs==p]
stopifnot(all(is.finite(weights)),all(abs(colSums(weights)-1)<1e-12))
coherence <- bind_rows(lapply(program_levels,function(p) {
  g <- gene_order[programs==p]
  within <- vapply(matrices,function(m) mean(m[g,g][upper.tri(m[g,g])]),numeric(1))
  tibble(program=p,n_genes=length(g),joint_mean_r=within[["Joint"]],tahoe_mean_r=within[["Tahoe"]],
    cytosig_egf_mean_r=within[["CytoSig_EGF"]],egf_effect=mean(egf$egf_median_effect[match(g,egf$gene)]),
    vem_effect=mean(effects[meta$split=="Training"&meta$compound=="vemurafenib",g]),
    cobi_effect=if("cobimetinib" %in% drugs) mean(effects[meta$split=="Training"&meta$compound=="cobimetinib",g]) else NA_real_,
    n_vem_genes=sum(membership$selected_vemurafenib[programs==p]),
    n_cobi_genes=sum(membership$selected_cobimetinib[programs==p]),n_egf_genes=sum(membership$selected_egf[programs==p]),
    genes=paste(g,collapse=", "))
}))
if(length(program_levels)) coherence <- coherence %>% mutate(candidate_role=case_when(
  tahoe_mean_r>=.4 & cytosig_egf_mean_r>=.4 ~ "Coherent in both sources",
  tahoe_mean_r>=.4 ~ "Tahoe-supported core",cytosig_egf_mean_r>=.4 ~ "EGF-supported core",
  TRUE ~ "Joint core; source support weak"))
## Enrichment uses the actual 270-gene clustering pool as its background.
## These tests describe selection provenance; they do not infer causal edges.
enrichment <- bind_rows(lapply(program_levels,function(p) bind_rows(lapply(
  c("selected_vemurafenib","selected_cobimetinib","selected_egf"),function(source) {
    inside <- candidate_membership$program==p; positive <- candidate_membership[[source]]
    tab <- matrix(c(sum(inside&positive),sum(inside&!positive),sum(!inside&positive),sum(!inside&!positive)),2,byrow=TRUE)
    ft <- fisher.test(tab,alternative="greater")
    tibble(program=p,source=sub("selected_","",source),n_overlap=tab[1,1],n_program=sum(inside),
      n_source=sum(positive),n_background=nrow(candidate_membership),odds_ratio=unname(ft$estimate),nominal_p=ft$p.value)
  })))) %>% mutate(bh_q=p.adjust(nominal_p,"BH"))
fwrite(enrichment,file.path(out,"program_source_enrichment.tsv"),sep="\t",na="NA")
anchor_program <- function(g) unname(programs[g])
branch_assignment <- tibble(
  program=c(anchor_program("AREG"),anchor_program("AKAP12"),anchor_program("JUNB"),anchor_program("IFIT2")),
  label=c("AREG / ERRFI1","FOSL1 / AKAP12","EGR1 / JUNB","IFIT2 / OASL"),
  placement=c("Receptor feedback candidate","Shared response candidate","Shared response candidate","Pathway position unresolved"),
  evidence=c("Contains receptor ligand AREG and feedback regulator ERRFI1; EGF-induced; selected drug-associated genes present",
    "EGF-induced and coherently drug-responsive; one selected gene per drug; weaker Tahoe coherence",
    "Coherent in Tahoe and CytoSig EGF; all genes selected via EGF; no selected PRISM genes",
    "Nine genes selected through cobimetinib PRISM; no selected EGF genes; negligible mean EGF effect"),
  limitation=c("Cannot place the whole transcript module causally upstream of MEK/ERK",
    "Broad response readout; functional survival role not established for the whole program",
    "Shared molecular readout candidate; functional PRISM support must be tested at program level",
    "Cobimetinib association does not establish a BRAF-specific arm or an upstream position"),
  placement_basis="Gene content and discovery source evidence; hypothesis, not causal reconstruction",
  mechanism_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC8009031/")
stopifnot(setequal(branch_assignment$program,program_levels),!anyNA(branch_assignment$program))
fwrite(branch_assignment,file.path(out,"frozen_branch_assignments.tsv"),sep="\t")

predictions <- coherence %>% transmute(program,candidate_role,
  expected_egf_direction=case_when(n_egf_genes==0L~"exploratory",egf_effect>0~"increase",TRUE~"decrease"),
  expected_cobimetinib_direction=case_when(!is.finite(cobi_effect)~"exploratory",cobi_effect<0~"decrease",TRUE~"increase"),
  primary_window="Trapezoidal 0-8 h mean displacement; each biological day is one unit",
  sign_rule="Positive score means increased raw expression; all program weights positive",
  missing_gene_rule="Map exact symbols and unique HGNC stable IDs/history; no ambiguous mappings, replacements or refits",
  caution="GSE127988 was inspected previously; this revised fit is a retrospective molecular check")
fwrite(membership,file.path(out,"frozen_program_genes.tsv"),sep="\t",na="NA")
fwrite(coherence,file.path(out,"frozen_program_annotation.tsv"),sep="\t",na="NA")
fwrite(predictions,file.path(out,"frozen_predictions.tsv"),sep="\t",na="NA")
fwrite(membership,file.path(out,"gene_inclusion_evidence.tsv"),sep="\t",na="NA")
freeze_paths <- file.path(out,c("frozen_program_genes.tsv","frozen_program_annotation.tsv","frozen_predictions.tsv","frozen_branch_assignments.tsv"))
freeze_hashes <- tools::md5sum(freeze_paths)
fwrite(tibble(path=freeze_paths,md5=unname(freeze_hashes)),file.path(out,"freeze_manifest.tsv"),sep="\t")
annotation <- data.frame(Program=factor(programs,levels=c(program_levels,"Unassigned")),row.names=gene_order)

## A readable gene hierarchy, with a visible cut and explicit gene provenance.
tree_geometry <- function(hc) {
  n <- length(hc$order); ypos <- setNames(rev(seq_len(n)),hc$labels[hc$order]); nodes <- list(); seg <- list()
  for(i in seq_len(n-1L)) {
    children <- lapply(hc$merge[i,],function(j) if(j<0) c(x=0,y=ypos[hc$labels[-j]]) else nodes[[j]])
    yy <- vapply(children,function(a) unname(a[2]),numeric(1)); xx <- vapply(children,function(a) unname(a[1]),numeric(1))
    nodes[[i]] <- c(x=hc$height[i],y=mean(yy))
    seg[[i]] <- tibble(x=c(xx,hc$height[i]),xend=rep(hc$height[i],3),y=c(yy,min(yy)),yend=c(yy,max(yy)))
  }
  list(segments=bind_rows(seg),leaves=tibble(label=hc$labels[hc$order],y=rev(seq_len(n))))
}
tg <- tree_geometry(tree)
leaf <- tg$leaves %>% left_join(transmute(membership,label=gene,program,selected_vemurafenib,selected_cobimetinib,selected_egf),by="label")
p_tree <- ggplot() + geom_segment(data=tg$segments,aes(x,y,xend=xend,yend=yend),linewidth=.65,color="#506070") +
  geom_vline(xintercept=cut_height,linetype=2,color="#B44B35") +
  geom_point(data=leaf,aes(0,y,color=program),size=2.5) +
  geom_text(data=leaf,aes(-.03,y,label=label,color=program),hjust=0,size=4,show.legend=FALSE) +
  geom_text(data=leaf,aes(-.67,y,label=ifelse(selected_vemurafenib,"V","")),size=4,color="#286E92") +
  geom_text(data=leaf,aes(-.75,y,label=ifelse(selected_cobimetinib,"C","")),size=4,color="#A35935") +
  geom_text(data=leaf,aes(-.83,y,label=ifelse(selected_egf,"E","")),size=4,color="#7757A4") +
  scale_color_manual(values=palette) + scale_x_reverse(limits=c(max(tree$height)+.05,-.94),breaks=c(0,.2,.4,.6,.8,1)) +
  scale_y_continuous(NULL,breaks=NULL,expand=expansion(add=c(.5,1))) +
  labs(x="Joint correlation distance (1 - r)",color=NULL,
    title="Gene hierarchy | Frozen 52-gene program set",
    subtitle=paste0(length(gene_order)," genes; ",length(program_levels)," modules from the cut"),
    caption="Dashed line: average-linkage height 0.6. Whole groups require >=8 genes; excluded candidates are audited separately.\nV = robust vemurafenib-PRISM association; C = cobimetinib; E = replicated EGF response.") +
  theme_classic(base_size=13) + theme(axis.line.y=element_blank(),legend.position="bottom",
    plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
ggsave(file.path(out,"gene_hierarchy.pdf"),p_tree,width=10.4,height=max(6,length(gene_order)*.19+2),device=cairo_pdf)
ggsave(file.path(out,"gene_hierarchy.png"),p_tree,width=10.4,height=max(6,length(gene_order)*.19+2),dpi=150)

write_matrix <- function(m,filename) fwrite(data.table(gene=rownames(m),m,check.names=FALSE),
  file.path(out,filename),sep="\t",na="NA")
save_heatmap <- function(mat,filename,title,range,legend,column_annotation=NULL,row_annotation=NULL,
    labelled=TRUE,column_labels=colnames(mat),width=10,height=8,row_labels=rownames(mat),
    row_tree=FALSE,column_tree=FALSE) {
  if(!is.null(row_annotation)) stopifnot(identical(rownames(mat),rownames(row_annotation)))
  if(!is.null(column_annotation)) stopifnot(identical(colnames(mat),rownames(column_annotation)))
  cairo_pdf(file.path(out,filename),width=width,height=height)
  ht <- draw_heatmap(mat,scale_rows=FALSE,cluster_rows=row_tree,cluster_columns=column_tree,
    row_dend_reorder=FALSE,column_dend_reorder=FALSE,
    row_dend_width=grid::unit(22,"mm"),column_dend_height=grid::unit(22,"mm"),
    row_anno_df=row_annotation,row_anno_cols=names(row_annotation),
    col_anno_df=column_annotation,col_anno_cols=names(column_annotation),
    palette_overrides=list(Program=palette,Arm=c(Vem="#247A96",`Vem + Cobi`="#B95B32",DMSO="#777777")),
    ann_row_args=list(show_legend=FALSE),hm_name=legend,value_range=range,
    value_colors=c("#2166AC","white","#B2182B"),show_row_names=labelled,show_column_names=!is.null(column_labels),
    row_labels=row_labels,column_labels=column_labels,row_names_gp=grid::gpar(fontsize=12),
    column_names_gp=grid::gpar(fontsize=10),column_title_gp=grid::gpar(fontsize=15,fontface="bold"),
    title=title,na_col="#BDBDBD",heatmap_legend_param=list(labels_gp=grid::gpar(fontsize=11),title_gp=grid::gpar(fontsize=11,fontface="bold")))
  extra <- if(anyNA(mat)) list(ComplexHeatmap::Legend(labels="Unavailable",legend_gp=grid::gpar(fill="#BDBDBD"))) else list()
  drawn <- ComplexHeatmap::draw(ht,merge_legends=TRUE,annotation_legend_list=extra,padding=grid::unit(c(5,5,5,5),"mm"))
  stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn)),if(inherits(row_tree,"hclust")) row_tree$order else seq_len(nrow(mat))),
    identical(as.integer(ComplexHeatmap::column_order(drawn)),if(inherits(column_tree,"hclust")) column_tree$order else seq_len(ncol(mat))))
  dev.off()
}
for(s in names(matrices)) {
  a <- annotation[match(selected,rownames(annotation)),,drop=FALSE]
  save_heatmap(matrices[[s]][selected,selected],paste0(tolower(s),"_correlation.pdf"),
    paste0(s," | Same joint gene tree in every source"),c(-1,0,1),"Pearson r",a,a,
    column_labels=selected,row_labels=selected,row_tree=tree,column_tree=tree,width=13,height=13)
}
core_genes <- gene_order[programs!="Unassigned"]
core_annotation <- annotation[core_genes,,drop=FALSE]
core_split <- factor(programs[core_genes],levels=program_levels)
core_maps <- lapply(seq_along(matrices),function(i) {
  s <- names(matrices)[i];m <- matrices[[i]][core_genes,core_genes]
  write_matrix(m,paste0("core_",tolower(s),"_correlations.tsv.gz"))
  save_heatmap(m,paste0("core_",tolower(s),"_correlation.pdf"),paste0(s," | Frozen ",length(core_genes),"-gene programs"),
    c(-1,0,1),"Pearson r",core_annotation,core_annotation,width=12,height=12)
  draw_heatmap(m,scale_rows=FALSE,cluster_rows=FALSE,cluster_columns=FALSE,
    row_order=seq_along(core_genes),column_order=seq_along(core_genes),
    row_anno_df=if(i==1) core_annotation else NULL,row_anno_cols=if(i==1) "Program" else NULL,
    palette_overrides=list(Program=palette),ann_row_args=list(show_legend=FALSE),
    show_row_names=i==1,row_names_side="left",show_column_names=TRUE,
    row_names_gp=grid::gpar(fontsize=11),column_names_gp=grid::gpar(fontsize=9),
    row_title_gp=grid::gpar(fontsize=12,fontface="bold"),column_title_gp=grid::gpar(fontsize=15,fontface="bold"),
    hm_name=paste0("r",i),title=s,value_range=c(-1,0,1),value_colors=c("#2166AC","white","#B2182B"),
    show_heatmap_legend=i==3,heatmap_legend_param=list(title="Pearson r",labels_gp=grid::gpar(fontsize=11)),
    row_gap=grid::unit(2,"mm"),column_gap=grid::unit(2,"mm"))
})
core_comparison <- core_maps[[1]]+core_maps[[2]]+core_maps[[3]]
draw_core_comparison <- function() {
  drawn <- ComplexHeatmap::draw(core_comparison,merge_legends=TRUE,padding=grid::unit(c(5,5,5,5),"mm"))
  for(i in 1:3) {
    name <- paste0("r",i)
    stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn,name=name)),seq_along(core_genes)),
      identical(as.integer(ComplexHeatmap::column_order(drawn,name=name)),seq_along(core_genes)))
    ComplexHeatmap::decorate_heatmap_body(name,{
      for(b in head(cumsum(as.integer(table(core_split))),-1)) {
        grid::grid.lines(c(0,1),rep(1-b/length(core_genes),2),gp=grid::gpar(col="white",lwd=2))
        grid::grid.lines(rep(b/length(core_genes),2),c(0,1),gp=grid::gpar(col="white",lwd=2))
      }
    })
  }
}
cairo_pdf(file.path(out,"core_correlation_comparison.pdf"),width=24,height=10.4)
draw_core_comparison();dev.off()
png(file.path(out,"core_correlation_comparison.png"),width=3600,height=1560,res=150,type="cairo")
draw_core_comparison();dev.off()
curve <- bind_rows(lapply(2:min(8,length(gene_order)-1L),function(k) {
  groups <- cutree(tree,k)
  tibble(rank=k,mean_silhouette=mean(cluster::silhouette(groups,as.dist(1-matrices$Joint))[,"sil_width"]),smallest_group=min(table(groups)))
}))
fwrite(curve,file.path(out,"rank_silhouette.tsv"),sep="\t")
pr <- ggplot(curve,aes(rank,mean_silhouette))+geom_line(color="#247A96")+geom_point(size=2.5)+
  labs(x="Number of clusters (including small groups)",y="Mean silhouette",title="Cluster-count diagnostic",
    subtitle="Final cores use correlation strength and minimum size; rank is not forced")+theme_classic(base_size=13)
ggsave(file.path(out,"rank_silhouette.pdf"),pr,width=7,height=4.7,device=cairo_pdf)
for(s in c("Tahoe","CytoSig_EGF")) {
  x <- if(s=="Tahoe") t(effects[training,gene_order]) else t(study_effect[,gene_order])
  z <- sweep(sweep(x,1,rowMeans(x,na.rm=TRUE)),1,apply(x,1,sd,na.rm=TRUE),"/")
  a <- t(vapply(program_levels,function(p) colMeans(z[programs==p,,drop=FALSE],na.rm=TRUE),numeric(ncol(z))))
  group <- max.col(t(a),ties.method="first")
  ord <- unlist(lapply(seq_along(program_levels),function(i) {
    ids <- which(group==i)
    if(length(ids)<3) ids else ids[hclust(dist(t(a[,ids,drop=FALSE])),method="average")$order]
  }))
  ca <- data.frame(Program=factor(program_levels[group[ord]],levels=program_levels),row.names=colnames(z)[ord])
  if(s=="Tahoe") {
    mm <- meta[match(colnames(z)[ord],meta$observation_id),];ca$Drug <- mm$compound;ca$Dose <- factor(mm$dose)
  }
  save_heatmap(z[,ord],paste0(tolower(s),"_gene_samples.pdf"),paste0(s," | Fixed gene order; samples grouped by core"),
    c(-2,0,2),"Gene z-score",ca,annotation,column_labels=if(s=="Tahoe") NULL else colnames(z)[ord],width=12.5,height=12.5)
}
## Program hierarchy uses between-core correlations, annotated by source support
## and drug/EGF gene provenance. It is a descriptive relationship, not ancestry.
if(length(program_levels)>=2L) {
  between <- matrix(0,length(program_levels),length(program_levels),dimnames=list(program_levels,program_levels))
  for(p in program_levels) for(q in program_levels) between[p,q] <-
    if(p==q) 1 else mean(matrices$Joint[names(programs)[programs==p],names(programs)[programs==q]])
  ptree <- hclust(as.dist(1-between),method="average");pg <- tree_geometry(ptree)
  tips <- pg$leaves %>% left_join(coherence,by=c("label"="program"))
  pp <- ggplot()+geom_segment(data=pg$segments,aes(x,y,xend=xend,yend=yend),linewidth=.8,color="#506070")+
    geom_point(data=tips,aes(0,y,color=label),size=3)+
    geom_text(data=tips,aes(-.04,y,label=paste0(label," | ",n_genes," genes\n",candidate_role," | V=",n_vem_genes,", C=",n_cobi_genes,", E=",n_egf_genes)),
      hjust=0,size=4)+scale_color_manual(values=palette,guide="none")+
    scale_x_reverse(limits=c(max(ptree$height)+.05,-1.8),breaks=c(0,.2,.4,.6,.8,1))+
    scale_y_continuous(NULL,breaks=NULL,expand=expansion(add=.6))+
    labs(x="Between-program joint correlation distance",title="Program hierarchy | Core relationships",
      caption="Branches use average between-core correlation. Source support requires within-core mean r >= 0.4.\nV/C/E count selection provenance; they do not establish causal signaling direction.")+
    theme_classic(base_size=13)+theme(axis.line.y=element_blank(),plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
  ggsave(file.path(out,"program_hierarchy.pdf"),pp,width=12,height=max(4,length(program_levels)*.8+1.6),device=cairo_pdf)
  ggsave(file.path(out,"program_hierarchy.png"),pp,width=12,height=max(4,length(program_levels)*.8+1.6),dpi=150)
}
pearson_test <- function(x,y) {
  ok <- is.finite(x)&is.finite(y);n <- sum(ok)
  if(n<3||sd(x[ok])==0||sd(y[ok])==0) return(tibble(n=n,pearson_r=NA_real_,nominal_p=NA_real_))
  a <- cor.test(x[ok],y[ok]);tibble(n=n,pearson_r=unname(a$estimate),nominal_p=a$p.value)
}
score_matrix <- effects[,gene_order] %*% weights
scores <- bind_rows(lapply(program_levels,function(p) meta %>% transmute(observation_id,depmap_id,
  cell_line=tahoe_cell_name,compound,dose_uM=dose,split,prism_sensitization,program=p,score=score_matrix[,p])))
associations <- scores %>% group_by(compound,dose_uM,split,program) %>%
  group_modify(~pearson_test(.x$score,-.x$prism_sensitization)) %>% ungroup() %>%
  group_by(split) %>% mutate(bh_q=p.adjust(nominal_p,"BH")) %>% ungroup()
fwrite(scores,file.path(out,"tahoe_program_scores.tsv"),sep="\t")
fwrite(associations,file.path(out,"program_prism_associations.tsv"),sep="\t",na="NA")
cairo_pdf(file.path(out,"program_prism_scatterplots.pdf"),width=11.6,height=8.4,onefile=TRUE)
for(p in program_levels) {
  d <- scores %>% filter(program==p,is.finite(prism_sensitization))
  a <- associations %>% filter(program==p) %>% mutate(label=sprintf("%s n=%d r=%.2f p=%.2g",split,n,pearson_r,nominal_p),v=if_else(split=="Training",1.2,2.6))
  print(ggplot(d,aes(score,-prism_sensitization,color=split))+geom_point(size=1.8)+
    geom_smooth(method="lm",formula=y~x,se=FALSE,linewidth=.7)+
    geom_text(data=a,aes(x=-Inf,y=Inf,label=label,vjust=v),inherit.aes=FALSE,hjust=-.02,size=3)+
    facet_grid(compound~dose_uM,scales="free")+scale_color_manual(values=c(Training="#247A96",Holdout="#B95B32"))+
    labs(x="Mean gene-expression change (positive = increased expression)",y="PRISM resistance (-sensitization)",
      title=paste0(p," | ",coherence$n_genes[coherence$program==p]," frozen genes"),subtitle="Dose labels in uM; training p values reuse discovery data")+
    theme_classic(base_size=12)+theme(aspect.ratio=1,legend.position="bottom"))
}
dev.off()
## Mechanistic Y: pathway direction comes from Gerosa et al.; the expression
## modules are candidate readouts, with no inferred program-to-program arrows.
render_mechanistic_y <- function(coherence,branches,filename) {
  ink <- "#26394A"; blue <- "#287A91"; orange <- "#B96535"; gray <- "#6D7680"
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale=c(0,15),yscale=c(0,12)))
  txt <- function(x,y,label,size=12,col=ink,font=1,just="centre")
    grid::grid.text(label,x=grid::unit(x,"native"),y=grid::unit(y,"native"),just=just,
      gp=grid::gpar(fontsize=size,col=col,fontface=font,lineheight=1.15))
  line <- function(x,y,xe,ye,color=ink,dashed=FALSE,arrow=FALSE,lwd=2)
    grid::grid.lines(x=grid::unit(c(x,xe),"native"),y=grid::unit(c(y,ye),"native"),
      arrow=if(arrow) grid::arrow(length=grid::unit(2.5,"mm"),type="closed") else NULL,
      gp=grid::gpar(col=color,lty=if(dashed) 3 else 1,lwd=lwd))
  box <- function(x,y,w,h,label,fill="#F2F5F7",border=ink,size=14) {
    grid::grid.roundrect(x=grid::unit(x,"native"),y=grid::unit(y,"native"),width=grid::unit(w,"native"),
      height=grid::unit(h,"native"),r=grid::unit(2,"mm"),gp=grid::gpar(fill=fill,col=border,lwd=1.4))
    txt(x,y,label,size,border,2)
  }
  card <- function(g,x,y,w,h,header) {
    p <- branches$program[branches$label==g]; c <- coherence[coherence$program==p,]
    color <- palette[[p]]
    box(x,y,w,h,"",fill=adjustcolor(color,alpha.f=.065),border=color)
    txt(x,y+h/2-.28,paste0(p," | ",header," | ",c$n_genes," genes"),13,color,2)
    txt(x,y+.28,g,12,ink)
    txt(x,y-.08,sprintf("Within-program r: Tahoe %.2f  /  EGF %.2f",c$tahoe_mean_r,c$cytosig_egf_mean_r),11.5)
    txt(x,y-.43,paste0("Selected genes: V ",c$n_vem_genes,"  /  C ",c$n_cobi_genes,"  /  E ",c$n_egf_genes),11.5)
    txt(x,y-h/2+.22,branches$placement[branches$program==p],11,gray)
  }
  txt(.45,11.65,"RAF–EGF convergence | 52 genes in four frozen programs",21,ink,2,"left")
  txt(.45,11.2,"Solid arrows: established signaling. Dotted links: candidate placement of expression readouts.",12.5,gray,just="left")
  box(3.4,10.3,3.2,.72,"EGF / EGFR",fill="#E9F2F6",border=blue,size=16)
  box(3.4,8.85,3.6,.75,"RAS-GTP / RAF dimers",fill="#E9F2F6",border=blue,size=15)
  line(3.4,9.94,3.4,9.24,blue,arrow=TRUE)
  box(11.25,10.3,3.2,.72,"BRAF V600E",fill="#FCF0E8",border=orange,size=16)
  box(11.8,8.95,2.6,.55,"Vemurafenib",fill="white",border=orange,size=13)
  line(11.8,9.23,11.8,9.89,orange);line(11.55,9.89,12.05,9.89,orange,lwd=3)
  # The two upstream inputs visibly converge onto the common MEK node.
  line(3.4,8.47,6.4,7.39,blue,arrow=TRUE,lwd=2.7)
  line(10.6,9.94,8.1,7.39,orange,arrow=TRUE,lwd=2.7)
  box(7.25,7.0,2.8,.75,"MEK",fill="#EEF2EA",size=17)
  box(7.25,5.6,2.8,.75,"ERK",fill="#EEF2EA",size=17)
  line(7.25,6.62,7.25,5.98,arrow=TRUE,lwd=2.7)
  box(11.8,7.0,2.6,.55,"Cobimetinib",fill="white",border=orange,size=13)
  line(10.5,7.0,8.68,7.0,orange);line(8.68,6.77,8.68,7.23,orange,lwd=3)
  # A receptor-associated readout; the dotted line is explicitly not causation.
  line(1.8,10.3,.45,10.3,blue,dashed=TRUE,lwd=1.4)
  line(.45,10.3,.45,6.5,blue,dashed=TRUE,lwd=1.4)
  line(.45,6.5,.7,6.5,blue,dashed=TRUE,lwd=1.4)
  card("AREG / ERRFI1",3.0,6.5,4.6,2.05,"Receptor feedback")
  card("IFIT2 / OASL",11.8,4.65,5.35,2.05,"Drug-associated")
  # Both response modules are downstream readout candidates, not parent/child.
  line(7.25,5.22,7.25,4.2,dashed=TRUE,lwd=1.4)
  line(7.25,4.2,3.85,3.16,dashed=TRUE,lwd=1.4)
  line(7.25,4.2,10.15,3.16,dashed=TRUE,lwd=1.4)
  card("EGR1 / JUNB",3.85,2.1,5.9,2.05,"Shared response")
  card("FOSL1 / AKAP12",10.15,2.1,5.9,2.05,"Shared response")
  txt(.45,.61,"V / C / E = selection through vemurafenib–PRISM / cobimetinib–PRISM / CytoSig EGF. Counts may overlap.",11,gray,just="left")
  txt(.45,.29,"Mechanism: Gerosa et al., Cell Systems 2020, doi:10.1016/j.cels.2020.10.002. Program placements remain hypotheses.",11,gray,just="left")
  grid::popViewport()
}
for(format in c("pdf","png","svg")) {
  destination <- file.path(out,paste0("raf_egf_mechanistic_y.",format))
  switch(format,pdf=cairo_pdf(destination,width=15,height=12),
    png=png(destination,width=2250,height=1800,res=150,type="cairo"),
    svg=svg(destination,width=15,height=12))
  render_mechanistic_y(coherence,branch_assignment,destination);dev.off()
}

if("--discovery-only" %in% commandArgs(trailingOnly=TRUE)) quit(save="no")
## HELD-OUT EXPRESSION IS FIRST READ BELOW THIS LINE.
geo_dir <- file.path(root, "artifacts/22_raf_egf_bypass/gse127988")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
geo_urls <- c(
  GSE127988_family.soft.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/soft/GSE127988_family.soft.gz",
  GSE127988_RPM_Gerosa_et_al.csv.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/suppl/GSE127988_RPM_Gerosa_et_al.csv.gz")
for (n in names(geo_urls)) if (!file.exists(file.path(geo_dir, n)))
  download.file(geo_urls[[n]], file.path(geo_dir, n), mode = "wb", quiet = TRUE)
fwrite(tibble(file = names(geo_urls), url = unname(geo_urls),
  md5 = unname(tools::md5sum(file.path(geo_dir, names(geo_urls))))),
  file.path(out, "external_input_manifest.tsv"), sep = "\t")
soft <- readLines(gzfile(file.path(geo_dir, names(geo_urls)[1])), warn = FALSE)
starts <- which(startsWith(soft, "^SAMPLE = "))
geo_meta <- bind_rows(lapply(seq_along(starts), function(i) {
  block <- soft[starts[i]:(if (i == length(starts)) length(soft) else starts[i+1] - 1L)]
  field <- function(prefix) {
    hit <- block[startsWith(block, prefix)]
    stopifnot(length(hit) == 1L)
    substring(hit, nchar(prefix) + 1L)
  }
  title <- field("!Sample_title = ")
  tibble(gsm = field("^SAMPLE = "), title = title,
    sample_number = as.integer(sub(".*\\(S([0-9]+)\\)$", "\\1", title)),
    vem_uM = as.numeric(field("!Sample_characteristics_ch1 = vemurafenib for 24 hrs (um): ")),
    cobi_uM = as.numeric(field("!Sample_characteristics_ch1 = cobimetinib for 24 hrs (um): ")),
    egf_h = as.numeric(field("!Sample_characteristics_ch1 = egf at 100 ng/ml (h): ")),
    biological_set = as.integer(field("!Sample_characteristics_ch1 = set id: ")))
})) %>% arrange(sample_number) %>% mutate(sample_id = paste0("S", sample_number),
  experiment = if_else(biological_set %in% 1:2, "Vemurafenib dose", "EGF time course"),
  arm = case_when(vem_uM == 0 ~ "DMSO", cobi_uM > 0 ~ "Vem + Cobi", TRUE ~ "Vem"))
stopifnot(nrow(geo_meta) == 96L, identical(geo_meta$sample_number, 1:96),
  setequal(geo_meta$biological_set, 1:4))
csv_path <- file.path(geo_dir, names(geo_urls)[2])
headers <- fread(csv_path, nrows = 6L, header = FALSE, colClasses = "character", check.names = FALSE)
individual <- which(as.character(headers[4, ]) == "Individual")
sample_numbers <- as.integer(as.character(headers[5, ..individual]))
stopifnot(identical(sample_numbers, 1:96),
  all(as.numeric(as.character(headers[1, ..individual])) == geo_meta$vem_uM),
  all(as.numeric(as.character(headers[2, ..individual])) == geo_meta$cobi_uM),
  all(as.numeric(as.character(headers[3, ..individual])) == geo_meta$egf_h))
geo_meta$geo_reported_set <- geo_meta$biological_set
geo_meta$csv_reported_set <- as.integer(as.character(headers[6, ..individual]))
discrepancies <- geo_meta %>% filter(geo_reported_set != csv_reported_set)
## The publisher matrix assigns S15/S16 DMSO wells to dose-experiment set 1;
## GEO SOFT assigns them to set 2. Use the matrix's explicit set labels, audit
## both, and require complete agreement for the held-out time-course experiment.
stopifnot(all(discrepancies$sample_number %in% c(15L, 16L)),
  all(discrepancies$geo_reported_set == 2L), all(discrepancies$csv_reported_set == 1L))
fwrite(discrepancies %>% mutate(resolution = "Use published expression-matrix set ID; dose-series DMSO only"),
  file.path(out, "gse127988_metadata_discrepancies.tsv"), sep = "\t")
geo_meta$biological_set <- geo_meta$csv_reported_set
raw_rpm <- fread(csv_path, skip=7L, header=FALSE, select=c(1L,2L,individual),
  col.names=c("external_ensembl","external_symbol",geo_meta$sample_id),check.names=FALSE)
stopifnot(nrow(raw_rpm)==18261L,all(as.matrix(raw_rpm[,-(1:2),with=FALSE])>=0))
map_rows <- lapply(seq_along(gene_order),function(i) {
  g <- gene_order[i]; h <- hi[i]
  exact <- which(raw_rpm$external_symbol==g)
  stable <- if(!is.na(h)) which(sub("\\..*","",raw_rpm$external_ensembl) %in%
    strsplit(ifelse(is.na(hgnc$ensembl_gene_id[h]),"",hgnc$ensembl_gene_id[h]),"|",fixed=TRUE)[[1]]) else integer()
  historical <- if(!is.na(h)) which(raw_rpm$external_symbol %in% c(hgnc$symbol[h],
    strsplit(ifelse(is.na(hgnc$prev_symbol[h]),"",hgnc$prev_symbol[h]),"|",fixed=TRUE)[[1]])) else integer()
  ids <- if(length(exact)) exact else if(length(stable)==1L) stable else if(length(historical)==1L) historical else integer()
  list(ids=ids,audit=tibble(gene=g,status=if(length(exact)) "Exact symbol" else if(length(stable)==1L) "HGNC Ensembl ID" else if(length(historical)==1L) "Unique previous symbol" else "Absent or ambiguous in published matrix",
    external_symbol=paste(raw_rpm$external_symbol[ids],collapse=";"),
    external_ensembl=paste(raw_rpm$external_ensembl[ids],collapse=";"),ensembl_rows=length(ids)))
})
mapping <- bind_rows(lapply(map_rows,`[[`,"audit"))
fwrite(mapping,file.path(out,"gse127988_identifier_audit.tsv"),sep="\t")
## Preserve frozen gene identities. Resolve only exact or unambiguous HGNC
## identifiers; no outcome-dependent mapping or replacement genes.
rpm <- rbindlist(lapply(seq_along(map_rows),function(i) {
  ids <- map_rows[[i]]$ids
  if(!length(ids)) return(NULL)
  cbind(data.table(gene=gene_order[i]),as.data.table(as.list(colSums(as.matrix(raw_rpm[ids,-(1:2),with=FALSE])))))
}),use.names=TRUE)
available_genes <- gene_order[gene_order %in% rpm$gene]
external_weights <- weights[available_genes, , drop = FALSE]
coverage_weight <- colSums(abs(external_weights))
external_weights <- sweep(external_weights, 2L, ifelse(coverage_weight > 0, coverage_weight, 1), "/")
gene_coverage <- membership %>% mutate(measured_gse127988 = gene %in% available_genes) %>%
  left_join(mapping, by = "gene")
program_coverage <- gene_coverage %>% group_by(program) %>% summarise(
  n_frozen = n(), n_measured = sum(measured_gse127988), fraction_measured = mean(measured_gse127988), .groups = "drop") %>% filter(program %in% program_levels)
fwrite(gene_coverage, file.path(out, "gse127988_gene_coverage.tsv"), sep = "\t", na = "NA")
fwrite(program_coverage, file.path(out, "gse127988_program_coverage.tsv"), sep = "\t")
gene_expression <- log2(1 + as.matrix(rpm[match(available_genes, gene), -1L, with = FALSE]) / 100)
rownames(gene_expression) <- available_genes
colnames(gene_expression) <- geo_meta$sample_id
absolute_score <- t(gene_expression) %*% external_weights
absolute_score[, coverage_weight == 0] <- NA_real_
sample_scores <- bind_rows(lapply(program_levels, function(p) geo_meta %>% mutate(
  program = p, absolute_score = absolute_score[, p])))
## Set ID identifies the independent experiment day. All technical wells of a
## condition (including the four DMSO wells in each time-course set) are averaged.
bio <- sample_scores %>% group_by(program, experiment, biological_set, arm, vem_uM, cobi_uM, egf_h) %>%
  summarise(score = mean(absolute_score), n_technical_wells = n(), .groups = "drop")
baseline <- bio %>% filter(arm == "DMSO") %>% select(program, experiment, biological_set, dmso_score = score)
stopifnot(!anyDuplicated(baseline[c("program", "experiment", "biological_set")]))
bio <- bio %>% left_join(baseline, by = c("program", "experiment", "biological_set")) %>%
  mutate(dmso_displacement = score - dmso_score)
timecourse <- bio %>% filter(experiment == "EGF time course", arm != "DMSO")
time_zero <- timecourse %>% filter(egf_h == 0) %>% select(program, biological_set, arm, time_zero_score = score)
timecourse <- timecourse %>% left_join(time_zero, by = c("program", "biological_set", "arm")) %>%
  mutate(egf_displacement = score - time_zero_score)
stopifnot(nrow(timecourse) == length(program_levels) * 2L * 2L * 7L,
  all(timecourse$n_technical_wells == 2L), !anyNA(timecourse$time_zero_score))
fwrite(geo_meta, file.path(out, "gse127988_sample_metadata.tsv"), sep = "\t")
fwrite(sample_scores, file.path(out, "gse127988_sample_scores.tsv"), sep = "\t", na = "NA")
fwrite(bio, file.path(out, "gse127988_biological_scores.tsv"), sep = "\t", na = "NA")
fwrite(timecourse, file.path(out, "gse127988_timecourse_scores.tsv"), sep = "\t", na = "NA")

## Time-integrated contrasts, one value per biological day. No genes, wells or
## time points are counted as independent biological replicates.
auc <- function(time, value) {
  o <- order(time); time <- time[o]; value <- value[o]
  stopifnot(identical(as.numeric(time), c(0, .5, 1, 2, 3, 4, 8)))
  sum(diff(time) * (head(value, -1) + tail(value, -1)) / 2) / 8
}
area <- timecourse %>% group_by(program, biological_set, arm) %>% summarise(
  egf_auc_mean = auc(egf_h, egf_displacement), total_auc_mean = auc(egf_h, dmso_displacement), .groups = "drop")
v <- area %>% filter(arm == "Vem")
c <- area %>% filter(arm == "Vem + Cobi")
paired <- inner_join(v, c, by = c("program", "biological_set"), suffix = c("_vem", "_combo"))
contrasts <- bind_rows(
  transmute(v, program, biological_set, contrast = "EGF during vemurafenib", effect = egf_auc_mean),
  transmute(c, program, biological_set, contrast = "EGF during Vem + Cobi", effect = egf_auc_mean),
  transmute(paired, program, biological_set, contrast = "Cobimetinib x EGF interaction",
    effect = egf_auc_mean_combo - egf_auc_mean_vem),
  transmute(paired, program, biological_set, contrast = "Cobi change in total program level",
    effect = total_auc_mean_combo - total_auc_mean_vem))
test_mean <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2L || sd(x) == 0) return(tibble(n_biological = n, mean_effect = mean(x),
    ci_low = NA_real_, ci_high = NA_real_, nominal_p = NA_real_, sign_test_p = NA_real_))
  a <- t.test(x, mu = 0)
  tibble(n_biological = n, mean_effect = mean(x), ci_low = a$conf.int[1], ci_high = a$conf.int[2],
    nominal_p = a$p.value, sign_test_p = if (any(x != 0)) binom.test(sum(x > 0), sum(x != 0), p = .5)$p.value else 1)
}
tests <- contrasts %>% group_by(program, contrast) %>% group_modify(~test_mean(.x$effect)) %>% ungroup() %>%
  mutate(bh_q = p.adjust(nominal_p, method = "BH")) %>% left_join(predictions, by = "program") %>%
  left_join(program_coverage, by = "program") %>%
  mutate(expected_direction = case_when(contrast == "EGF during vemurafenib" ~ expected_egf_direction,
    contrast == "Cobi change in total program level" ~ expected_cobimetinib_direction, TRUE ~ "exploratory"),
    direction_agrees = case_when(expected_direction == "increase" ~ mean_effect > 0,
      expected_direction == "decrease" ~ mean_effect < 0, TRUE ~ NA),
    interpretation = "Exploratory biological-day t test (n=2); exact sign test also reported; no pseudoreplication")
fwrite(contrasts, file.path(out, "gse127988_day_level_contrasts.tsv"), sep = "\t")
fwrite(tests, file.path(out, "gse127988_validation_tests.tsv"), sep = "\t", na = "NA")

arm_colors <- c(Vem = "#247A96", `Vem + Cobi` = "#B95B32")
for (measure in c("egf_displacement", "dmso_displacement")) {
  d <- timecourse %>% mutate(displacement = .data[[measure]],
    panel = paste0(program, " | ", coherence$candidate_role[match(program, coherence$program)],
      " | ", program_coverage$n_measured[match(program, program_coverage$program)], "/",
      program_coverage$n_frozen[match(program, program_coverage$program)], " genes matched"))
  means <- d %>% group_by(panel, arm, egf_h) %>% summarise(displacement = mean(displacement), .groups = "drop")
  p <- ggplot(d, aes(egf_h, displacement, color = arm)) + geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = .4) +
    geom_line(aes(group = interaction(arm, biological_set)), linewidth = .5, alpha = .4) +
    geom_point(aes(shape = factor(biological_set)), size = 2.4) +
    geom_line(data = means, aes(group = arm), linewidth = .9) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") + scale_color_manual(values = arm_colors) +
    scale_x_continuous(breaks = c(0, .5, 1, 2, 3, 4, 8)) +
    labs(x = "Hours after EGF addition", y = "Frozen program displacement", color = NULL, shape = "Biological set",
      title = if (measure == "egf_displacement") "EGF response under RAF / MEK inhibition" else "Program state relative to vehicle",
      subtitle = if (measure == "egf_displacement") "Relative to the same drug arm at time zero" else "Relative to same-day DMSO",
      caption = "Points: two biological days, each averaging duplicate wells. Thick lines: day means.\nPositive program displacement means increased gene expression; no sign reversals.\nVemurafenib 1 uM; cobimetinib 1 uM; EGF 100 ng/mL.") +
    theme_classic(base_size = 12) + theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
      legend.position = "bottom", plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9))
  ggsave(file.path(out, paste0("gse127988_", measure, ".pdf")), p, width = 8.8, height = 11.3, device = cairo_pdf)
}

## Validation heatmaps keep discovery gene order; chronological columns are
## essential for the perturbation test. No validation-driven sample clustering.
conditions <- geo_meta %>% filter(experiment == "EGF time course", arm != "DMSO") %>%
  distinct(arm, egf_h) %>% arrange(factor(arm, levels = c("Vem", "Vem + Cobi")), egf_h) %>%
  mutate(condition_id = paste0(arm, " | ", egf_h, " h"))
gene_displacement <- matrix(NA_real_, length(gene_order), nrow(conditions), dimnames = list(gene_order, conditions$condition_id))
for (j in seq_len(nrow(conditions))) {
  day_effects <- vapply(3:4, function(day) {
    target <- geo_meta$sample_id[geo_meta$biological_set == day & geo_meta$arm == conditions$arm[j] & geo_meta$egf_h == conditions$egf_h[j]]
    reference <- geo_meta$sample_id[geo_meta$biological_set == day & geo_meta$arm == conditions$arm[j] & geo_meta$egf_h == 0]
    rowMeans(gene_expression[, target, drop = FALSE]) - rowMeans(gene_expression[, reference, drop = FALSE])
  }, numeric(length(available_genes)))
  gene_displacement[available_genes, j] <- rowMeans(day_effects)
}
write_matrix(gene_displacement, "gse127988_gene_displacements.tsv.gz")
ca <- data.frame(Arm = conditions$arm, row.names = conditions$condition_id)
for (labelled in c(FALSE, TRUE)) save_heatmap(gene_displacement,
  paste0("gse127988_gene_displacement", if (labelled) "_labelled" else "", ".pdf"),
  "GSE127988 | EGF displacement | Discovery gene order", c(-1, 0, 1), "Log2 change", ca,
  labelled = labelled, column_labels = paste0(conditions$egf_h, " h"), row_annotation = annotation, width = 10.5, height = 12.5)
program_displacement <- t(gene_displacement[available_genes, , drop = FALSE]) %*% external_weights
program_displacement <- t(program_displacement)
pa <- data.frame(Program = factor(program_levels, levels = program_levels), row.names = program_levels)
save_heatmap(program_displacement, "gse127988_program_heatmap.pdf", "Core programs | EGF response in A375",
  c(-max(abs(program_displacement), na.rm = TRUE), 0, max(abs(program_displacement), na.rm = TRUE)),
  "Program change", ca, row_annotation = pa, labelled = TRUE,
  column_labels = paste0(conditions$egf_h, " h"), width = 9.5, height = 4.2,
  row_labels = paste0(program_levels, " (", program_coverage$n_measured, "/", program_coverage$n_frozen, " genes)"))
stopifnot(identical(unname(tools::md5sum(freeze_paths)), unname(freeze_hashes)))

## Render the COMPLETE experimental design, including the dose experiment.
design_counts <- geo_meta %>% group_by(experiment,vem_uM,cobi_uM,egf_h) %>%
  summarise(n_wells=n(),n_biological_days=n_distinct(biological_set),.groups="drop")
fwrite(design_counts,file.path(out,"gse127988_complete_design.tsv"),sep="\t")
dose_levels <- sort(unique(geo_meta$vem_uM[geo_meta$experiment=="Vemurafenib dose"]))
dose_design <- matrix(4,1,length(dose_levels),dimnames=list("Vemurafenib dose (uM)",as.character(dose_levels)))
time_design <- matrix(NA_real_,3,7,dimnames=list(c("DMSO","Vem 1 uM","Vem 1 uM + Cobi 1 uM"),c("0","0.5","1","2","3","4","8")))
time_design[1,1] <- 8;time_design[2:3,] <- 4
design_heatmap <- function(m,title) draw_heatmap(m,scale_rows=FALSE,cluster_rows=FALSE,cluster_columns=FALSE,
  show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=12),column_names_gp=grid::gpar(fontsize=12),
  hm_name="Wells",value_range=c(0,4,8),value_colors=c("#F7FBFF","#6BAED6","#08519C"),
  na_col="#DDDDDD",title=title,column_title_gp=grid::gpar(fontsize=14,fontface="bold"),
  cell_fun=function(j,i,x,y,w,h,fill) grid::grid.text(if(is.na(m[i,j])) "-" else m[i,j],x,y,
    gp=grid::gpar(fontsize=13,col=if(!is.na(m[i,j])&&m[i,j]>=6) "white" else "black")))
cairo_pdf(file.path(out,"gse127988_complete_design.pdf"),width=11,height=6.6)
grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1,heights=grid::unit(c(1,1.7),"null"))))
for(i in 1:2) {
  grid::pushViewport(grid::viewport(layout.pos.row=i))
  ht <- if(i==1) design_heatmap(dose_design,"A | 24-hour vemurafenib dose series: 32 samples") else
    design_heatmap(time_design,"B | Cobimetinib x EGF time (hours): 64 samples")
  ComplexHeatmap::draw(ht,newpage=FALSE,padding=grid::unit(c(7,5,5,5),"mm"));grid::popViewport()
}
grid::popViewport();grid::grid.text("Numbers = sequenced wells. Two biological days per experiment. Gray cells were not measured.",
  x=.5,y=grid::unit(4,"mm"),gp=grid::gpar(fontsize=11));dev.off()

dose_scores <- bio %>% filter(experiment=="Vemurafenib dose") %>%
  mutate(dose_label=factor(as.character(vem_uM),levels=as.character(dose_levels)))
dose_means <- dose_scores %>% group_by(program,dose_label) %>% summarise(dmso_displacement=mean(dmso_displacement),.groups="drop")
dp <- ggplot(dose_scores,aes(dose_label,dmso_displacement))+
  geom_hline(yintercept=0,color="#AAAAAA",linewidth=.4)+
  geom_line(aes(group=biological_set),color="#247A96",alpha=.45)+
  geom_point(aes(shape=factor(biological_set)),color="#247A96",size=2.5)+
  geom_line(data=dose_means,aes(group=1),color="#247A96",linewidth=1)+
  facet_wrap(~program,ncol=1,scales="free_y")+
  labs(x="Vemurafenib (uM); 0 = vehicle",y="Core expression change from same-day vehicle",shape="Biological set",
    title="The separate vemurafenib dose experiment",
    subtitle=paste0("All eight doses shown; fixed ",length(gene_order),"-gene programs"),
    caption="Points: technical-well means within each biological day. Line: mean across days.\nPositive scores mean increased raw expression. Dose levels are spaced as categories.")+
  theme_classic(base_size=13)+theme(legend.position="bottom",strip.background=element_blank(),strip.text=element_text(face="bold"),
    plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
ggsave(file.path(out,"gse127988_vemurafenib_dose_series.pdf"),dp,width=8.3,height=10.3,device=cairo_pdf)

interaction_curve <- timecourse %>% filter(arm=="Vem") %>% select(program,biological_set,egf_h,vem=egf_displacement) %>%
  inner_join(timecourse %>% filter(arm=="Vem + Cobi") %>% select(program,biological_set,egf_h,combo=egf_displacement),
    by=c("program","biological_set","egf_h")) %>% mutate(interaction=combo-vem)
interaction_mean <- interaction_curve %>% group_by(program,egf_h) %>% summarise(interaction=mean(interaction),.groups="drop")
ip <- ggplot(interaction_curve,aes(egf_h,interaction))+geom_hline(yintercept=0,color="#888888",linewidth=.4)+
  geom_line(aes(group=biological_set),color="#7757A4",alpha=.4)+geom_point(aes(shape=factor(biological_set)),color="#7757A4",size=2.5)+
  geom_line(data=interaction_mean,color="#7757A4",linewidth=1)+facet_wrap(~program,ncol=1,scales="free_y")+
  scale_x_continuous(breaks=c(0,.5,1,2,3,4,8))+
  labs(x="Hours after EGF",y="Cobimetinib's change in the EGF-induced response",shape="Biological set",
    title="Factorial interaction | Cobimetinib x EGF",
    subtitle="[Vem+Cobi+EGF - Vem+Cobi at 0 h] - [Vem+EGF - Vem at 0 h]",
    caption="Negative = cobimetinib reduces EGF induction. Independent replication remains two biological days.\nThe dose experiment is separate; there is no EGF-only or cobimetinib-only series.")+
  theme_classic(base_size=13)+theme(legend.position="bottom",strip.background=element_blank(),strip.text=element_text(face="bold"),
    plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
ggsave(file.path(out,"gse127988_factorial_interaction.pdf"),ip,width=9,height=10.3,device=cairo_pdf)
fwrite(interaction_curve,file.path(out,"gse127988_factorial_interaction.tsv"),sep="\t")

## Show the diagnosed selection failure and a retained drug-associated gene.
examples <- tibble(compound=c("vemurafenib","vemurafenib","cobimetinib","vemurafenib"),
  dose=c(.5,.5,.05,.5),gene=c("OR52D1","KRT76","FAM83G","AKAP12"))
example_points <- example_labels <- list()
for(i in seq_len(nrow(examples))) {
  e <- examples[i,];idx <- which(meta$split=="Training"&meta$compound==e$compound&meta$dose==e$dose&is.finite(meta$prism_sensitization))
  xx <- effects[idx,e$gene]; yy <- -meta$prism_sensitization[idx];dom <- which.max((xx-mean(xx))^2)
  all <- pearson_test(xx,yy);drop <- pearson_test(xx[-dom],yy[-dom])
  panel <- paste0(e$gene," | ",e$compound," ",e$dose," uM")
  example_points[[i]] <- tibble(panel=panel,x=xx,y=yy,highlight=seq_along(idx)==dom,line=meta$tahoe_cell_name[idx])
  example_labels[[i]] <- tibble(panel=panel,label=paste0("All: r=",sprintf("%.2f",all$pearson_r),"\nWithout ",meta$tahoe_cell_name[idx[dom]],": r=",
    if(is.finite(drop$pearson_r)) sprintf("%.2f",drop$pearson_r) else "undefined"))
}
ep <- bind_rows(example_points);el <- bind_rows(example_labels)
xp <- ggplot(ep,aes(x,y))+geom_point(aes(color=highlight),size=2.1)+
  geom_smooth(method="lm",formula=y~x,se=FALSE,color="#247A96",linewidth=.8)+
  geom_text(data=el,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=3.3)+
  facet_wrap(~panel,ncol=2,scales="free")+scale_color_manual(values=c(`FALSE`="#506070",`TRUE`="#D27631"),
    labels=c(`FALSE`="Other training lines",`TRUE`="Largest variance contribution"),name=NULL)+
  labs(x="Gene-expression displacement",y="PRISM resistance (-sensitization)",
    title="Why correlation rank alone failed",
    caption="The earlier selector admitted tiny changes present in one line. Removing that line can make r undefined.\nRevised drug candidates must remain measurable and significant after every single-line deletion.")+
  theme_classic(base_size=13)+theme(aspect.ratio=1,legend.position="bottom",plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=10))
ggsave(file.path(out,"gene_selection_failure_examples.pdf"),xp,width=10.8,height=9.8,device=cairo_pdf)


## Complete gene-by-sample view, including both experiments and every well.
all_sample_delta <- gene_expression
for(i in seq_len(nrow(geo_meta))) {
  mm <- geo_meta[i,];refs <- geo_meta$sample_id[geo_meta$experiment==mm$experiment &
    geo_meta$biological_set==mm$biological_set & geo_meta$arm=="DMSO"]
  all_sample_delta[,i] <- gene_expression[,i]-rowMeans(gene_expression[,refs,drop=FALSE])
}
all_order <- order(factor(geo_meta$experiment,levels=c("Vemurafenib dose","EGF time course")),
  factor(geo_meta$arm,levels=c("DMSO","Vem","Vem + Cobi")),geo_meta$vem_uM,geo_meta$egf_h,geo_meta$biological_set,geo_meta$sample_number)
all_ca <- data.frame(Experiment=geo_meta$experiment[all_order],Arm=geo_meta$arm[all_order],
  `Vem uM`=factor(geo_meta$vem_uM[all_order]),`EGF h`=factor(geo_meta$egf_h[all_order]),
  Day=factor(geo_meta$biological_set[all_order]),row.names=geo_meta$sample_id[all_order],check.names=FALSE)
write_matrix(all_sample_delta[gene_order,all_order],"gse127988_all_sample_gene_displacements.tsv.gz")
save_heatmap(all_sample_delta[gene_order,all_order],"gse127988_all_sample_genes.pdf",
  "GSE127988 | All 96 wells | Change from same-day vehicle",c(-1,0,1),"Log2 change",all_ca,
  annotation,column_labels=colnames(all_sample_delta)[all_order],width=15,height=12.8)

fwrite(tibble(input=names(paths),path=unname(paths)),file.path(out,"input_manifest.tsv"),sep="\t")
fwrite(tibble(check=c("30-60 frozen genes","Four whole modules; 8/22/13/9 genes","Fixed positive weights sum to one",
  "All 52 genes mapped in GSE127988","Frozen files unchanged after external projection","All three correlation panels retain identical row/column order"),
  passed=c(length(gene_order)>=30&length(gene_order)<=60,identical(as.integer(coherence$n_genes),c(8L,22L,13L,9L)),
    all(abs(colSums(weights)-1)<1e-12),length(available_genes)==length(gene_order),
    identical(unname(tools::md5sum(freeze_paths)),unname(freeze_hashes)),TRUE)),file.path(out,"verification.tsv"),sep="\t")
writeLines(c("# RAF–EGF pilot: 52 genes and an explicit mechanistic Y", "",
  "Run from /work: `Rscript /work/tahoe_prism_v2/scripts/pilot_tahoe_cytosig_programs.R`.","",
  "The final model contains four whole programs: P1 8 genes, P2 22 genes, P3 13 genes and P4 9 genes (52 total). The former 12-gene result is superseded. A 9-gene program is part of this larger model; it is not the entire analysis.","",
  "Start with [the rendered Y](raf_egf_mechanistic_y.pdf), [matched correlation heatmaps](core_correlation_comparison.pdf), [all 52 genes in the hierarchy](gene_hierarchy.pdf), and [the external gene heatmap](gse127988_gene_displacement_labelled.pdf). The Y is also available as PNG and editable SVG. Figures use the existing /work/src/R/lib_heatmap.R wrapper; gene order is frozen across sources. The PDF plots have been rendered and checked.","",
  "## What the Y means", "",
  "The solid signaling arrows show the established EGF–EGFR–RAS/RAF-dimer route and oncogenic BRAF route converging on MEK–ERK. Vemurafenib acts on the BRAF arm and cobimetinib acts at MEK. Gerosa et al. independently establish receptor-driven ERK activity and persistence under RAF/MEK inhibition. These mechanistic arrows are not estimated from RNA correlations.",
  "Dotted links place expression modules as hypotheses: P1 (AREG/ERRFI1) is a receptor-feedback candidate; P2 (FOSL1/AKAP12) and P3 (EGR1/JUNB) are shared-response candidates. P4 (IFIT2/OASL) is enriched for cobimetinib-selected genes and is left without a pathway position. There is no evidence that P1 as a whole is upstream of P2/P3, nor that either downstream module causes the other. The data do not establish a complete three-program causal Y or a BRAF-specific transcriptional arm.",
  "Source coherence and source activation are different. P4 genes co-vary in both datasets, but their mean EGF response is near zero. P2 has weaker Tahoe coherence (r=0.37), below the stated 0.40 source-support label. P3 has no directly selected PRISM genes; it remains a molecular response candidate, not a proven survival program. The Y labels and frozen_branch_assignments.tsv make these limitations explicit.","",
  "## Selection and freezing", "",
  "All 19,339 genes measured in Tahoe and CytoSig are eligible before selection. All doses (0.05, 0.5, 5 uM) and the original 26-training/7-holdout split are retained. There is no old 3,000-gene or drug-activity gate. Both drugs have adequate matched PRISM overlap (at least 21 training lines per dose).", 
  "The prior strongest-|r| screen admitted one-line artifacts: for example OR52D1/KRT76 had tiny changes concentrated in LOX-IMVI, and deletion of that line made correlation undefined. Drug candidates now require nominal Pearson p<0.05, abundance >=10 CPM in at least half the matched training lines (and at least eight), and finite same-sign p<0.05 after every single-line deletion. Rank by the weakest deletion |r|, breaking ties deterministically. Retain up to 100 genes per drug: 74 vemurafenib and 100 cobimetinib candidates. This is a gene detectability/influence rule; it does not discard drug conditions.",
  "EGF selection uses one median per independent study (51 profiles, 23 studies), a nominal two-sided sign test p<0.05, and measurable Tahoe abundance in at least half the training observations. Retain the top 100 absolute study-median EGF effects. The combined pool contains 272 unique genes. RAB43 and SCHIP1 have no finite EGF-study observations and cannot supply correlations; they are explicitly audited rather than imputed. The clustering pool is 270 genes.",
  "Average the Tahoe-training and CytoSig-study Pearson correlation matrices with equal weight. Use average-linkage hclust on 1-r and the existing height=0.6 coherence cut. To implement the user's revised 30–60-total preference, retain whole groups with >=8 genes: 52 genes in four modules. This minimum size was chosen during discovery revision, not prespecified in the original pilot; it is recorded in selection_parameters.tsv. It was not chosen from GSE outcomes. The other 218 candidates remain in the audit and do not enter scoring. The original tree is pruned without reclustering or changing retained merge heights.",
  "The rank-versus-silhouette curve covers 2–8 clusters in the frozen set and is diagnostic; it does not select the result. Each final program has equal positive gene weights summing to one. A positive score means increased expression. Program-to-PRISM scatterplots show every drug/dose, with separate training and holdout trends. Selection enrichment uses the 270-gene pool as background and reports one-sided Fisher tests with BH correction; these conditional descriptive tests do not establish causal hierarchy.",
  "All final gene names, HGNC annotations, selection evidence, signs and weights are in gene_inclusion_evidence.tsv. Complete candidate decisions are in candidate_gene_inclusion_evidence.tsv. Program membership, source annotation, predictions and branch hypotheses are frozen and hashed before this run reads the external matrix. GSE127988 and the Tahoe holdout had already been inspected in prior iterations; this revision is a retrospective check with discovery-only gene selection, not fresh independent validation.","",
  "## Gene mapping and complete external design", "",
  "All 52 frozen genes map by exact symbol to the published GSE127988 matrix. No replacement genes or score renormalization are needed. The old 88 apparent misses comprised two resolvable historical symbols (CCN2/CTGF, FBXL21P/FBXL21) and 86 genes absent from the published processed matrix. Its processing description excludes all-zero genes; absence cannot identify the reason for each particular gene. The original audit remains in ../23_raf_egf_repair/old_gene_identifier_audit.tsv.",
  "GSE127988 contains two experiments, both rendered in gse127988_complete_design.pdf. The 32-sample dose experiment has vehicle plus seven vemurafenib doses. The 64-sample time experiment has two backgrounds (vemurafenib 1 uM; vemurafenib 1 uM + cobimetinib 1 uM), seven EGF times (0, 0.5, 1, 2, 3, 4, 8 h), and vehicle controls. It is a drug-background-by-EGF-time design, not a complete vemurafenib x cobimetinib x EGF presence factorial. No EGF-only or cobimetinib-only series exists, and later time points lack time-matched EGF-free controls.",
  "Score individual RPM columns with log2(1+RPM/100) and the frozen equal weights; average duplicate wells within each biological day. Analyze the two experiments separately. S15/S16 have inconsistent GEO versus matrix set labels; use the published matrix's set IDs and preserve both in the discrepancy audit. The EGF experiment's set labels agree.",
  "EGF induction uses the same drug arm at time zero. The cobimetinib x EGF contrast is [combo(t)-combo(0)] - [vem(t)-vem(0)]. Total program state also uses same-day DMSO. Integrate 0–8 h by trapezoids, yielding one value per biological day. Report n=2 t tests/intervals, exact sign tests and BH correction across the 16 external contrasts. Time points, genes and wells are never counted as independent biological replicates.","",
  "## Observed result", "",
  "P1/P2/P3 increase after EGF under vemurafenib, and adding cobimetinib reduces their EGF-induced component. P4 shows little EGF displacement. P3's total level is not reduced by cobimetinib even though its EGF-induced component is reduced; those are different contrasts. None of the 16 external contrasts passes BH q<0.05 (smallest q about 0.094). Thus the molecular directions are suggestive, not conclusive validation.",
  "The whole-program PRISM association is weak for most training drug/dose combinations. P4 is associated with cobimetinib response at 0.5 and 5 uM in training; no program-drug-dose holdout association passes the within-holdout BH correction. Vemurafenib 0.5 uM yields nominal holdout correlations for P2/P3, but these are dose-specific and do not establish scalable function. PRISM supplies gene-selection evidence, Gerosa supplies mechanistic/functional evidence, and this RNA-seq series supplies molecular displacement. The full pilot claim is not established.","",
  paste0(coherence$program,": ",coherence$genes,". Mean within-program r: Tahoe ",sprintf("%.2f",coherence$tahoe_mean_r),"; EGF ",sprintf("%.2f",coherence$cytosig_egf_mean_r),"."),"",
  "Sources: [Gerosa et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8009031/), [GEO GSE127988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988), [HGNC](https://hgnc.genenames.org/download/)."),file.path(out,"README.md"))
figures <- sort(list.files(out,pattern="\\.pdf$"))
fwrite(tibble(file=figures,title=sub("\\.pdf$","",figures)),file.path(out,"FIGURE_INDEX.tsv"),sep="\t")
print(coherence);print(program_coverage)
message("Completed 52-gene discovery, rendered mechanistic Y, fixed-order heatmaps and all external experiments.")
