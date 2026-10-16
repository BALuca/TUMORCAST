#!/usr/bin/env Rscript
# Prepare all source-PDEX cell lines with functional coverage for stages 11/12.
# The shared gene universe and aggregate-plates-first normalization match 09.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(cluster)
  library(RcppML)
  library(glmnet)
  library(arrow)
})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
output_dir <- file.path(work_dir, 'results/program_pilot/11_panel_functional_support_programs_nmf')
input_dir <- file.path(output_dir, 'data/molecular_inputs')
audit_dir <- file.path(output_dir, 'data/panel_preparation')
for (path in c(input_dir, audit_dir)) dir.create(path, recursive=TRUE, showWarnings=FALSE)
data.table::setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
minimum_conditions <- 20L
minimum_drugs <- 10L

# PDEX contains 50 lines. Public Tahoe identity metadata describes a larger set;
# identities without profiles in this pinned PDEX source are explicitly excluded.
identity <- unique(fread(file.path(work_dir,'almanac_conditionals/diagnosis/public_tahoe_cell_metadata.tsv'))[
  , .(depmap_id=Cell_ID_DepMap, cell_line=Cell_ID_Cellosaur, depmap_name=cell_name)])
rowgroups <- fread(file.path(work_dir,'almanac_conditionals/diagnosis/pdex_rowgroups.tsv'))
almanac_map <- fread(file.path(work_dir,'almanac_conditionals/repaired/results/cell_mapping.tsv'))
prism_ids <- fread(file.path(cache_dir,'prism/secondary-screen-replicate-collapsed-logfold-change.csv'),select=1)[[1]]
identity[, `:=`(molecular_source=cell_line %chin% c(rowgroups$cell_min,rowgroups$cell_max),
  prism_present=nzchar(depmap_id) & depmap_id %chin% prism_ids,
  almanac_present=nzchar(depmap_id) & depmap_id %chin% almanac_map[nzchar(depmap_id),depmap_id])]
identity[, CELLNAME := almanac_map$CELLNAME[match(depmap_id, almanac_map$depmap_id)]]
fwrite(identity,file.path(audit_dir,'cell_identity_coverage.tsv'),sep='\t')
cell_map <- identity[molecular_source & (prism_present | almanac_present)]
check(!anyDuplicated(cell_map$depmap_id), 'Ambiguous public cell identity')
# Preserve the historical seven-cell names and their NMF seed positions.
original_map <- fread(local_cell_mapping_path)
idx <- match(cell_map$depmap_id, original_map$depmap_id)
cell_map[!is.na(idx), depmap_name := original_map$depmap_name[idx[!is.na(idx)]]]
cell_map[, prior_order := match(depmap_id,original_map$depmap_id)]
setorder(cell_map, prior_order, depmap_id, na.last=TRUE)
gene_table <- fread(file.path(tahoe_dir,'genes.tsv'))[order(gene_index)]
audit <- fread(gene_audit_path)
genes <- sort(audit[protein_coding & median_TPM>=5 & fraction_ge_1_TPM>=0.8 & gene %chin% gene_table$feature,gene])
gene_index <- match(genes,gene_table$feature)
conditions <- fread(file.path(tahoe_dir,'perturbations.tsv'))[order(condition_index)]
raw_metadata <- fread(file.path(cache_dir,'response_effect_metadata.tsv'))[available==1]
old_metadata <- fread(tahoe_metadata_path)
plates <- fread(file.path(tahoe_dir,'pdex_plate_audit.tsv.gz'))
plate_keys <- plates[,.(contributing_plates=paste(sort(unique(plate)),collapse=',')),by=.(depmap_id,perturbation_id)]
references <- read_npy_matrix(file.path(tahoe_dir,'pdex_plate_reference.npy'))
ref_md <- fread(file.path(tahoe_dir,'pdex_plate_reference_metadata.tsv'))[order(array_row)]
cached_cells <- fread(file.path(tahoe_dir,'cells.tsv'))
legacy <- qread(file.path(local_output_dir,'data/controlled_molecular_inputs.qs'))

# Read-only cached Parquet ranges supply the four additional ALMANAC cell lines.
# Missing row groups are fetched into a separate sparse patch file, never into
# the existing source cache. Neighboring validated range boundaries delimit
# complete missing groups, avoiding a download of the 40-GB source file.
chunk_files <- list.files(file.path(work_dir,'almanac_conditionals/diagnosis/pdex_projection_chunks'),full.names=TRUE)
chunk_index <- rbindlist(lapply(chunk_files,function(path) {
  z <- jsonlite::fromJSON(path)
  data.table(group=z$rowgroup,start=min(z$ranges$start),end=max(z$ranges$end),
    cells=list(unique(z$profiles$cell_line)))
}))
setorder(chunk_index,group)
source_path <- file.path(work_dir,'almanac_conditionals/diagnosis/pdex_selected_sparse.parquet')
source_reader <- ParquetFileReader$create(source_path)
patch_path <- file.path(input_dir,'pdex_missing_rowgroups.parquet')
patch_manifest_path <- file.path(input_dir,'pdex_fetched_ranges.tsv')
source_url <- 'https://huggingface.co/datasets/tahoebio/tahoe-de-rhaister/resolve/c7963cf334bec0683225d41c9586d900ca6303a2/pdex/all_plates_pdex.parquet'
source_size <- file.info(source_path)$size
if (!file.exists(patch_path)) {
  footer_path <- file.path(work_dir,'almanac_conditionals/diagnosis/pdex_metadata.parquet')
  footer <- readBin(footer_path,'raw',file.info(footer_path)$size)
  con <- file(patch_path,'w+b')
  writeBin(charToRaw('PAR1'),con)
  seek(con,source_size-length(footer)+4,origin='start',rw='write')
  writeBin(footer[-(1:4)],con)
  close(con)
}
fetch_group <- function(group_id) {
  # Explicit subset indexing avoids data.table argument/column shadowing.
  before <- chunk_index[which(chunk_index$group < group_id)][.N]
  after <- chunk_index[which(chunk_index$group > group_id)][1]
  check(nrow(before)==1 && nrow(after)==1, 'Missing boundary group cannot be delimited safely')
  start <- before$end
  end <- after$start
  fetched <- if(file.exists(patch_manifest_path)) fread(patch_manifest_path) else data.table(start=numeric(),end=numeric())
  if (!any(fetched$start<=start & fetched$end>=end)) {
    message('Fetching missing PDEX interval: ',start,'--',end)
    h <- curl::new_handle(range=sprintf('%.0f-%.0f',start,end-1), timeout=120, followlocation=TRUE)
    response <- curl::curl_fetch_memory(paste0(source_url,'?panel=',format(Sys.time(),'%Y%m%d%H%M%S')),h)
    check(response$status_code==206 && length(response$content)==end-start,'Invalid ranged PDEX response')
    con <- file(patch_path,'r+b')
    seek(con,start,origin='start',rw='write')
    writeBin(response$content,con)
    close(con)
    fetched <- rbind(fetched,data.table(start=start,end=end),fill=TRUE)
    fwrite(fetched,patch_manifest_path,sep='\t')
  }
  ParquetFileReader$create(patch_path)$ReadRowGroup(as.integer(group_id),c(0L,1L,2L,3L,11L,12L))
}

metadata_parts <- list()
for (id in cell_map$depmap_id) {
  path <- file.path(input_dir,paste0(id,'.qs'))
  recipient <- cell_map$depmap_name[match(id,cell_map$depmap_id)]
  if(file.exists(path)) {
    prepared <- qread(path)
    check(identical(colnames(prepared$response),genes) && identical(prepared$normalizer,
      'paired_trimmed_median_log_ratio_cp10k_v1'),paste('Incompatible molecular cache',id))
    metadata_parts[[id]] <- prepared$metadata
    message('Reusing molecular inputs: ',recipient)
    next
  }
  message('Preparing molecular inputs: ',recipient)
  if (id %chin% cached_cells$depmap_id) {
    md <- copy(raw_metadata[depmap_id==id])
    md[, observation_id := paste(depmap_id,perturbation_id,sep='__')]
    old <- old_metadata[match(md$observation_id,observation_id)]
    md[, contributing_plates := plate_keys$contributing_plates[match(paste(id,perturbation_id),
      paste(plate_keys$depmap_id,plate_keys$perturbation_id))]]
    reference_rows <- which(ref_md$depmap_id==id)
    ref_counts <- references[reference_rows,,drop=FALSE]
    plate_ids <- ref_md$plate[reference_rows]
    if (id %chin% names(legacy$null)) {
      response <- legacy$response[match(md$observation_id,rownames(legacy$response)),,drop=FALSE]
      check(identical(rownames(response),md$observation_id),'Legacy row mismatch')
      md[, `:=`(treated_total=old$treated_total, reference_total=old$reference_total,
        treated_correction_factor=old$treated_correction_factor)]
    } else {
      shard <- read_npy_matrix(file.path(cache_dir,unique(md$shard_path)))
      response <- matrix(NA_real_,nrow(md),length(genes),dimnames=list(md$observation_id,genes))
      md[, `:=`(treated_total=NA_real_,reference_total=NA_real_,treated_correction_factor=NA_real_)]
      for(i in seq_len(nrow(md))) {
        selected_plates <- as.integer(strsplit(md$contributing_plates[i],',',fixed=TRUE)[[1]])
        ref <- colMeans(ref_counts[match(selected_plates,plate_ids),,drop=FALSE])
        treated <- shard[md$shard_row[i]+1L,]+ref
        tolerance <- max(1e-6,32*2^-23*max(abs(treated),abs(ref),1))
        check(min(treated)>=-tolerance,'Material negative treated count')
        treated <- pmax(treated,0)
        ts <- if(!is.na(old$treated_total[i])) old$treated_total[i] else sum(treated)
        rs <- if(!is.na(old$reference_total[i])) old$reference_total[i] else sum(ref)
        factor <- if(!is.na(old$treated_correction_factor[i])) old$treated_correction_factor[i] else composition_factor(treated,ref)
        response[i,] <- log2(1+treated[gene_index]/ts*10000/factor)-log2(1+ref[gene_index]/rs*10000)
        md[i,`:=`(treated_total=ts,reference_total=rs,treated_correction_factor=factor)]
      }
      rm(shard)
    }
  } else {
    cvcl <- cell_map$cell_line[match(id,cell_map$depmap_id)]
    # Pure row groups plus all cached mixed groups known to contain this cell.
    pure <- rowgroups[cell_min==cvcl & cell_max==cvcl,group]
    mixed <- chunk_index[vapply(cells,function(z) cvcl %in% z,logical(1)),group]
    groups <- sort(unique(c(pure,mixed)))
    treated_sum <- matrix(0,nrow(conditions),nrow(gene_table))
    profile_plates <- vector('list',nrow(conditions))
    reference_list <- list()
    for(group in groups) {
      tab <- as.data.table(if(group %in% chunk_index$group)
        source_reader$ReadRowGroup(as.integer(group),c(0L,1L,2L,3L,11L,12L)) else fetch_group(group))
      tab <- tab[cell_line==cvcl]
      keys <- unique(tab[,.(target,plate)])
      for(j in seq_len(nrow(keys))) {
        g <- tab[target==keys$target[j] & plate==keys$plate[j]]
        check(nrow(g)==nrow(gene_table) && identical(g$feature,gene_table$feature),'PDEX profile gene order mismatch')
        # Single-compound target syntax is fixed by the source, not evaluated.
        target <- keys$target[j]
        parts <- regmatches(target,regexec("^\\[\\(['\"](.*)['\"], ([0-9.eE+-]+), ['\"]uM['\"]\\)\\]$",target))[[1]]
        if(length(parts)!=3L) next
        drug <- gsub("\\\\'","'",parts[2])
        dose <- as.numeric(parts[3])
        ci <- which(conditions$drug==drug & abs(conditions$dose-dose)<1e-9)
        if(!length(ci)) next
        check(length(ci)==1L,'Ambiguous source condition')
        plate <- keys$plate[j]
        check(!plate %in% profile_plates[[ci]],'Duplicate cell/plate/condition in Parquet')
        treated_sum[ci,] <- treated_sum[ci,]+g$target_mean
        profile_plates[[ci]] <- c(profile_plates[[ci]],plate)
        key <- as.character(plate)
        if(is.null(reference_list[[key]])) reference_list[[key]] <- g$ref_mean else
          check(max(abs(reference_list[[key]]-g$ref_mean))<1e-6,'Inconsistent plate vehicle profile')
      }
    }
    observed <- which(lengths(profile_plates)>0)
    md <- copy(conditions[observed])
    md[, `:=`(depmap_id=id,cell_line=cvcl,observation_id=paste(id,perturbation_id,sep='__'),
      n_plates=lengths(profile_plates)[observed],contributing_plates=vapply(profile_plates[observed],
        function(z) paste(sort(z),collapse=','),character(1)),treated_total=NA_real_,reference_total=NA_real_,
      treated_correction_factor=NA_real_)]
    plate_ids <- as.integer(names(reference_list))
    ref_counts <- do.call(rbind,reference_list)
    response <- matrix(NA_real_,nrow(md),length(genes),dimnames=list(md$observation_id,genes))
    for(i in seq_len(nrow(md))) {
      ci <- observed[i]
      treated <- treated_sum[ci,]/length(profile_plates[[ci]])
      ref <- colMeans(ref_counts[match(profile_plates[[ci]],plate_ids),,drop=FALSE])
      factor <- composition_factor(treated,ref)
      response[i,] <- log2(1+treated[gene_index]/sum(treated)*10000/factor)-log2(1+ref[gene_index]/sum(ref)*10000)
      md[i,`:=`(treated_total=sum(treated),reference_total=sum(ref),treated_correction_factor=factor)]
    }
    rm(treated_sum,tab)
  }
  md[, `:=`(source_drug=drug,drug=panel_name(drug),recipient=recipient)]
  md[, focal := drug %chin% local_drugs]
  check(nrow(ref_counts)>=10L && all(is.finite(ref_counts)) && all(ref_counts>=0),'Invalid plate controls')
  vehicle <- ref_counts/rowSums(ref_counts)*10000
  baseline <- colMeans(log2(1+vehicle[,gene_index,drop=FALSE]))
  null <- matrix(NA_real_,nrow(ref_counts),length(genes),dimnames=list(paste(id,plate_ids,sep='__'),genes))
  for(i in seq_len(nrow(ref_counts))) {
    ref <- colMeans(ref_counts[-i,,drop=FALSE])
    factor <- composition_factor(ref_counts[i,],ref)
    null[i,] <- log2(1+ref_counts[i,gene_index]/sum(ref_counts[i,])*10000/factor)-log2(1+ref[gene_index]/sum(ref)*10000)
  }
  check(all(is.finite(response)) && !anyDuplicated(md$observation_id),'Invalid normalized expression')
  if(id %chin% names(legacy$null)) check(max(abs(null-legacy$null[[id]]))<1e-12,'Original control normalization changed')
  prepared <- list(response=response,null=setNames(list(null),id),baseline=setNames(list(baseline),id),metadata=md,
    genes=genes,normalizer='paired_trimmed_median_log_ratio_cp10k_v1')
  prepared$signature <- digest::digest(list(id,md,response,null),algo='sha256')
  qsave(prepared,path,preset='high')
  metadata_parts[[id]] <- md
  message('Molecular inputs ready: ',recipient,'; ',nrow(md),' conditions')
  rm(prepared,response,null,ref_counts)
  invisible(gc())
}
metadata <- rbindlist(metadata_parts,fill=TRUE)
sensitivity <- align_panel_outcomes(metadata,cell_map,audit_dir)
coverage <- sensitivity[,.(prism_conditions=sum(!is.na(prism_abundance_pct)),
  prism_drugs=uniqueN(drug[!is.na(prism_abundance_pct)]),almanac_conditions=sum(!is.na(almanac_abundance_pct)),
  almanac_drugs=uniqueN(drug[!is.na(almanac_abundance_pct)])),by=.(depmap_id,recipient)]
coverage[, included := (prism_conditions>=minimum_conditions & prism_drugs>=minimum_drugs) |
  (almanac_conditions>=minimum_conditions & almanac_drugs>=minimum_drugs)]
coverage[, criterion := sprintf('At least %d matched conditions from %d drugs in either assay',minimum_conditions,minimum_drugs)]
fwrite(coverage,file.path(audit_dir,'cell_assay_coverage.tsv'),sep='\t')
cell_map <- cell_map[depmap_id %chin% coverage[included==TRUE,depmap_id]]
metadata <- metadata[depmap_id %chin% cell_map$depmap_id]
sensitivity <- sensitivity[depmap_id %chin% cell_map$depmap_id]
fwrite(cell_map,file.path(output_dir,'data/panel_cell_mapping.tsv'),sep='\t')
fwrite(metadata,file.path(audit_dir,'molecular_metadata.tsv.gz'),sep='\t')
fwrite(sensitivity,file.path(audit_dir,'sensitivity_alignment.tsv'),sep='\t')
qsave(list(cell_map=cell_map,metadata=metadata,sensitivity=sensitivity,genes=genes,
  input_paths=file.path(input_dir,paste0(cell_map$depmap_id,'.qs')),
  coverage=coverage),file.path(input_dir,'panel_index.qs'),preset='high')
writeLines(capture.output(sessionInfo()),file.path(audit_dir,'sessionInfo.txt'))
message('Full panel prepared: ',nrow(cell_map),' eligible cell lines; ',nrow(metadata),' expression conditions')
