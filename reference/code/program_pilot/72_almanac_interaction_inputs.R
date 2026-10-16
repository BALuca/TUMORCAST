#!/usr/bin/env Rscript
# Reuse matched ALMANAC/Tahoe queries and score the frozen original40 bank.
suppressPackageStartupMessages({library(data.table); library(qs); library(Matrix); library(digest)})
root <- '/work/results/program_pilot/72_almanac_full_program_interactions'
prior <- '/work/results/program_pilot/21_almanac_two_partner_ward_programs'
weight_path <- '/work/results/program_pilot/61_individual_responses_combination_transfer/models/original40_weights.qs'
raw_path <- '/work/almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv'
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
for (folder in c('data','models','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
prepared <- qread(file.path(prior,'models/matched_combinations.qs'))
d <- as.data.table(prepared$combinations)
stopifnot(nrow(d)==18832L,uniqueN(d$pair)==679L,!anyDuplicated(d$combination_id))
raw <- fread(file.path(prior,'data/raw_bliss_reconstruction.tsv.gz'))
# QS and TSV preserve slightly different binary representations of decimal doses.
# Canonical dose strings identify the same experimental setting without nearest
# dose matching; require unique condition keys and complete reconstruction.
for (field in c('dose_a','dose_b')) {
  d[,(paste0(field,'_key')) := sprintf('%.10g',get(field))]
  raw[,(paste0(field,'_key')) := sprintf('%.10g',get(field))]
}
keys <- c('depmap_id','pair','dose_a_key','dose_b_key','SCREENER')
stopifnot(!anyDuplicated(d[,..keys]))
raw <- merge(raw,d[,c(keys,'combination_id'),with=FALSE],by=keys,all=FALSE,sort=FALSE)
stopifnot(!anyDuplicated(raw$COMBODRUGSEQ),uniqueN(raw$combination_id)==nrow(d))
raw[, log_eligible := survival>0 & survival_a>0 & survival_b>0]
fwrite(raw[,.(raw_records=.N,excluded_nonpositive=sum(!log_eligible)),by=SCREENER],
  file.path(root,'data/log_endpoint_exclusions.tsv'),sep='\t')
raw <- raw[log_eligible==TRUE]
raw[, observed_g := log2(survival/(survival_a*survival_b))]
outcomes <- raw[,.(observed_g=mean(observed_g),individual_log2=mean(log2(survival_a*survival_b)),
  observed_ratio=mean(survival),expected_ratio=mean(survival_a*survival_b),
  raw_records_log=.N,study_records=uniqueN(paste(STUDY,TESTDATE))),by=combination_id]
d <- merge(d,outcomes,by='combination_id',all=FALSE,sort=FALSE)
setorder(d,combination_id)
d[, row := .I]
stopifnot(all(is.finite(d$observed_g)))

# Recover original sample identifiers before canonicalizing the partner order.
identity <- fread(raw_path,select=c('COMBODRUGSEQ','NSC1','NSC2','SAMPLE1','SAMPLE2'))
identity <- identity[match(raw$COMBODRUGSEQ,COMBODRUGSEQ)]
stopifnot(identical(identity$COMBODRUGSEQ,raw$COMBODRUGSEQ))
sample_a <- ifelse(identity$NSC1==raw$NSC1,identity$SAMPLE1,identity$SAMPLE2)
sample_b <- ifelse(identity$NSC2==raw$NSC2,identity$SAMPLE2,identity$SAMPLE1)
stopifnot(!anyNA(c(sample_a,sample_b)))
context <- paste(raw$SCREENER,raw$STUDY,raw$TESTDATE,raw$depmap_id,sep='|')
refs <- unique(rbind(data.table(combination_id=raw$combination_id,
  ref=paste(context,raw$NSC1,sample_a,sprintf('%.12g',raw$dose_a),sep='|')),
  data.table(combination_id=raw$combination_id,
  ref=paste(context,raw$NSC2,sample_b,sprintf('%.12g',raw$dose_b),sep='|'))))
refs[, row := match(combination_id,d$combination_id)]
stopifnot(!anyNA(refs$row))
rm(identity,raw)

W <- qread(weight_path)
stopifnot(nrow(W)==8842L,ncol(W)==40L)
scores <- list()
inputs <- c(file.path(prior,'models/matched_combinations.qs'),
  file.path(prior,'data/raw_bliss_reconstruction.tsv.gz'),weight_path,raw_path)
for (recipient in unique(d$recipient)) {
  path <- file.path(prior,paste0('cell_',gsub('[^A-Za-z0-9]+','_',recipient)),'models/partner_queries.qs')
  stopifnot(file.exists(path))
  molecular <- qread(path)
  stopifnot(all(rownames(W) %in% colnames(molecular$displacement)),
    nrow(molecular$displacement)==nrow(molecular$queries))
  score <- as.matrix(molecular$displacement[,rownames(W),drop=FALSE] %*% W)
  rownames(score) <- as.character(molecular$queries$query_id)
  scores[[recipient]] <- score
  inputs <- c(inputs,path)
}
Q <- do.call(rbind,unname(scores))
stopifnot(!anyDuplicated(rownames(Q)),all(is.finite(Q)),
  all(as.character(c(d$query_a,d$query_b)) %in% rownames(Q)))
A <- Q[as.character(d$query_a),,drop=FALSE]
B <- Q[as.character(d$query_b),,drop=FALSE]
qsave(list(metadata=d,A=A,B=B,reference_map=refs,query_scores=Q,
  weight_sha256=digest(weight_path,file=TRUE,algo='sha256')),
  file.path(root,'data/prepared_inputs.qs'))
fwrite(d,file.path(root,'data/matched_conditions.tsv.gz'),sep='\t')
fwrite(refs,file.path(root,'data/single_reference_map.tsv.gz'),sep='\t')
fwrite(d[,.(conditions=.N,pairs=uniqueN(pair),drugs=uniqueN(c(drug_a,drug_b))),by=recipient],
  file.path(root,'data/coverage_by_cell.tsv'),sep='\t')
fwrite(data.table(path=inputs,sha256=vapply(inputs,digest,character(1),file=TRUE,algo='sha256')),
  file.path(root,'data/input_manifest.tsv'),sep='\t')
print(d[,.(conditions=.N,pairs=uniqueN(pair),cells=uniqueN(depmap_id),drugs=uniqueN(c(drug_a,drug_b)))])
capture.output(sessionInfo(),file=file.path(root,'preparation_sessionInfo.txt'))
writeLines('Inputs complete; all original40 genes measured; exact reference map saved.',file.path(root,'STATUS.txt'))
