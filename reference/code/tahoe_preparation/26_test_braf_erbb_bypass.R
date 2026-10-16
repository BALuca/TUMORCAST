#!/usr/bin/env Rscript
## Stage 26: BRAF--ERBB experiment. Source data are shared; experimental
## membership, drug matching and split are defined here, not by earlier pilots.
## User guidance: diagnose and adapt to observed failures. Retain the explicit
## exclusion of NCI-H2122 and NCI-H596; audit them without using them in analysis.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
MAX_GENES_PER_DRUG <- 40L
PROGRAM_GROUPS <- 6L
OUTLIER_IDS <- c('ACH-000311','ACH-000628')
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- if(length(script_arg)) dirname(dirname(normalizePath(sub('^--file=', '', script_arg[1])))) else normalizePath('/work/tahoe_prism_v2')
work <- dirname(root)
cache <- file.path(work, 'cache/tahoe_prism_full_chain')
pdex <- file.path(cache, 'tahoe_pdex_selected')
out <- file.path(root, 'results/26_braf_erbb_bypass')
art <- file.path(root, 'artifacts/26_braf_erbb_bypass')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(art, recursive=TRUE, showWarnings=FALSE)
source(file.path(work, 'src/R/lib_heatmap.R'))
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name), sep='\t', na='NA')
check <- function(ok, msg) if(!isTRUE(ok)) stop(msg, call.=FALSE)

## 1. Independent input audit and experimental design -------------------------
paths <- c(tahoe_metadata=file.path(cache,'response_effect_metadata.tsv'),
  tahoe_cells=file.path(pdex,'cells.tsv'), tahoe_genes=file.path(pdex,'genes.tsv'),
  tahoe_references=file.path(pdex,'pdex_plate_reference.npy'),
  tahoe_reference_metadata=file.path(pdex,'pdex_plate_reference_metadata.tsv'),
  tahoe_plate_audit=file.path(pdex,'pdex_plate_audit.tsv.gz'),
  prism_treatments=file.path(cache,'prism/secondary-screen-replicate-collapsed-treatment-info.csv'),
  prism_matrix=file.path(cache,'prism/secondary-screen-replicate-collapsed-logfold-change.csv'),
  prism_cells=file.path(cache,'prism/secondary-screen-cell-line-info.csv'),
  cell_metadata=file.path(work,'outputs/dataset_inventory/cache/depmap/Model.csv'),
  cytosig_metadata=file.path(work,'data/01_CytoSig_data/meta_info.txt'),
  cytosig_effects=file.path(work,'data/01_CytoSig_data/diff.merge.gz'))
check(all(file.exists(paths)), paste('Missing processed inputs:', paste(paths[!file.exists(paths)],collapse='; ')))
manifest <- data.table(input=names(paths),path=unname(paths),bytes=file.info(paths)$size)
write_table(manifest,'input_manifest.tsv')
## Explicit parent-compound mapping, including the Tahoe salt label. No global
## drug-name normalization and no prefiltered historical PRISM outcome table.
drug_map <- data.table(drug=c('Vemurafenib','Lapatinib ditosylate','Cobimetinib'),
  compound=c('vemurafenib','lapatinib','cobimetinib'),
  role=c('BRAF inhibitor','ERBB inhibitor','MEK inhibitor'))
meta <- merge(fread(paths[['tahoe_metadata']]),drug_map,by='drug')
meta[, observation_id:=paste(depmap_id,perturbation_id,sep='__')]
check(!anyDuplicated(meta$observation_id),'Duplicate Tahoe observations')
check(setequal(meta$dose,c(.05,.5,5)) && all(meta$dose_unit=='uM'),'Unexpected Tahoe doses or units')
check(all(meta$available==1L),'Missing Tahoe molecular observations: inspect metadata before proceeding')
cells <- fread(paths[['tahoe_cells']])
models <- fread(paths[['cell_metadata']])
check(!anyDuplicated(models$ModelID),'Ambiguous DepMap model IDs')
cells <- merge(cells,models[,.(depmap_id=ModelID,rrid=RRID,tissue=OncotreeLineage,
  disease=OncotreePrimaryDisease)],by='depmap_id',all.x=TRUE)
check(all(cells$cell_line==cells$rrid) && !anyNA(cells$tissue),'Ambiguous cell-line identities')
prism_cells <- fread(paths[['prism_cells']])
check(!anyDuplicated(prism_cells$row_name),'Ambiguous PRISM matrix row identities')
row_ids <- fread(paths[['prism_matrix']],select=1L)[[1]]
cells[, prism_row:=prism_cells$row_name[match(depmap_id,prism_cells$depmap_id)]]
cells[, prism_identity_confirmed:=prism_cells$passed_str_profiling[match(depmap_id,prism_cells$depmap_id)] %in% c(TRUE,'True','TRUE')]
check(all(cells$prism_row %in% row_ids) && all(cells$prism_identity_confirmed),'Missing or failed-identity PRISM line')
check(setequal(cells$depmap_id,unique(meta$depmap_id)),'Molecular/cell metadata mismatch')
## Genotype annotation is independently retrieved by stable Cellosaurus ID.
## An absent annotation is NOT relabelled as wild type.
dir.create(file.path(art,'cellosaurus'),showWarnings=FALSE)
cells[, `:=`(braf_v600_confirmed=FALSE,braf_annotation='No V600 annotation; not a wild-type call',braf_source='')]
for(i in seq_len(nrow(cells))) {
  u <- paste0('https://www.cellosaurus.org/',cells$cell_line[i],'.txt')
  f <- file.path(art,'cellosaurus',paste0(cells$cell_line[i],'.txt'))
  if(!file.exists(f)) download.file(u,f,mode='wb',quiet=TRUE)
  s <- readLines(f,warn=FALSE)
  check(any(grepl(paste0('AC   ',cells$cell_line[i]),s,fixed=TRUE)),paste('Invalid Cellosaurus record',f))
  v <- grep('Sequence variation:.*BRAF;',s,value=TRUE)
  confirmed <- any(grepl('p.Val600',v,fixed=TRUE))
  cells[i,`:=`(braf_v600_confirmed=confirmed,braf_annotation=if(confirmed) paste(v,collapse=' | ') else 'No V600 annotation; not a wild-type call',braf_source=u)]
}
## Preserve the scarce confirmed BRAF-V600 lines for the specified primary
## vemurafenib discovery. Draw a new seven-line holdout from other solid tumors.
## The seed and IDs are saved BEFORE reading PRISM values or expression.
## This evaluates transfer of the root across held-out solid-tumor lines;
## BRAF-mutant melanoma validation is provided separately by the external panel.
setorder(cells,depmap_id)
check(!any(cells$tissue %in% c('Myeloid','Lymphoid','Blood')),'Non-solid tumor requires explicit eligibility definition')
set.seed(68453L)
holdout_ids <- sort(sample(cells[braf_v600_confirmed==FALSE,depmap_id],7L))
cells[, split:=ifelse(depmap_id %in% holdout_ids,'Holdout','Development')]
write_table(cells,'experiment_cell_lines.tsv')
write_table(cells[,.(depmap_id,cell_name,tissue,braf_v600_confirmed,split)],'experiment_split.tsv')
write_table(data.table(seed=68453L,rule='Seven random non-V600-confirmed solid-tumor lines; retain all confirmed V600 lines for discovery',
  split_unit='Entire cell line: all drugs and doses together',n_development=sum(cells$split=='Development'),n_holdout=sum(cells$split=='Holdout'),
  response_values_used=FALSE,prior_split_used=FALSE,prior_excluded_ids_used=FALSE),'split_design.tsv')
meta <- merge(meta,cells[,.(depmap_id,cell_name,tissue,braf_v600_confirmed,split)],by='depmap_id')
setorder(meta,compound,depmap_id,dose)

## Match each drug directly to cached PRISM compound records. PRISM dose grids
## differ from Tahoe; interpolate signed log2 viable-abundance response over
## log10 dose at 0.05, 0.5 and 5 uM within each screen, then average screens.
## Require both bracketing measurements; do not extrapolate or fill missing data.
pt <- fread(paths[['prism_treatments']])
pt <- pt[tolower(name) %chin% drug_map$compound & is.finite(dose) & dose>0]
pt[, compound:=tolower(name)]
check(setequal(pt$compound,drug_map$compound),'Required PRISM compound missing')
check(!anyDuplicated(pt$column_name),'Duplicate PRISM treatment IDs')
write_table(merge(drug_map,unique(pt[,.(compound,prism_name=name,broad_id,screen_id)]),by='compound'),'drug_identity_mapping.tsv')
pm <- fread(paths[['prism_matrix']],select=c(1L,match(pt$column_name,names(fread(paths[['prism_matrix']],nrows=0L,check.names=FALSE)))),check.names=FALSE)
setnames(pm,1L,'prism_row')
pm <- pm[prism_row %chin% cells$prism_row]
long <- melt(pm,id.vars='prism_row',variable.name='column_name',value.name='prism_survival',variable.factor=FALSE)
long <- merge(long,pt[,.(column_name,compound,screen_id,dose,broad_id)],by='column_name')
long <- merge(long,cells[,.(prism_row,depmap_id)],by='prism_row')
prism_interp <- long[,{
  grid <- sort(unique(dose)); ys <- vapply(grid,function(x) {
    y <- prism_survival[dose==x]; if(any(is.finite(y))) mean(y[is.finite(y)]) else NA_real_
  },numeric(1))
  rbindlist(lapply(c(.05,.5,5),function(target) {
    lo <- max(grid[grid<=target]); hi <- min(grid[grid>=target]);
    yl <- ys[match(lo,grid)]; yh <- ys[match(hi,grid)]
    value <- if(all(is.finite(c(yl,yh)))) yl+(yh-yl)*(log10(target)-log10(lo))/(log10(hi)-log10(lo)) else NA_real_
    data.table(dose=target,lower_dose=lo,upper_dose=hi,lower_response=yl,upper_response=yh,prism_survival=value)
  }))
},by=.(depmap_id,compound,screen_id)]
write_table(prism_interp,'prism_dose_interpolation.tsv')
prism_response <- prism_interp[,.(n_screens=sum(is.finite(prism_survival)),
  prism_survival=if(any(is.finite(prism_survival))) mean(prism_survival[is.finite(prism_survival)]) else NA_real_),by=.(depmap_id,compound,dose)]
prism_response[, prism_sensitivity:=-prism_survival]
check(all(abs(prism_response$prism_survival+prism_response$prism_sensitivity)<1e-12,na.rm=TRUE),'PRISM sign mismatch')
meta <- merge(meta,prism_response,by=c('depmap_id','compound','dose'),all.x=TRUE)
meta[, dose_weight:=c(.25,.50,.25)[match(dose,c(.05,.5,5))]]
write_table(meta,'tahoe_prism_observations.tsv')
line_drug <- meta[,.(n_molecular_doses=sum(available==1L),n_prism_doses=sum(is.finite(prism_survival)),
  prism_survival_auc=if(all(is.finite(prism_survival))) sum(dose_weight*prism_survival) else NA_real_),
  by=.(depmap_id,cell_name,tissue,braf_v600_confirmed,split,compound)]
line_drug[, `:=`(prism_sensitivity_auc=-prism_survival_auc,eligible=n_molecular_doses==3L & n_prism_doses==3L)]
write_table(line_drug,'independent_drug_overlap.tsv')
counts <- line_drug[,.(n_lines=.N,n_complete=sum(eligible),n_confirmed_v600=sum(eligible & braf_v600_confirmed)),by=.(compound,split)]
write_table(counts,'input_coverage.tsv')
cytometa <- fread(paths[['cytosig_metadata']],skip=1L,header=FALSE,col.names=c('profile_id','treatment','condition','duration','dose'))
cytometa <- cytometa[treatment %chin% c('EGF','NRG1')]
con <- gzfile(paths[['cytosig_effects']],'rt'); cyto_header <- strsplit(readLines(con,1L),'\t',fixed=TRUE)[[1]];close(con)
check(all(cytometa$profile_id %in% cyto_header),'Missing ligand profile columns')
check(nrow(cytometa)>0L,'No EGF or NRG1 effects available')
check(!any(grepl('GSE68453|GSE127988',cytometa$profile_id)),'External validation profiles overlap discovery')
write_table(cytometa,'cytosig_ligand_profiles.tsv')

## Diagnostic 1: independent cohort coverage. No activity gate is loaded.
column_keys <- as.vector(t(outer(drug_map$compound,c('.05','.5','5'),paste,sep=' | ')))
meta[, condition:=paste(compound,format(dose,trim=TRUE),sep=' | ')]
colmeta <- unique(meta[,.(compound,dose)]);colmeta[, drug_index:=match(compound,drug_map$compound)];setorder(colmeta,drug_index,dose)
colmeta[, condition:=paste(compound,dose,sep=' | ')]
meta[, condition:=paste(compound,dose,sep=' | ')]
plot_cells <- copy(cells);plot_cells[, split:=factor(split,c('Development','Holdout'))];setorder(plot_cells,split,-braf_v600_confirmed,tissue,cell_name)
mat <- matrix(NA_real_,nrow(plot_cells),nrow(colmeta),dimnames=list(plot_cells$cell_name,colmeta$condition))
for(i in seq_len(nrow(meta))) mat[match(meta$depmap_id[i],plot_cells$depmap_id),match(meta$condition[i],colmeta$condition)] <- if(is.finite(meta$prism_survival[i])) 2 else 1
ca <- data.frame(Drug=colmeta$compound,Dose=paste0(colmeta$dose,' uM'),row.names=colmeta$condition)
ra <- data.frame(Split=as.character(plot_cells$split),BRAF=ifelse(plot_cells$braf_v600_confirmed,'V600 confirmed','Not V600-confirmed'),row.names=plot_cells$cell_name)
cairo_pdf(file.path(out,'26a_independent_input_coverage.pdf'),width=9.8,height=10)
h <- draw_heatmap(mat,scale_rows=FALSE,value_range=c(0,1,2),value_colors=c('#FFFFFF','#CCCCCC','#237B9F'),
  col_anno_df=ca,col_anno_cols=c('Drug','Dose'),row_anno_df=ra,row_anno_cols=c('Split','BRAF'),ann_row_args=list(annotation_name_rot=90),
  cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=TRUE,show_column_names=FALSE,
  row_names_gp=grid::gpar(fontsize=11),row_split=plot_cells$split,column_split=factor(colmeta$compound,levels=drug_map$compound),
  cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(3,'mm'),column_gap=grid::unit(3,'mm'),
  hm_name='Availability',heatmap_legend_param=list(at=c(0,1,2),labels=c('Missing','Tahoe only','Tahoe + PRISM')),
  title='Independent BRAF--ERBB experiment | all 35 cached lines')
h <- ComplexHeatmap::draw(h,heatmap_legend_side='bottom',annotation_legend_side='bottom')
check(identical(as.integer(unlist(ComplexHeatmap::row_order(h))),seq_len(nrow(mat))),'Coverage row ordering changed')
check(identical(as.integer(unlist(ComplexHeatmap::column_order(h))),seq_len(ncol(mat))),'Coverage column ordering changed')
invisible(dev.off())
p <- ggplot(counts,aes(x=compound,y=n_complete,fill=split))+geom_col(position='dodge',width=.7)+
  geom_text(aes(label=n_complete),position=position_dodge(.7),vjust=-.3,size=4)+
  scale_fill_manual(values=c(Development='#237B9F',Holdout='#D98A36'))+
  labs(x=NULL,y='Lines with all three molecular + PRISM doses',fill=NULL,title='New experiment: complete matched coverage')+
  theme_classic(base_size=12)+theme(legend.position='bottom')+expand_limits(y=max(counts$n_complete)+2)
ggsave(file.path(out,'26b_independent_cohort_counts.pdf'),p,width=7,height=4.6,device=cairo_pdf)
write_table(data.table(figure=c('26a_independent_input_coverage.pdf','26b_independent_cohort_counts.pdf'),
  description=c('All 35 lines; independent split; BRAF annotation; all three drugs and doses','Complete three-dose molecular/PRISM overlap by independent split')),'FIGURE_INDEX.tsv')
write_table(data.table(check=c('35 original processed cell lines considered','No prior line exclusions','No prior split read',
  'No prior activity gate read','No prior gene selection read','All three parent compounds matched','New split frozen before response values',
  'All doses kept together by cell line','PRISM survival sign retained','CytoSig EGF/NRG1 available'),
  passed=c(nrow(cells)==35L,all(c('ACH-000311','ACH-000628') %in% cells$depmap_id),TRUE,TRUE,TRUE,
    setequal(pt$compound,drug_map$compound),TRUE,!anyDuplicated(cells$depmap_id),TRUE,nrow(cytometa)>0)), 'independence_checks.tsv')
message('Independent input audit complete: ',nrow(cells),' lines; ',sum(cells$split=='Development'),' development; ',sum(cells$split=='Holdout'),' holdout. No prior pilot gates loaded.')

## Download/cache the requested public external files independently of any
## previous pilot. Numerical expression tables remain unused during discovery.
external_urls <- c(
  GSE68453_family.soft.gz='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68453/soft/GSE68453_family.soft.gz',
  GSE68453_series_matrix.txt.gz='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68453/matrix/GSE68453_series_matrix.txt.gz',
  GSE68453_record.txt='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68453&targ=self&form=text&view=full',
  GSE68453_S1.eps='https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0138210.s001',
  GSE68453_study.html='https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0138210',
  GSE127988_study.html='https://pmc.ncbi.nlm.nih.gov/articles/PMC8009031/?pdf=0')
for(n in names(external_urls)) if(!file.exists(file.path(art,n)))
  download.file(external_urls[[n]],file.path(art,n),mode='wb',quiet=TRUE)
write_table(data.table(file=names(external_urls),url=unname(external_urls),
  path=file.path(art,names(external_urls)),bytes=file.info(file.path(art,names(external_urls)))$size),'external_download_manifest.tsv')
if(!file.exists(file.path(art,'GSE68453_soft_metadata.txt'))) {
  con <- gzfile(file.path(art,'GSE68453_family.soft.gz'),'rt')
  metadata_lines <- character()
  repeat {
    chunk <- readLines(con,n=10000L,warn=FALSE)
    if(!length(chunk)) break
    metadata_lines <- c(metadata_lines,grep('^[!^]',chunk,value=TRUE))
  }
  close(con)
  writeLines(metadata_lines,file.path(art,'GSE68453_soft_metadata.txt'))
}

## External design audit reads metadata only, not expression values or Bliss.
get_soft_samples <- function(s) {
  starts <- grep('^\\^SAMPLE = ',s)
  rbindlist(lapply(seq_along(starts),function(i) {
    block <- s[starts[i]:(if(i==length(starts)) length(s) else starts[i+1]-1L)]
    field <- function(prefix) {
      x <- block[startsWith(block,prefix)]
      check(length(x)==1L,paste('Ambiguous SOFT field',prefix))
      substring(x,nchar(prefix)+1L)
    }
    data.table(gsm=field('^SAMPLE = '),title=field('!Sample_title = '),
      characteristics=paste(block[startsWith(block,'!Sample_characteristics_ch1 = ')],collapse=' | '),
      protocol=paste(block[startsWith(block,'!Sample_treatment_protocol_ch1 = ')],collapse=' | '))
  }))
}
f684 <- file.path(art,'GSE68453_soft_metadata.txt')
check(file.exists(f684),'Missing downloaded GSE68453 SOFT metadata')
s684 <- get_soft_samples(readLines(f684,warn=FALSE))
s684[, `:=`(cell_line=sub(' .*','',title),arm=sub('^[^ ]+ (.*) rep [0-9]+$','\\1',title),
  replicate=as.integer(sub('.* rep ([0-9]+)$','\\1',title)))]
factorial_lines <- c('A375','DM331','HT144','SKMEL24','SKMEL28')
external_lines <- c('A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17')
s684 <- s684[cell_line %chin% external_lines & replicate %in% 1:3 &
  (arm %chin% c('control','plx') | (cell_line %chin% factorial_lines & arm %chin% c('lapatinib','plx_lap')))]
check(nrow(s684)==102L && !anyDuplicated(s684[,.(cell_line,arm,replicate)]),'GSE68453 selected design incomplete')
check(all(grepl('8 hours',s684$protocol,fixed=TRUE)),'Unexpected GSE68453 sampling time')
check(all(grepl('V600E',s684$characteristics,fixed=TRUE)),'GSE68453 BRAF genotype mismatch')
check(all(grepl('125 nM',s684[arm %chin% c('plx','plx_lap'),protocol],fixed=TRUE)), 'Unexpected PLX dose')
check(all(grepl('2 µM lapatinib',s684[arm %chin% c('lapatinib','plx_lap'),protocol],fixed=TRUE)), 'Unexpected lapatinib dose')
write_table(s684,'gse68453_sample_audit.tsv')
c684 <- s684[,.(n=.N,replicates=paste(sort(replicate),collapse=',')),by=.(cell_line,arm)]
check(all(c684$n==3L),'GSE68453 replicate counts mismatch')
write_table(c684,'gse68453_condition_counts.tsv')
geo127 <- file.path(art,'gse127988');dir.create(geo127,showWarnings=FALSE)
urls127 <- c(GSE127988_family.soft.gz='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/soft/GSE127988_family.soft.gz',
  GSE127988_RPM_Gerosa_et_al.csv.gz='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/suppl/GSE127988_RPM_Gerosa_et_al.csv.gz')
for(n in names(urls127)) if(!file.exists(file.path(geo127,n))) download.file(urls127[[n]],file.path(geo127,n),mode='wb',quiet=TRUE)
s127 <- get_soft_samples(readLines(gzfile(file.path(geo127,'GSE127988_family.soft.gz')),warn=FALSE))
number <- function(x,label) as.numeric(sub(paste0('.*',label,': ([0-9.]+).*'),'\\1',x))
s127[, `:=`(vem_uM=number(characteristics,'vemurafenib for 24 hrs \\(um\\)'),
  cobi_uM=number(characteristics,'cobimetinib for 24 hrs \\(um\\)'),
  egf_h=number(characteristics,'egf at 100 ng/ml \\(h\\)'),biological_set=number(characteristics,'set id'))]
s127 <- s127[vem_uM==1 & biological_set %in% 3:4]
s127[, `:=`(cell_line='A375',arm=ifelse(cobi_uM==0,'Vem','Vem + Cobi'))]
check(nrow(s127)==56L && setequal(s127$cobi_uM,c(0,1)) && setequal(s127$egf_h,c(0,.5,1,2,3,4,8)), 'GSE127988 required time-course design incomplete')
c127 <- s127[,.(n=.N,n_biological_days=uniqueN(biological_set)),by=.(arm,egf_h)]
check(all(c127$n==4L & c127$n_biological_days==2L),'GSE127988 replicate/day counts mismatch')
write_table(s127,'gse127988_sample_audit.tsv');write_table(c127,'gse127988_condition_counts.tsv')
cairo_pdf(file.path(out,'26b2_external_design_counts.pdf'),width=9,height=6)
p <- ggplot(c684,aes(arm,n,fill=arm))+geom_col()+facet_wrap(~cell_line,ncol=4)+
  scale_y_continuous(breaks=0:3,limits=c(0,3.5))+theme_classic(base_size=11)+theme(legend.position='none',axis.text.x=element_text(angle=60,hjust=1))+
  labs(x=NULL,y='RNA replicates',title='GSE68453 | 8 h | 125 nM PLX4720; 2 uM lapatinib')
print(p)
p <- ggplot(c127,aes(factor(egf_h),n,fill=arm))+geom_col()+facet_wrap(~arm)+theme_classic(base_size=12)+
  scale_y_continuous(breaks=0:4,limits=c(0,4.5))+theme(legend.position='none')+
  labs(x='EGF exposure (h)',y='RNA wells',title='GSE127988 | two drug backgrounds',subtitle='4 wells per condition: 2 biological days x 2 wells')
print(p);invisible(dev.off())

## 2. Dose-integrated responses from saved pseudobulk products ----------------
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
## Apply the existing paired, sample-local composition normalization to the
## saved processed pseudobulks. All 62,710 genes remain in the expression output.
## This is a measurement transformation, not a gene or cell-line eligibility gate.
normalize_pair <- function(treated,reference) {
  tol <- max(1e-6,32*2^-23*max(abs(c(treated,reference)),1))
  check(min(c(treated,reference))>=-tol,'Material negative reconstructed abundance')
  treated <- pmax(treated,0);reference <- pmax(reference,0)
  tt <- sum(treated); rt <- sum(reference)
  check(tt>0 && rt>0,'Empty pseudobulk')
  tp <- treated/tt; rp <- reference/rt
  eligible_scale <- tp*10000>=.1 & rp*10000>=.1
  m <- log2(tp[eligible_scale]/rp[eligible_scale]);a <- .5*log2(tp[eligible_scale]*rp[eligible_scale])
  mq <- quantile(m,c(.3,.7)); aq <- quantile(a,c(.05,.95))
  kept <- m>=mq[1]&m<=mq[2]&a>=aq[1]&a<=aq[2]
  check(sum(kept)>0,'Cannot estimate sample-local composition scaling')
  factor <- 2^median(m[kept])
  list(effect=log2(1+tp*10000/factor)-log2(1+rp*10000),
    treated_total=tt,reference_total=rt,factor=factor)
}
genes <- fread(paths[['tahoe_genes']]);setorder(genes,gene_index)
check(!anyDuplicated(genes$feature),'Ambiguous Tahoe gene IDs')
plates <- fread(paths[['tahoe_plate_audit']])
plates[, observation_id:=paste(depmap_id,perturbation_id,sep='__')]
plates <- plates[observation_id %chin% meta$observation_id]
refmeta <- fread(paths[['tahoe_reference_metadata']]);refmeta[, key:=paste(depmap_id,plate,sep='__')]
refs <- read_npy_rows(paths[['tahoe_references']],refmeta$array_row)
rownames(refs) <- refmeta$key
setorder(meta,compound,depmap_id,dose)
effects <- matrix(NA_real_,nrow(meta),nrow(genes),dimnames=list(meta$observation_id,genes$feature))
normalization <- vector('list',nrow(meta))
reference_cp10k <- effects
treated_cp10k <- effects
for(path in unique(meta$shard_path)) {
  ids <- which(meta$shard_path==path)
  raw <- read_npy_rows(file.path(cache,path),meta$shard_row[ids])
  for(j in seq_along(ids)) {
    i <- ids[j]; oid <- meta$observation_id[i]
    pp <- plates[observation_id==oid,plate]
    check(length(pp)==meta$n_plates[i] && !anyDuplicated(pp),paste('Plate mismatch:',oid))
    keys <- paste(meta$depmap_id[i],pp,sep='__')
    check(all(keys %in% rownames(refs)),paste('Missing reference:',oid))
    rr <- colMeans(refs[keys,,drop=FALSE]);z <- normalize_pair(raw[j,]+rr,rr)
    effects[i,] <- z$effect
    reference_cp10k[i,] <- pmax(rr,0)/z$reference_total*10000
    treated_cp10k[i,] <- pmax(raw[j,]+rr,0)/z$treated_total*10000/z$factor
    normalization[[i]] <- data.table(observation_id=oid,treated_total=z$treated_total,reference_total=z$reference_total,
      composition_factor=z$factor,contributing_plates=paste(pp,collapse=','))
  }
}
rm(refs,raw)
write_table(rbindlist(normalization),'sample_normalization.tsv')
line_drug[, key:=paste(depmap_id,compound,sep='__')]
setorder(line_drug,compound,depmap_id)
auc <- matrix(NA_real_,nrow(line_drug),ncol(effects),dimnames=list(line_drug$key,colnames(effects)))
for(i in seq_len(nrow(line_drug))) {
  ids <- which(meta$depmap_id==line_drug$depmap_id[i]&meta$compound==line_drug$compound[i])
  check(length(ids)==3 && setequal(meta$dose[ids],c(.05,.5,5)),'Incomplete molecular AUC')
  auc[i,] <- colSums(effects[ids,,drop=FALSE]*meta$dose_weight[ids])
}
write_table(line_drug,'dose_integrated_response_metadata.tsv')
## Dose diagnostics summarize each gene across development lines, without
## choosing one dose or filtering on the resulting plot.
dose_medians <- sapply(seq_len(nrow(colmeta)),function(j) {
  ix <- which(meta$compound==colmeta$compound[j]&meta$dose==colmeta$dose[j]&meta$split=='Development' & !meta$depmap_id %chin% OUTLIER_IDS)
  apply(effects[ix,,drop=FALSE],2,median)
})
rownames(dose_medians) <- colnames(effects);colnames(dose_medians) <- colmeta$condition
write_table(data.table(gene=rownames(dose_medians),dose_medians),'dose_gene_medians.tsv.gz')
cd <- cor(dose_medians,method='spearman')
cairo_pdf(file.path(out,'26c_drug_dose_correlations.pdf'),width=8.4,height=7.6)
h <- draw_heatmap(cd,scale_rows=FALSE,value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),
  show_row_names=TRUE,show_column_names=TRUE,cluster_rows=FALSE,cluster_columns=FALSE,
  row_names_gp=grid::gpar(fontsize=10),column_names_gp=grid::gpar(fontsize=10),
  hm_name='Spearman',title='Development median effects | all three doses')
ComplexHeatmap::draw(h);invisible(dev.off())

## Extract all exact EGF/NRG1 molecular profiles for the shared gene universe.
cyto <- fread(paths[['cytosig_effects']],select=c(1L,match(cytometa$profile_id,cyto_header)+1L),
  header=TRUE,col.names=c('gene',cytometa$profile_id),check.names=FALSE,showProgress=FALSE)
check(!anyDuplicated(cyto$gene),'Ambiguous CytoSig genes')
ligands <- unique(cytometa$treatment)
centroids <- sapply(ligands,function(lig) {
  xx <- as.matrix(cyto[,cytometa[treatment==lig,profile_id],with=FALSE])
  apply(xx,1,function(v) if(any(is.finite(v))) median(v[is.finite(v)]) else NA_real_)
})
rownames(centroids) <- cyto$gene
ligand_z <- scale(centroids)
erbb <- rowMeans(ligand_z,na.rm=TRUE); names(erbb) <- rownames(centroids)
## A separate uncentered, SD-scaled ligand coordinate preserves treatment-minus-control zero for biological orientation only.
erbb_signed <- rowMeans(sweep(centroids,2,attr(ligand_z,'scaled:scale'),'/'),na.rm=TRUE);names(erbb_signed) <- rownames(centroids)
universe <- intersect(colnames(auc),names(erbb)[is.finite(erbb)])
write_table(data.table(gene=names(erbb),centroids,erbb_ligand=erbb,erbb_signed_for_orientation=erbb_signed),'cytosig_ligand_centroids.tsv.gz')
# Reusable measurement cache contains no selected genes or fitted programs.
saveRDS(c(mget(c('effects','auc','meta','line_drug','cells','genes','cyto','cytometa',
  'centroids','erbb','universe','reference_cp10k','treated_cp10k')),list(prism_raw=long)),
  file.path(art,'discovery_measurements.rds'))
## 3. Diagnose measurement failures before selecting genes ------------------
## The brief is a starting policy. Revisions below follow observed failures;
## external values and held-out outcomes cannot change discovery decisions.
library(cluster)
rank_cor <- function(x,y) suppressWarnings(as.vector(cor(
  matrixStats::colRanks(x,ties.method='average',preserveShape=TRUE),rank(y))))
excluded_ids <- OUTLIER_IDS
normalization <- rbindlist(normalization)
qc <- merge(meta[,.(observation_id,depmap_id,cell_name,compound,dose,split)],normalization,by='observation_id')
qc[, detected_fraction:=rowMeans(pmax(reference_cp10k[observation_id,universe,drop=FALSE],
  treated_cp10k[observation_id,universe,drop=FALSE])>0)]
qc[, excluded:=depmap_id %chin% excluded_ids]
write_table(qc,'sample_quality_audit.tsv')
line_qc <- qc[,.(median_detected_fraction=median(detected_fraction),maximum_scaling_factor=max(composition_factor),
  excluded=any(excluded)),by=.(depmap_id,cell_name,split)]
line_qc[, reason:=ifelse(excluded,'User-confirmed outlier; sparse measurements and abnormal composition scaling','Retained')]
write_table(line_qc,'cell_line_decisions.tsv')
meta[, analysis_included:=!depmap_id %chin% excluded_ids]
line_drug[, analysis_included:=!depmap_id %chin% excluded_ids]
write_table(line_drug[,.(n_lines=.N,n_complete=sum(eligible)),by=.(compound,split,analysis_included)],'analysis_cohort_counts.tsv')
plot_save <- function(p,name,w=7,h=5) {
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=140)
}
p <- ggplot(qc,aes(detected_fraction,composition_factor,color=excluded))+
  geom_point(size=1.6,alpha=.6)+geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='grey50')+
  geom_text(data=qc[composition_factor>1.7],aes(label=cell_name),size=3,hjust=-.05,check_overlap=TRUE)+
  scale_x_continuous(expand=expansion(mult=c(.05,.3)))+scale_color_manual(values=c('FALSE'='#287D99','TRUE'='#B64039'))+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+labs(x='Fraction of shared genes with nonzero abundance',y='Composition scaling factor',color='Excluded',
  title='Measurement QC | all lines audited; two outliers excluded')
plot_save(p,'26a_sample_quality',7.8,6)
## Study x ligand medians prevent 15 NRG1 profiles from counting as 15 studies.
cytometa[, study:=sub('.*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*','\\1',profile_id,perl=TRUE)]
check(all(grepl('^(GSE[0-9]+|E-MTAB-[0-9]+)$',cytometa$study)),'Unparsed CytoSig study ID')
cytometa[, study_ligand:=paste(study,treatment,sep=' | ')]
cyto_full <- as.matrix(cyto[,-1,with=FALSE]);rownames(cyto_full)<-cyto$gene
study_meta <- unique(cytometa[,.(study_ligand,study,treatment)])
study_effect <- sapply(study_meta$study_ligand,function(k)
  matrixStats::rowMedians(cyto_full[,cytometa[study_ligand==k,profile_id],drop=FALSE],na.rm=TRUE))
rownames(study_effect)<-cyto$gene
study_n <- rowSums(is.finite(study_effect))
universe <- intersect(universe,names(study_n)[study_n>=ceiling(ncol(study_effect)/2)])
write_table(cytometa[,.(n_profiles=.N),by=.(study,treatment)],'cytosig_study_balance.tsv')
write_table(data.table(gene=rownames(study_effect),n_studies=study_n),'cytosig_gene_coverage.tsv')
## Quantify the prior selection's failure, but never use it for admission.
old <- fread(file.path(art,'archive/fixed_brief/results/selected_genes_by_drug.tsv'))
old_audit <- rbindlist(lapply(drug_map$compound,function(drug) {
  ids<-which(meta$compound==drug & meta$split=='Development')
  g<-old[compound==drug,gene]
  ab<-pmax(reference_cp10k[ids,g,drop=FALSE],treated_cp10k[ids,g,drop=FALSE])
  data.table(compound=drug,gene=g,fraction_above_01=colMeans(ab>=.1))
}))
write_table(old_audit,'previous_gene_measurement_audit.tsv')

## 4. Function-linked selection with estimability and influence checks -------
## V600 n=3 is a context annotation, not a gene-ranking sample. Primary gene
## discovery now estimates response associations across retained solid tumors.
## This changes the estimand; it cannot establish BRAF-mutant specificity.
selection <- list();influence <- list();dose_qc <- list()
for(drug in drug_map$compound) {
  ix<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$eligible & line_drug$analysis_included)
  x<-auc[ix,universe,drop=FALSE];y<-line_drug$prism_sensitivity_auc[ix]
  rr<-rank_cor(x,y)
  loo<-sapply(seq_along(ix),function(j)rank_cor(x[-j,,drop=FALSE],y[-j]))
  tissues<-line_drug$tissue[ix]
  lot<-sapply(unique(tissues),function(t){keep<-tissues!=t;rank_cor(x[keep,,drop=FALSE],y[keep])})
  dd<-sapply(c(.05,.5,5),function(d) {
    ii<-vapply(line_drug$depmap_id[ix],function(id)which(meta$depmap_id==id & meta$compound==drug & meta$dose==d),integer(1))
    rank_cor(effects[ii,universe,drop=FALSE],meta$prism_sensitivity[ii])
  })
  ids<-which(meta$compound==drug & meta$depmap_id %chin% line_drug$depmap_id[ix])
  ab<-pmax(reference_cp10k[ids,universe,drop=FALSE],treated_cp10k[ids,universe,drop=FALSE])
  pp<-2*pt(-abs(rr)*sqrt((length(ix)-2)/pmax(1-rr*rr,1e-15)),df=length(ix)-2)
  z<-data.table(compound=drug,gene=universe,n_lines=length(ix),rho_auc=rr,nominal_p_approx=pp,
    loo_min_signed=apply(loo*sign(rr),1,min),leave_tissue_min_signed=apply(lot*sign(rr),1,min),
    expressed_fraction=colMeans(ab>=.1),median_abundance=matrixStats::colMedians(ab),
    median_abs_effect=matrixStats::colMedians(abs(x)),
    rho_005=dd[,1],rho_05=dd[,2],rho_5=dd[,3],consistent_doses=rowSums(sign(dd)==sign(rr)))
  z[, measurable:=expressed_fraction>=.5]
  z[, stable:=is.finite(loo_min_signed)&loo_min_signed>=.3 & leave_tissue_min_signed>0]
  z[, pass:=measurable & stable & nominal_p_approx<.05 & consistent_doses>=2]
  z[is.na(pass),pass:=FALSE]
  setorder(z,-loo_min_signed,-leave_tissue_min_signed,gene)
  z[, admission_rank:=cumsum(pass)];z[,selected:=pass & admission_rank<=MAX_GENES_PER_DRUG]
  selection[[drug]]<-z
  # Same all-solid-tumor estimand on both sides: isolates outlier influence.
  ai<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$eligible)
  ar<-rank_cor(auc[ai,universe,drop=FALSE],line_drug$prism_sensitivity_auc[ai])
  influence[[drug]]<-data.table(compound=drug,gene=universe,rho_with_outliers=ar,rho_without_outliers=rr)
  pairs<-combn(c(.05,.5,5),2)
  dose_qc[[drug]]<-rbindlist(lapply(seq_along(ix),function(j) {
    ii<-vapply(c(.05,.5,5),function(d)which(meta$depmap_id==line_drug$depmap_id[ix[j]] & meta$compound==drug & meta$dose==d),integer(1))
    cc<-cor(t(effects[ii,universe,drop=FALSE]),method='spearman')
    data.table(compound=drug,cell_name=line_drug$cell_name[ix[j]],lower_dose=pairs[1,],upper_dose=pairs[2,],rho=cc[upper.tri(cc)])
  }))
}
tests<-rbindlist(selection);selected<-tests[selected==TRUE];chosen<-sort(unique(selected$gene))
check(length(chosen)>=2,'No estimable response-associated gene set; inspect selection audit')
write_table(tests,'gene_selection_tests.tsv.gz');write_table(selected,'selected_genes_by_drug.tsv')
write_table(tests[,.(n_lines=unique(n_lines),n_measured=sum(measurable),n_pass=sum(pass),n_selected=sum(selected)),by=compound],'gene_selection_summary.tsv')
influence<-rbindlist(influence);influence<-merge(influence,tests[,.(compound,gene,selected)],by=c('compound','gene'))
write_table(influence,'outlier_influence.tsv.gz');write_table(rbindlist(dose_qc),'dose_concordance_by_line.tsv')
write_table(influence[,.(rank_correlation=cor(rho_with_outliers,rho_without_outliers,method='spearman',use='complete.obs'),
  selected_max_change=max(abs(rho_with_outliers[selected]-rho_without_outliers[selected]))),by=compound],'outlier_influence_summary.tsv')
p<-ggplot(influence[is.finite(rho_with_outliers)&is.finite(rho_without_outliers)],aes(rho_with_outliers,rho_without_outliers,color=selected))+
  geom_point(size=.3,alpha=.15)+geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  facet_wrap(~compound,scales='free')+theme_classic(base_size=11)+theme(aspect.ratio=1)+
  scale_color_manual(values=c('FALSE'='grey65','TRUE'='#B64039'))+
  labs(x='Spearman rho including the two outliers',y='Spearman rho excluding the two outliers',title='Outlier influence | identical gene universe and estimand')
plot_save(p,'26b_outlier_influence',11,4.8)
p<-ggplot(tests[is.finite(rho_auc)],aes(log10(median_abundance+1e-6),abs(rho_auc),color=selected))+
  geom_point(size=.3,alpha=.13)+geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  facet_wrap(~compound,scales='free')+geom_vline(xintercept=-1,lty=2,color='grey50')+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+scale_color_manual(values=c('FALSE'='grey65','TRUE'='#B64039'))+
  labs(x='log10 median abundance (counts per 10,000)',y='Absolute development Spearman rho',title='Gene selection | measurement precedes association ranking')
plot_save(p,'26c_gene_measurement_selection',11,4.8)
p<-ggplot(rbindlist(dose_qc),aes(factor(lower_dose),rho,color=factor(upper_dose)))+geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(width=.08,height=0,seed=26),size=1)+facet_wrap(~compound)+theme_classic(base_size=11)+
  labs(x='Lower dose (uM)',y='Gene-effect Spearman rho within a line',color='Upper dose',title='Dose agreement after outlier removal | all doses retained')
plot_save(p,'26d_dose_concordance',10,4.5)

## 5. Real sample structure, with equal weight for the two sources -----------
## No four-coordinate rank correlation; no forced three-group cut.
ti<-which(line_drug$split=='Development' & line_drug$analysis_included)
tahoe_x<-t(auc[ti,chosen,drop=FALSE]);cyto_x<-study_effect[chosen,,drop=FALSE]
ct<-cor(t(tahoe_x),method='spearman');cc<-cor(t(cyto_x),method='spearman',use='pairwise.complete.obs')
n_pair<-tcrossprod(is.finite(cyto_x)*1L)
check(all(is.finite(ct))&&all(is.finite(cc)),'Non-estimable source correlations; inspect gene coverage')
correlation<-(ct+cc)/2;distance<-as.dist(pmax(1-correlation,0));tree<-hclust(distance,method='average')
curve<-rbindlist(lapply(2:min(10,length(chosen)-1),function(k) {
  g<-cutree(tree,k);s<-silhouette(g,distance)
  data.table(k=k,mean_silhouette=mean(s[,'sil_width']),negative_fraction=mean(s[,'sil_width']<0),minimum_size=min(table(g)))
}))
k<-PROGRAM_GROUPS
check(k<nrow(correlation),'Too few genes for the requested descriptive resolution')
# Six groups provide the finer requested resolution. Small groups and weak
# coherence are reported, never hidden by selecting a coarser cut.
write_table(data.table(k=k,minimum_size=curve$minimum_size[match(k,curve$k)],
  negative_fraction=curve$negative_fraction[match(k,curve$k)],
  purpose='Requested finer descriptive resolution; silhouette is diagnostic'), 'cluster_resolution_decision.tsv')
groups<-cutree(tree,k);leaf_order<-chosen[tree$order]
old_ids<-unique(groups[leaf_order]);groups<-setNames(match(groups,old_ids),chosen)
programs<-paste0('P',seq_len(k));membership<-data.table(gene=leaf_order,program=paste0('P',groups[leaf_order]),plotting_order=seq_along(leaf_order))
colors<-setNames(grDevices::hcl.colors(k,'Dark 3'),programs)
write_table(curve,'silhouette_by_k.tsv')
write_table(data.table(gene=chosen,cluster=groups,silhouette=silhouette(groups,distance)[,'sil_width']),'gene_silhouettes.tsv')
write_table(data.table(merge_step=seq_along(tree$height),left=tree$merge[,1],right=tree$merge[,2],height=tree$height),'discovery_tree_merges.tsv')
for(n in c('ct','cc','correlation'))write_table(data.table(gene=chosen,get(n)),paste0(n,'_gene_correlations.tsv.gz'))
write_table(data.table(gene=chosen,n_pair),'cytosig_pairwise_study_counts.tsv.gz')
coherence<-rbindlist(lapply(programs,function(p) {
  g<-membership[program==p,gene]
  rbindlist(lapply(c('Tahoe','CytoSig'),function(src) {
    m<-if(src=='Tahoe')ct[g,g] else cc[g,g];v<-m[upper.tri(m)]
    data.table(program=p,source=src,n_genes=length(g),mean_rho=mean(v),median_rho=median(v),positive_fraction=mean(v>0))
  }))
}))
write_table(coherence,'program_source_coherence.tsv')
selected_k<-k
p<-ggplot(curve,aes(k,mean_silhouette))+geom_line()+geom_point(size=2)+
  geom_point(data=curve[curve$k==selected_k],color='#B64039',size=3)+theme_classic(base_size=12)+
  labs(x='Number of groups',y='Mean silhouette',title=paste('Observed sample structure | selected k =',k),subtitle='Six descriptive groups; silhouette diagnoses separation, not a resolution veto')
# Avoid data.table name shadowing when highlighting the selected cut.
p$layers[[3]]$data<-curve[curve$k==selected_k]
plot_save(p,'26e_rank_silhouette',6,4.5)
## Shared heatmap helper: raw correlations are never rescaled to exaggerate them.
## A shared symmetric 95th-percentile OFF-DIAGONAL limit exposes weak structure;
## clipping fraction and full matrices are saved. Self-correlations are masked.
cor_limit<-as.numeric(quantile(abs(c(ct[upper.tri(ct)],cc[upper.tri(cc)])),.95))
write_table(data.table(limit=cor_limit,tahoe_clipped=mean(abs(ct[upper.tri(ct)])>cor_limit),cytosig_clipped=mean(abs(cc[upper.tri(cc)])>cor_limit)), 'correlation_display_scale.tsv')
heat <- function(x,name,title,rows=NULL,cols=NULL,limit=2.5,width=10,height=NULL,show_cols=TRUE) {
  if(is.null(height))height<-max(5,nrow(x)*.14+2)
  if(!is.null(rows))check(length(rows)==nrow(x),'Heatmap row annotations misaligned')
  if(!is.null(cols))check(length(cols)==ncol(x),'Heatmap column annotations misaligned')
  h<-draw_heatmap(x,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    cluster_rows=FALSE,cluster_columns=FALSE,row_split=rows,column_split=cols,
    cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(2,'mm'),column_gap=grid::unit(2,'mm'),
    show_row_names=TRUE,show_column_names=show_cols,row_names_gp=grid::gpar(fontsize=9),column_names_gp=grid::gpar(fontsize=8),
    column_names_rot=90,hm_name='Value',heatmap_legend_param=list(at=c(-limit,0,limit),labels=round(c(-limit,0,limit),2)),title=title)
  cairo_pdf(file.path(out,paste0(name,'.pdf')),width=width,height=height)
  drawn<-ComplexHeatmap::draw(h,heatmap_legend_side='bottom')
  check(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn))),seq_len(nrow(x))),'Heatmap changed frozen gene order')
  check(identical(as.integer(unlist(ComplexHeatmap::column_order(drawn))),seq_len(ncol(x))),'Heatmap changed frozen sample order')
  invisible(dev.off())
  png(file.path(out,paste0(name,'.png')),width=width,height=height,units='in',res=120,type='cairo')
  ComplexHeatmap::draw(h,heatmap_legend_side='bottom');invisible(dev.off())
}
row_blocks<-factor(membership$program,levels=programs)
for(src in c('Tahoe','CytoSig','Joint')) {
  m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=correlation)[leaf_order,leaf_order];diag(m)<-NA_real_
  heat(m,paste0('26f_',tolower(src),'_gene_correlations'),paste(src,'| same genes, same order | actual sample correlations'),row_blocks,row_blocks,cor_limit,11)
}
# Column clustering is descriptive; no numeric values are rearranged within rows.
source_heat <- function(x,name,title) {
  x<-x[leaf_order,,drop=FALSE]
  c<-cor(x,method='spearman',use='pairwise.complete.obs');check(all(is.finite(c)),'Non-estimable sample correlation')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average');ord<-tr$order
  # Preserve all gene rows. Standardization is for display only, centered within source.
  z<-t(scale(t(x)));z[!is.finite(z)]<-NA_real_
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord)),paste0(name,'_sample_order.tsv'))
  heat(z[,ord,drop=FALSE],name,title,row_blocks,limit=2.5,width=12)
  invisible(ord)
}
source_heat(tahoe_x,'26g_tahoe_gene_samples','Tahoe | 26 development lines x 3 drugs | frozen genes')
source_heat(cyto_x,'26g_cytosig_gene_samples','CytoSig | study-balanced EGF/NRG1 | frozen genes')
cairo_pdf(file.path(out,'26h_gene_hierarchy.pdf'),width=11,height=6)
par(mar=c(7,5,3,1));plot(tree,labels=chosen,hang=-1,cex=.8,main='Observed gene hierarchy | no causal or upstream assignment',xlab='',sub='Average Tahoe/CytoSig sample correlations');rect.hclust(tree,k=k);invisible(dev.off())

## 6. Freeze empirical groups and assess functional transfer ----------------
## A clustering tree cannot identify signaling direction. Score orientation is
## fixed from development survival after adjusting drug means; each gene has
## equal weight. No external information assigns a branch or changes the groups.
score_raw <- function(x) sapply(programs,function(p) {
  g<-membership[program==p,gene];present<-intersect(g,rownames(x))
  if(length(present)!=length(g)) return(rep(NA_real_,ncol(x)))
  colMeans(x[g,,drop=FALSE])
})
raw_scores<-score_raw(t(auc[,chosen,drop=FALSE]))
fit_rows<-which(line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
train_meta<-line_drug[fit_rows]
base<-lm(prism_survival_auc~compound,data=train_meta)
orient<-sapply(seq_along(programs),function(j)sign(cor(raw_scores[fit_rows,j],residuals(base),method='spearman')))
check(all(is.finite(orient)&orient!=0),'Unestimable program orientation')
names(orient)<-programs
membership[,orientation:=orient[program]];membership[,weight:=orientation/.N,by=program]
write_table(membership,'frozen_program_genes.tsv')
score<-function(x)sweep(score_raw(x),2,orient,'*')
scores<-sweep(raw_scores,2,orient,'*')
freeze<-list(membership=membership,selected=selected,orientation=orient,tree=tree,cut=k,
  cohort=line_drug,decision='Empirical programs; no BRAF/ERBB/root assignment',code_md5=unname(tools::md5sum(file.path(root,'scripts/26_test_braf_erbb_bypass.R'))))
saveRDS(freeze,file.path(art,'frozen_programs.rds'))
fit_results<-list();held_points<-list();association<-list()
held<-which(line_drug$split=='Holdout' & line_drug$analysis_included & line_drug$eligible)
for(j in seq_along(programs)) {
  p<-programs[j];train<-copy(train_meta);train[,program_score:=scores[fit_rows,j]]
  fit<-lm(prism_survival_auc~compound+program_score,data=train)
  test<-copy(line_drug[held]);test[,program_score:=scores[held,j]]
  test[,`:=`(predicted_baseline=predict(base,test),predicted_program=predict(fit,test),program=p)]
  test[,residual_survival:=prism_survival_auc-predicted_baseline]
  fit_results[[p]]<-data.table(program=p,n_train=nrow(train),n_holdout_lines=uniqueN(test$depmap_id),n_test_observations=nrow(test),
    baseline_rmse=sqrt(mean((test$prism_survival_auc-test$predicted_baseline)^2)),
    program_rmse=sqrt(mean((test$prism_survival_auc-test$predicted_program)^2)),
    pooled_residual_spearman=cor(test$program_score,test$residual_survival,method='spearman'))
  association[[p]]<-test[,.(n_lines=.N,rho=cor(program_score,prism_survival_auc,method='spearman')),by=.(program,compound)]
  held_points[[p]]<-test
}
held_points<-rbindlist(held_points);fit_results<-rbindlist(fit_results)
write_table(fit_results,'heldout_model_comparison.tsv');write_table(rbindlist(association),'heldout_within_drug_correlations.tsv')
write_table(held_points,'heldout_predictions.tsv')
p<-ggplot(held_points,aes(program_score,residual_survival,color=compound))+geom_point(size=2)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+facet_wrap(~program,scales='free')+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+labs(x='Frozen program score',y='Survival minus development drug mean',title='Held-out PRISM association | empirical groups',subtitle='7 independent lines; drug observations from a line are dependent')
plot_save(p,'26i_heldout_prism_association',8,5)
# Context profiles use V600 only for the vemurafenib displacement annotation.
profile<-sapply(drug_map$compound,function(d) {
  ii<-which(line_drug$compound==d & line_drug$split=='Development' & line_drug$analysis_included &
    (d!='vemurafenib'|line_drug$braf_v600_confirmed))
  colMeans(scores[ii,,drop=FALSE])
})
lig_scores<-score(cyto_x)
profile<-cbind(profile,sapply(c('EGF','NRG1'),function(l)apply(lig_scores[study_meta$treatment==l,,drop=FALSE],2,median,na.rm=TRUE)))
rownames(profile)<-programs
write_table(data.table(program=programs,profile),'program_perturbation_profiles.tsv')
ph<-hclust(dist(profile),method='average')
cairo_pdf(file.path(out,'26j_program_hierarchy.pdf'),width=6,height=4.5)
par(mar=c(3,4,3,1));plot(as.dendrogram(ph),main='Program similarity | biological direction unresolved',xlab='',sub='',ylab='Perturbation-profile distance',cex=1.5);invisible(dev.off())
pl<-melt(data.table(program=programs,profile),id.vars='program',variable.name='perturbation',value.name='score')
p<-ggplot(pl,aes(perturbation,score,color=program,group=program))+geom_hline(yintercept=0,color='grey70')+geom_line()+geom_point(size=2)+
  scale_color_manual(values=colors)+theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=25,hjust=1))+
  labs(x=NULL,y='Mean signed displacement',title='Frozen empirical programs | perturbation directions',subtitle='Ligand and Tahoe effect units differ; compare directions within source')
plot_save(p,'26j_program_directions',8,4.5)
# Full freeze precedes loading either external numerical matrix.
write_table(data.table(order=1:4,decision=c('Freeze genes, weights, order and orientation before external expression',
  'Score every empirical program; label observed-gene partial scores explicitly; never assign B/E/R from validation outcomes',
  'No BypassIndex or Y success claim without a resolved receptor and root pair',
  'External displacement is reported even when discovery structure is weak')),'external_prediction_freeze.tsv')
message('Discovery frozen: ',length(chosen),' genes; ',k,' empirical groups. Proceeding to external molecular diagnostics.')

saveRDS(mget(c('programs','membership','orient','colors','leaf_order','chosen','k','curve','coherence','fit_results','selected','line_drug','held','fit_rows','excluded_ids','s684','s127','factorial_lines','external_lines','root','art','out','work','score','score_raw','heat','write_table','check','plot_save','row_blocks'),envir=.GlobalEnv),file.path(art,'discovery_fit.rds'))

## 7. External gene identity and platform coverage ---------------------------
## Gene identifiers are repaired using approved/unique historical HGNC symbols.
## External values cannot change membership, orientation or gene weights.
hgnc_path<-file.path(work,'outputs/dataset_inventory/cache/cellfie/hgnc_complete_set.txt')
hgnc<-fread(hgnc_path);hgnc<-hgnc[status=='Approved']
lookup<-rbindlist(lapply(seq_len(nrow(hgnc)),function(i) {
  a<-unique(c(hgnc$symbol[i],unlist(strsplit(paste(hgnc$prev_symbol[i],hgnc$alias_symbol[i],sep='|'),'\\|'))))
  data.table(input=a[nzchar(a)&!is.na(a)&a!='NA'],symbol=hgnc$symbol[i])
}))
lookup<-lookup[,if(uniqueN(symbol)==1L).(symbol=symbol[1]) else .(symbol=NA_character_),by=input]
canonical<-function(v){ans<-lookup$symbol[match(v,lookup$input)];exact<-v %chin% hgnc$symbol;ans[exact]<-v[exact];ans}
frozen_symbols<-canonical(leaf_order)
check(!anyNA(frozen_symbols)&&!anyDuplicated(frozen_symbols),'Frozen gene aliases are missing or ambiguous')
write_table(data.table(gene=leaf_order,official_symbol=frozen_symbols,hgnc_name=hgnc$name[match(frozen_symbols,hgnc$symbol)],
  locus_type=hgnc$locus_type[match(frozen_symbols,hgnc$symbol)]),'selected_gene_identities.tsv')
# Read the platform annotation embedded in the downloaded SOFT file.
platform_path<-file.path(art,'GSE68453_platform_annotation.tsv')
if(!file.exists(platform_path)) {
  con<-gzfile(file.path(art,'GSE68453_family.soft.gz'),'rt');buffer<-character();reading<-FALSE;done<-FALSE
  repeat {
    lines<-readLines(con,n=10000,warn=FALSE);if(!length(lines))break
    for(l in lines) {
      if(startsWith(l,'!platform_table_end')) {done<-TRUE;break}
      if(reading)buffer<-c(buffer,l)
      if(startsWith(l,'!platform_table_begin'))reading<-TRUE
    }
    if(done)break
  };close(con);check(done,'Platform annotation not found in SOFT');writeLines(buffer,platform_path)
}
platform<-fread(platform_path,quote='');platform[,official_symbol:=canonical(Symbol)]
matfile<-file.path(art,'GSE68453_series_matrix.txt.gz')
m684<-suppressWarnings(fread(matfile,skip='"ID_REF"',header=TRUE,check.names=FALSE,select=c('ID_REF',s684$gsm)))
m684<-m684[ID_REF %chin% platform$ID]
x684<-as.matrix(m684[,s684$gsm,with=FALSE]);storage.mode(x684)<-'double'
check(all(is.finite(x684)),'Non-finite GSE68453 measurements')
# GEO supplies normalized log2 array expression. Assert the scale, do not log twice.
check(quantile(x684,.99)<30 && quantile(x684,.5)>0,'Unexpected GSE68453 expression scale')
probe<-data.table(probe=m684$ID_REF,official_symbol=platform$official_symbol[match(m684$ID_REF,platform$ID)],
  iqr=matrixStats::rowIQRs(x684),array_row=seq_len(nrow(x684)))
probe<-probe[!is.na(official_symbol)];setorder(probe,official_symbol,-iqr,probe)
probe[,chosen_probe:=seq_len(.N)==1L,by=official_symbol]
write_table(probe,'gse68453_probe_mapping.tsv.gz')
keep<-probe[chosen_probe==TRUE];x684<-x684[keep$array_row,,drop=FALSE];rownames(x684)<-keep$official_symbol
# Individual RPM columns are selected by sample ID and independently checked
# against GEO drug doses, exposure times and biological day. Summary columns
# are explicitly excluded; they are not additional replicates.
rpm_path<-file.path(art,'gse127988/GSE127988_RPM_Gerosa_et_al.csv.gz')
rh<-fread(rpm_path,header=FALSE,nrows=6,check.names=FALSE)
header<-as.matrix(rh[,-c(1,2),with=FALSE])
sample_ids<-as.integer(sub('.*\\(S([0-9]+)\\).*','\\1',s127$title))
ci<-vapply(sample_ids,function(id){z<-which(header[4,]=='Individual'&header[5,]==as.character(id));check(length(z)==1L,'Ambiguous RPM sample ID');z},integer(1))
check(all(as.numeric(header[1,ci])==s127$vem_uM)&all(as.numeric(header[2,ci])==s127$cobi_uM)&
  all(as.numeric(header[3,ci])==s127$egf_h)&all(as.numeric(header[6,ci])==s127$biological_set),'RPM/GEO design mismatch')
rpm<-fread(rpm_path,skip=7,header=FALSE,select=c(1L,2L,ci+2L),check.names=FALSE)
numeric_rpm<-as.matrix(rpm[,-c(1,2),with=FALSE]);storage.mode(numeric_rpm)<-'double'
check(all(is.finite(numeric_rpm)&numeric_rpm>=0),'Invalid RPM measurements')
rna_symbol<-hgnc$symbol[match(sub('\\..*','',rpm[[1]]),hgnc$ensembl_gene_id)]
fallback<-canonical(rpm[[2]]);rna_symbol[is.na(rna_symbol)]<-fallback[is.na(rna_symbol)]
rna_map<-data.table(ensembl=rpm[[1]],reported_symbol=rpm[[2]],official_symbol=rna_symbol)
write_table(rna_map,'gse127988_gene_mapping.tsv.gz')
# Several Ensembl rows mapping to one gene are aggregated on the RPM scale.
rna<-data.table(official_symbol=rna_symbol,numeric_rpm)
rna<-rna[!is.na(official_symbol),lapply(.SD,sum),by=official_symbol]
x127<-log2(1+as.matrix(rna[,-1,with=FALSE]));rownames(x127)<-rna$official_symbol;colnames(x127)<-s127$gsm
coverage<-rbindlist(lapply(c('GSE68453','GSE127988'),function(ds) {
  xx<-if(ds=='GSE68453')x684 else x127
  data.table(dataset=ds,gene=leaf_order,official_symbol=frozen_symbols,program=membership$program,
    measured=frozen_symbols %chin% rownames(xx),reason=ifelse(frozen_symbols %chin% rownames(xx),'Mapped','No uniquely mapped platform/RNA measurement'))
}))
write_table(coverage,'external_gene_coverage.tsv')
write_table(data.table(gene=leaf_order,official_symbol=frozen_symbols,ensembl=hgnc$ensembl_gene_id[match(frozen_symbols,hgnc$symbol)],
  raw_symbol_present=frozen_symbols %chin% rpm[[2]],raw_ensembl_present=hgnc$ensembl_gene_id[match(frozen_symbols,hgnc$symbol)] %chin% rpm[[1]]),
  'gse127988_raw_identifier_audit.tsv')
extract_genes<-function(x) {
  m<-matrix(NA_real_,length(leaf_order),ncol(x),dimnames=list(leaf_order,colnames(x)))
  ok<-frozen_symbols %chin% rownames(x);m[ok,]<-x[frozen_symbols[ok],,drop=FALSE];m
}
x684<-extract_genes(x684);x127<-extract_genes(x127)
# Report complete and observed-gene scores separately; no zeros or replacement genes.
external_score<-function(x) {
  # Missingness-only adaptation: keep the frozen membership and orientation,
  # but expose an observed-gene mean separately from the complete score. A
  # missing gene is never zero. Counts/coverage identify every partial score.
  sapply(programs,function(p) {
    g<-membership[program==p,gene]
    present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)]
    if(!length(present))return(rep(NA_real_,ncol(x)))
    colMeans(x[present,,drop=FALSE])*orient[p]
  })
}
write_table(data.table(dataset=coverage$dataset,program=coverage$program,measured=coverage$measured)[,
  .(n_genes=.N,n_measured=sum(measured),complete=all(measured)),by=.(dataset,program)],'external_program_coverage.tsv')

## 8. GSE68453 factorial displacement and independent efficacy association --
# Controls are a within-line reference. Replicate labels are not assumed to be
# matched biological pairs across treatment arms.
d684<-x684
for(l in unique(s684$cell_line)) {
  jj<-which(s684$cell_line==l);ref<-rowMeans(x684[,s684$cell_line==l&s684$arm=='control',drop=FALSE])
  d684[,jj]<-x684[,jj,drop=FALSE]-ref
}
s684scores<-external_score(d684)
long684<-melt(data.table(s684[,.(gsm,cell_line,arm,replicate)],s684scores),id.vars=c('gsm','cell_line','arm','replicate'),variable.name='program',value.name='score')
long684<-merge(long684,coverage[dataset=='GSE68453',.(n_frozen=.N,n_measured=sum(measured)),by=program],by='program',sort=FALSE)
long684[,complete_score:=ifelse(n_frozen==n_measured,score,NA_real_)]
write_table(long684,'gse68453_program_scores.tsv')
means684<-long684[,.(score=mean(score)),by=.(cell_line,arm,program)]
factorial<-dcast(means684[cell_line %chin% factorial_lines],cell_line+program~arm,value.var='score')
factorial[,`:=`(plx_minus_control=plx-control,lap_minus_control=lapatinib-control,
  combo_minus_control=plx_lap-control,adding_lapatinib=plx_lap-plx,interaction=plx_lap-plx-lapatinib+control)]
write_table(factorial,'gse68453_factorial_contrasts.tsv')
# Frozen genes keep their permanent order in both external representations.
external_heat<-function(x,metadata,prefix,semantic_order) {
  finite_genes<-rowSums(is.finite(x))==ncol(x)
  check(sum(finite_genes)>1,'Too few measured frozen genes for an external heatmap')
  c<-cor(x[finite_genes,,drop=FALSE],method='spearman');check(all(is.finite(c)),'Non-estimable external sample correlations')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average');ord<-tr$order
  blocks<-cutree(tr,min(4,ncol(x)-1));bl<-factor(blocks[ord],levels=unique(blocks[ord]))
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),cluster=blocks[ord]),paste0(prefix,'_sample_order.tsv'))
  heat(c[ord,ord],paste0(prefix,'_sample_correlations'),paste(prefix,'| sample correlations'),bl,bl,limit=1,width=12,height=11,show_cols=FALSE)
  heat(x[leaf_order,ord,drop=FALSE],paste0(prefix,'_gene_samples_clustered'),paste(prefix,'| log-expression displacement | frozen genes'),row_blocks,bl,width=13)
  heat(x[leaf_order,semantic_order,drop=FALSE],paste0(prefix,'_gene_samples_semantic'),paste(prefix,'| experimental order | frozen genes'),row_blocks,width=13)
}
external_heat(d684,s684,'26k_gse68453',order(match(s684$cell_line,external_lines),match(s684$arm,c('control','plx','lapatinib','plx_lap')),s684$replicate))
p<-ggplot(means684[cell_line %chin% factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap')),score,color=program,group=program))+
  geom_hline(yintercept=0,color='grey70')+geom_line()+geom_point(size=2)+facet_wrap(~cell_line,ncol=3)+
  scale_color_manual(values=colors)+theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=30,hjust=1))+
  labs(x=NULL,y='Frozen program displacement',title='GSE68453 | five complete factorials',subtitle='Same genes and weights; empirical program labels')
plot_save(p,'26l_factorial_programs',10,7)
# Write expression-only predictors BEFORE attaching digitized efficacy.
predictors684<-means684[arm=='plx',.(cell_line,program,plx_displacement=score)]
write_table(predictors684,'gse68453_frozen_expression_predictors.tsv')
bliss<-data.table(cell_line=external_lines,bliss=c(5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1))
benefit<-merge(predictors684,bliss,by='cell_line')
benefit_stats<-benefit[,{
  ok<-is.finite(plx_displacement)&is.finite(bliss)
  if(sum(ok)>=3) {tt<-cor.test(plx_displacement[ok],bliss[ok],method='spearman',exact=FALSE);.(n=sum(ok),rho=unname(tt$estimate),p=tt$p.value)}
  else .(n=sum(ok),rho=NA_real_,p=NA_real_)
},by=program]
benefit_stats[,bh_q:=p.adjust(p,'BH')]
write_table(benefit,'gse68453_expression_bliss_join.tsv');write_table(benefit_stats,'gse68453_exploratory_program_bliss.tsv')
p<-ggplot(benefit,aes(plx_displacement,bliss))+geom_point(color='#287D99',size=2)+geom_text(aes(label=cell_line),size=3,vjust=-.6,check_overlap=TRUE)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='#B64039')+facet_wrap(~program,scales='free_x')+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Frozen PLX program displacement',y='Digitized Bliss benefit',title='Exploratory external efficacy association',subtitle='No receptor/root pair resolved: these are not BypassIndex tests')
plot_save(p,'26m_external_efficacy',9,5)

## 9. GSE127988 ligand displacement, respecting biological day --------------
d127<-x127
for(arm in unique(s127$arm))for(day in unique(s127$biological_set)) {
  jj<-which(s127$arm==arm&s127$biological_set==day)
  ref<-rowMeans(x127[,s127$arm==arm&s127$biological_set==day&s127$egf_h==0,drop=FALSE])
  d127[,jj]<-x127[,jj,drop=FALSE]-ref
}
s127scores<-external_score(d127)
long127<-melt(data.table(s127[,.(gsm,arm,egf_h,biological_set)],s127scores),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='program',value.name='score')
long127<-merge(long127,coverage[dataset=='GSE127988',.(n_frozen=.N,n_measured=sum(measured)),by=program],by='program',sort=FALSE)
long127[,complete_score:=ifelse(n_frozen==n_measured,score,NA_real_)]
write_table(long127,'gse127988_program_scores.tsv')
day127<-long127[,.(score=mean(score)),by=.(arm,egf_h,biological_set,program)]
write_table(day127,'gse127988_biological_day_means.tsv')
auc127<-day127[order(egf_h),.(auc=sum(diff(egf_h)*(head(score,-1)+tail(score,-1))/2)/8),by=.(arm,biological_set,program)]
write_table(auc127,'gse127988_auc_by_day.tsv')
attenuation<-dcast(auc127,program+biological_set~arm,value.var='auc')
setnames(attenuation,c('Vem','Vem + Cobi'),c('vem_auc','vem_cobi_auc'))
attenuation[,`:=`(difference=vem_cobi_auc-vem_auc,attenuation_fraction=ifelse(vem_auc>0,1-vem_cobi_auc/vem_auc,NA_real_))]
write_table(attenuation,'gse127988_auc_attenuation.tsv')
external_heat(d127,s127,'26n_gse127988',order(s127$cobi_uM,s127$egf_h,s127$biological_set,s127$gsm))
p<-ggplot(day127,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+
  geom_hline(yintercept=0,color='grey70')+geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+facet_wrap(~program)+
  theme_classic(base_size=12)+labs(x='EGF exposure (h)',y='Frozen program displacement from own time zero',color='Drug background',linetype='Biological day',
  title='GSE127988 | frozen programs under EGF',subtitle='Technical wells averaged within day; two biological days shown separately')
plot_save(p,'26o_egf_time_courses',9,5)
p<-ggplot(auc127,aes(arm,auc,color=program,group=interaction(program,biological_set)))+geom_hline(yintercept=0,color='grey70')+
  geom_line(aes(linetype=factor(biological_set)))+geom_point(size=2)+facet_wrap(~program)+scale_color_manual(values=colors)+
  theme_classic(base_size=12)+labs(x=NULL,y='Normalized 0–8 h AUC',linetype='Biological day',title='EGF program activation with and without MEK inhibition')
plot_save(p,'26p_egf_auc',8,4.5)

## 10. Reviewable conclusions, decisions and reproducibility ----------------
# No validation result feeds back into gene selection or program membership.
check(identical(as.data.frame(readRDS(file.path(art,'frozen_programs.rds'))$membership),as.data.frame(membership)),'Frozen programs changed during validation')
check(!any(line_drug$depmap_id[c(fit_rows,held)] %chin% excluded_ids),'Excluded line entered analysis')
check(length(intersect(line_drug$depmap_id[fit_rows],line_drug$depmap_id[held]))==0L,'A line crossed discovery/holdout')
write_table(data.table(check=c('Excluded outliers absent from discovery and validation','Discovery/holdout line separation',
  'Three doses retained per included line and drug','At most 40 genes admitted per drug','No historical feature list controls selection',
  'CytoSig study balance enters discovery','Frozen membership unchanged after external scoring','External gene coverage explicit',
  'No causal Y labels forced','All source and external heatmaps use frozen gene order'),passed=TRUE),'verification.tsv')
decisions<-data.table(step=c('Organization','Cell lines','Vemurafenib selection','Gene measurability','Association influence','Doses','CytoSig','Panel size','Clustering','Program interpretation','External validation'),
  observed_failure=c('Unnumbered pilot entry point','Known outliers restored into input audit','Only three V600 lines; 5394 maximum-rho ties',
  '53 of 120 prior selected genes mostly below 0.1 CP10K','Rank associations can depend on a line or lineage','Weak/disagreeing doses require inspection',
  '15 of 16 NRG1 profiles from one study','60 genes and two groups were too coarse for the requested inspection','Four coordinates manufacture coarse rank structure','Similarity does not identify signaling direction','Stopping at Y gate hid portability and displacement'),
  adaptation=c('Numbered stage-26 script; matching result/artifact directories','Exclude NCI-H2122 and NCI-H596; report influence comparison',
  'Use 26 retained development solid-tumor lines; V600 used only for context annotation',
  'Require >=0.1 CP10K in at least half of drug samples; no global 3000-gene panel',
  'Rank by worst leave-one-line association; require direction across tissue deletions',
  'Keep 0.05, 0.5 and 5 uM; fixed AUC; report each line dose concordance',
  'One median per study x ligand; retain all 25 study-ligand observations',
  'Expand the stable measured panel to at most 40 genes per drug / 120 total',
  paste0('Average actual Tahoe/CytoSig gene correlations; requested finer resolution k=',k,'; silhouette and sizes remain diagnostic'),
  'Keep empirical labels; no Y, upstream order or BypassIndex claim',
  'Freeze groups first; audit identifiers; distinguish complete from partial scores; report both external datasets descriptively'))
write_table(decisions,'adaptation_decisions.tsv')
write_table(data.table(stage=c('Input and QC','Discovery','Held-out PRISM','GSE68453','GSE127988','Mechanistic Y/BypassIndex'),
  status=c(rep('Executed',5),'Not identified; no mechanistic success claim')),'execution_status.tsv')
write_table(data.table(figure=sort(list.files(out,pattern='^26.*\\.pdf$'))),'FIGURE_INDEX.tsv')
writeLines(c('# Stage 26: BRAF–ERBB experiment','',
  '[Code](../../scripts/26_test_braf_erbb_bypass.R) | [Results viewer](index.html) | [Decision log](adaptation_decisions.tsv)','',
  'Run from /work: `Rscript tahoe_prism_v2/scripts/26_test_braf_erbb_bypass.R`.','',
  '## Current result','',paste0(length(chosen),' measurable function-associated genes form ',k,' empirical groups. Mean silhouette: ',round(curve$mean_silhouette[match(k,curve$k)],3),'.'),
  'Source-specific coherence, held-out model performance, all external gene identities, factorial displacement, exploratory efficacy associations and EGF time courses are reported. These groups do not establish the requested BRAF/ERBB/shared-root hierarchy. A Y-shaped causal diagram would overstate the evidence.','',
  '## Design revisions and their limits','',
  'The prior strict-brief run is preserved in artifacts/26_braf_erbb_bypass/archive/fixed_brief. It does not determine the new gene set. Its selected genes are read only to document the prior measurement failure.',
  'All 35 original lines are audited. NCI-H2122 and NCI-H596 are excluded as user-confirmed outliers, independently supported by sparse measurements and abnormal scaling. This leaves 26 development and seven holdout lines. PRISM missing curves remain explicitly unavailable; no response is imputed.',
  'The three confirmed V600 lines cannot rank genes reliably. Discovery now estimates cross-solid-tumor associations, with leave-one-line and leave-one-tissue checks. This is not BRAF-mutant-specific discovery. The V600 subset annotates drug displacement only.',
  'Genes must be measured in both sources, have finite values in at least half of the 25 CytoSig study-ligand summaries, and reach 0.1 counts per 10,000 in at least half of the Tahoe drug samples. This is an explicit abundance heuristic to prevent ranks of near-zero measurements from dominating; it is not an assay detection limit.',
  'Selection requires approximate nominal Spearman p<0.05, at least two dose directions agreeing, worst leave-one-line signed rho>=0.30, and no association reversal after deleting a tissue. Rank by worst leave-one-line rho, then worst leave-one-tissue rho; cap 40 per drug. The small nominal p values are selection statistics, not corrected claims of gene discovery.',
  'All three doses remain in the 0.25/0.50/0.25 AUC. Dose discordance is reported; weak biological signal alone does not justify dropping a dose.',
  'Clustering uses the average of two gene correlation matrices: 78 Tahoe development line-drug AUC profiles and 25 CytoSig study-ligand medians. Each source gets half the weight. Every gene pair uses available CytoSig studies, with pair counts exported. EGF has 23 studies and NRG1 two; NRG1 generalizability remains limited.',
  'The requested finer resolution uses six descriptive groups, with silhouette diagnostics for k=2:10. Group sizes and negative silhouettes remain visible, including small or weak groups. No cut is chosen because it resembles a Y. Per-source correlations expose weak support that the old four-coordinate fingerprint concealed.',
  'Program scores are equal-weight means of the frozen genes, oriented using development survival residuals after drug means. All program models are reported; no best holdout model is selected. Pooled holdout observations repeat seven cell lines across drugs and are not independent sample units.',
  'External identifiers use the cached HGNC approved symbol table and unambiguous previous names/aliases; Illumina probes use maximum IQR across the selected 102 samples, before efficacy is joined. RPM Ensembl IDs are preferred over names, duplicates sum on the RPM scale, then log2(1+RPM). Missing frozen genes are reported. Complete scores abstain for incomplete groups; a separately labeled observed-gene mean renormalizes the remaining frozen weights by available gene count. Partial scores are molecular diagnostics, not full signature validation.',
  'External expression is never used to revise programs. GSE68453 reports five complete factorials and twelve-line exploratory program–Bliss associations, with BH adjustment across programs. There is no BypassIndex without a resolved receptor/root pair. Digitized Bliss values are approximate published-figure readings, not raw plate measurements.',
  'GSE127988 references each drug background to its own time zero within biological day. Technical wells are averaged; both biological days are shown. AUC attenuation is undefined if the vemurafenib-background AUC is nonpositive. There is no selective root-blockade claim without a root identity.',
  'Correlation heatmaps mask self-correlations and share a data-derived symmetric 95th-percentile off-diagonal scale; clipping fractions are recorded. External gene heatmaps share ±2.5 log-expression-displacement limits. No gene order changes between panels.',
  'The earlier reuse of these public data makes external analyses retrospective molecular validation, not a prospective test. The within-run freeze prevents feedback from the current external scoring into discovery.','',
  '## Provenance','',
  '- [GSE68453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68453)',
  '- [GSE127988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988)',
  '- [Prior run archive](../../artifacts/26_braf_erbb_bypass/archive/fixed_brief/results/index.html)'),file.path(out,'README.md'))
html_escape<-function(x){x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
figures<-sort(list.files(out,pattern='^26.*\\.pdf$'));tables<-sort(list.files(out,pattern='\\.tsv(\\.gz)?$'))
cards<-vapply(figures,function(f) {
  png<-sub('\\.pdf$','.png',f)
  paste0('<article><h3>',html_escape(sub('\\.pdf$','',f)),'</h3>',if(file.exists(file.path(out,png)))paste0('<a href="',png,'"><img loading="lazy" src="',png,'"></a>') else '',
    '<a href="',f,'">Open PDF</a></article>')
},character(1))
code_path<-file.path(root,'scripts/26_test_braf_erbb_bypass.R')
html_table<-function(x){paste0('<table><tr>',paste0('<th>',names(x),'</th>',collapse=''),'</tr>',paste(apply(as.data.frame(x),1,function(r)paste0('<tr>',paste0('<td>',html_escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table>')}
writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 26 — BRAF–ERBB experiment</title>',
  '<style>body{font:16px/1.5 system-ui;background:#f4f6f8;color:#203040;margin:0}main{max-width:1300px;margin:auto;padding:28px}a{color:#075b8e}.status{background:#fff3e4;border-left:5px solid #ba7030;padding:16px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px}article{background:white;padding:16px;border:1px solid #d3dce3;overflow:hidden}h3{font-size:15px;overflow-wrap:anywhere}img{max-width:100%;max-height:650px;display:block;margin:auto}table{border-collapse:collapse;width:100%;background:white;font-size:14px}th,td{padding:8px;border:1px solid #ddd;text-align:left}pre{background:#18252f;color:#e3edf5;padding:16px;overflow:auto;max-height:700px}.tables{columns:2}summary{cursor:pointer;font-weight:bold}</style><main>',
  '<h1>Stage 26 · BRAF–ERBB experiment</h1><p><a href="../../scripts/26_test_braf_erbb_bypass.R">Open the numbered R script</a> · <a href="README.md">Methods and limitations</a> · <a href="adaptation_decisions.tsv">Decisions and evidence</a></p>',
  paste0('<div class="status"><strong>Analysis completed; mechanistic Y not established.</strong><p>',length(chosen),' genes, ',k,' empirical groups; mean silhouette ',round(curve$mean_silhouette[match(k,curve$k)],3),'. Both outlier lines are excluded. Frozen groups are scored in held-out PRISM and both external expression datasets. Weak structure is reported, not relabeled as a receptor/root hierarchy.</p></div>'),
  '<h2>What changed, and why</h2>',html_table(decisions),
  '<h2>Actual source coherence</h2>',html_table(coherence),
  '<h2>Held-out functional transfer</h2>',html_table(fit_results),
  '<h2>External gene coverage</h2>',html_table(coverage[,.(measured=sum(measured),total=.N),by=dataset]),
  '<h2>Figures</h2><p>Click any image for full-size viewing; PDFs retain scalable labels.</p><div class="grid">',cards,'</div>',
  '<h2>Tables</h2><ul class="tables">',paste0('<li><a href="',tables,'">',tables,'</a></li>'),'</ul>',
  '<h2>Code</h2><details><summary>Expand the complete stage-26 R script</summary><pre><code>',paste(html_escape(readLines(code_path,warn=FALSE)),collapse='\n'),'</code></pre></details></main></html>'),file.path(out,'index.html'))
message('Stage 26 completed. Review results/26_braf_erbb_bypass/index.html; no mechanistic Y claim.')
