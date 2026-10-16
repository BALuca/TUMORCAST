#!/usr/bin/env Rscript
# Descriptive comparison only: fixed single-drug signatures beside combination outcomes.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
prior <- '/work/results/program_pilot/135_v128_with_corrected_HGF_inputs/linear'
out <- file.path(root,'figures/single_signatures_beside_combination_benefit')
data_out <- file.path(root,'data/single_signatures_beside_combination_benefit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(data_out,recursive=TRUE,showWarnings=FALSE)
# v135 and v136 have identical molecular inputs, training outcomes and single-drug model.
for (file in c('combination_training_outcomes.tsv','external_combination_outcomes.tsv',
               'mema_nonchallenge_training_outcomes.tsv','external_ligand_programs.tsv',
               'external_predicted_drug_programs.tsv','mema_molecular_inputs.tsv')) {
  stopifnot(identical(fread(file.path(root,'data',file)),fread(file.path(prior,'data',file))))
}
model <- qread(file.path(root,'models/final.qs'))
old <- qread(file.path(prior,'models/final.qs'))
stopifnot(identical(model$template,old$template),identical(model$rna,old$rna),
          identical(model$context$sigma,old$context$sigma))
R <- model$template$correlations
S <- model$template$collapsed
stopifnot(identical(rownames(R),rownames(S)),nrow(R)==40,ncol(S)==2)
programs <- rownames(R)
row_order <- hclust(dist(R),method='ward.D2')$order
col_order <- hclust(dist(t(R)),method='ward.D2')$order
selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
stopifnot(all(colnames(R) %in% selection$drug))
colnames(R) <- paste0(sub(' ditosylate| \\(hydrochloride\\)','',colnames(R)), ' ',
  selection$dose[match(colnames(R),selection$drug)], ' µM')
colnames(S) <- c('Signature 1','Signature 2')
corr <- fread(file.path(prior,'data/combination_departure_signatures/program_departure_correlations.tsv'))
focus <- copy(corr[dataset %in% c('PD-184352','PLX-4720') & grepl('^HGF \\+',context)])
stopifnot(nrow(focus)==160L)
focus[,dataset:='HGF_focus']
corr <- rbind(corr,focus)
# Reuse already-calculated exact Pearson correlations; only plotting order changes.
fwrite(corr,file.path(data_out,'program_combination_benefit_correlations.tsv'),sep='\t')
fwrite(data.table(program=programs[row_order],order=seq_len(40)),
  file.path(data_out,'program_order_from_single_drug_correlations.tsv'),sep='\t')
fwrite(data.table(program=programs,R),file.path(data_out,'single_drug_sensitivity_correlations.tsv'),sep='\t')
fwrite(data.table(program=programs,S),file.path(data_out,'fixed_collapsed_signatures.tsv'),sep='\t')
ids <- c('HGF_focus','HER2','PD-184352','PLX-4720','ALMANAC','Growth_across_ligands','MEMA_relative_only')
titles <- c(HGF_focus='HGF training combinations | across six cells per inhibitor',
 HER2='Lapatinib + ligands | across cells at fixed drug/ligand/dose',
 'PD-184352'='MEK inhibitor PD-184352 + ligands | across six cells per ligand',
 'PLX-4720'='RAF inhibitor PLX-4720 + ligands | across six cells per ligand',
 ALMANAC='Drug + drug | across cells at fixed pair and dose',
 Growth_across_ligands='RAF/MEK + ligands | across 13 ligands within each cell and drug',
 MEMA_relative_only='Lapatinib + ligands | relative response across 26 ligands per cell')
notes <- c(HGF_focus='Both HGF combinations show positive measured growth departure in 6/6 cells. Columns correlate variation in that departure.',
 HER2='Combination endpoint: final-abundance departure at 72 h. EGF has only two cells: grey = undefined.',
 'PD-184352'='Combination endpoint: reported net-growth departure, not final cell abundance. HGF is a training ligand.',
 'PLX-4720'='Combination endpoint: reported net-growth departure, not final cell abundance. HGF is a training ligand.',
 ALMANAC='Combination endpoint: log2[V(A+B)/(V(A)V(B))]. A/B follow the displayed drug order.',
 Growth_across_ligands='Reported net-growth departure. Each column compares ligands within a fixed cell/drug; HGF is included.',
 MEMA_relative_only='Cell-centred log2[drug/vehicle under the same ligand]. No absolute interaction contrast. HGF/NRG1 excluded.')
heat <- function(M,name,title,first=FALSE,width=5*ncol(M),legend=FALSE) {
  draw_heatmap(M,scale_rows=FALSE,value_range=c(-1,0,1),
    value_colors=c('#2166AC','white','#B2182B'),hm_name=name,na_col='#D0D0D0',
    show_heatmap_legend=legend,heatmap_legend_param=list(title='Correlation / mean r'),
    show_row_names=first,row_names_side='left',row_names_gp=grid::gpar(fontsize=8),
    show_column_names=TRUE,column_names_rot=60,column_names_gp=grid::gpar(fontsize=7),
    column_names_max_height=grid::unit(65,'mm'),width=grid::unit(width,'mm'),
    height=grid::unit(140,'mm'),title=title,column_title_gp=grid::gpar(fontsize=9))
}
render <- function(id) {
  dd <- corr[dataset==id]
  parts <- unique(dd$partner)
  columns <- unique(dd$context)
  matrices <- lapply(parts,function(part) {
    wide <- dcast(dd[partner==part],program~context,value.var='r')
    M <- as.matrix(wide[match(programs,program),..columns])
    dimnames(M) <- list(programs,columns)
    M
  })
  joint <- do.call(rbind,matrices)
  valid <- which(colSums(is.finite(joint))==nrow(joint))
  order <- if(length(valid)>1) c(valid[hclust(dist(t(joint[,valid,drop=FALSE])),method='ward.D2')$order],
    setdiff(seq_along(columns),valid)) else seq_along(columns)
  ht <- heat(R[row_order,col_order,drop=FALSE],'single','Single-drug sensitivity\n(signature-discovery bank)',TRUE,legend=TRUE) +
    heat(S[row_order,,drop=FALSE],'signatures','Fixed\nsignatures',width=14)
  for(j in seq_along(parts)) {
    part <- parts[j]
    M <- matrices[[j]][row_order,order,drop=FALSE]
    nn <- dd[partner==part,.(n=first(n)),by=context]
    counts <- nn$n[match(colnames(M),nn$context)]
    labels <- if(id=='ALMANAC') paste('Drug',part) else if(part=='A') 'Drug' else 'Ligand'
    stopifnot(identical(rownames(M),programs[row_order]),!anyNA(counts))
    colnames(M) <- paste0(sub('Lapatinib ditosylate','Lapatinib',colnames(M),fixed=TRUE),' [n=',counts,']')
    panel_title <- if(ncol(M)<=4) paste0(labels,'\nvs extra\nsurvival') else paste0(labels,' programs vs\nadditional survival')
    ht <- ht + heat(M,paste0('combo_',part),panel_title,
      width=max(20,5*ncol(M)))
  }
  ComplexHeatmap::draw(ht,heatmap_legend_side='bottom',column_title=titles[[id]],
    column_title_gp=grid::gpar(fontsize=11),ht_gap=grid::unit(4,'mm'),
    padding=grid::unit(c(30,4,4,4),'mm'))
  footer <- paste0(notes[[id]],'\nLeft: positive = association with greater SINGLE-DRUG SENSITIVITY. Right: positive = association with EXTRA SURVIVAL / protection.',
    '\nIdentical program order from single-drug correlations; no row scaling. Drug displacements predicted; ligand displacements transported.',
    '\nDescriptive training-data associations, not held-out predictions or measured combination RNA. Lapatinib+HGF and Wilson outcomes excluded.')
  footer <- paste(unlist(lapply(strsplit(footer,'\n',fixed=TRUE)[[1]],strwrap,
    width=floor(widths[id]*18))),collapse='\n')
  grid::grid.text(footer,x=.5,y=grid::unit(3,'mm'),just=c('centre','bottom'),gp=grid::gpar(fontsize=6.5))
}
widths <- vapply(ids,function(id) {
  z <- corr[dataset==id]
  max(7.2,(5*ncol(R)+14+uniqueN(z$partner)*max(20,5*uniqueN(z$context))+48)/25.4)
},numeric(1))
book <- file.path(out,'single_signatures_and_combination_benefit.pdf')
cairo_pdf(book,width=max(widths),height=9.6,onefile=TRUE)
for(id in ids) render(id)
invisible(dev.off())
for(id in ids) {
  cairo_pdf(file.path(out,paste0(id,'.pdf')),width=widths[id],height=9.6)
  render(id)
  invisible(dev.off())
  png(file.path(out,paste0(id,'.png')),width=widths[id],height=9.6,units='in',res=140,type='cairo')
  render(id)
  invisible(dev.off())
}
external <- fread(file.path(root,'data/external_combination_outcomes.tsv'))
rescue <- external[,.(n=.N,positive_departure=sum(observed_departure_log2>0),
  median_departure_log2=median(observed_departure_log2)),by=.(study,drug,dose,gene)]
fwrite(rescue,file.path(data_out,'measured_rescue_by_combination.tsv'),sep='\t')
print(rescue[gene=='HGF'])
cat('Saved:',book,'\n')
