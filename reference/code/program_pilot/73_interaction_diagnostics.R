#!/usr/bin/env Rscript
# Diagnose existing fits without fitting another predictor or selecting test cases.
suppressPackageStartupMessages({library(data.table); library(qs); library(Matrix); library(ggplot2)})
root <- '/work/results/program_pilot/73_bt20_almanac_program_interactions'
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0('diagnostic_',name,'.tsv.gz')),sep='\t')
safe_cor <- function(x,y) if(sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y)
means <- fread(file.path(root,'data/figure_condition_predictions.tsv.gz'))
# Establish whether failure is learned amplitude, partner transfer, or both.
means[,group:=fifelse(actual_count<single_product-10,'Below independence by >10pp',
  fifelse(actual_count>single_product+10,'Above independence by >10pp','Within 10pp of independence'))]
strong <- means[method!='Independent singles',.(n=.N,measured_log2=median(observed),predicted_log2=median(predicted),
  count_bias_pp=median(predicted_count-actual_count),interaction_r=safe_cor(observed,predicted)),by=.(method,split,group)]
save_table(strong,'departure_amplitude')
partner <- means[method!='Independent singles',.(n=.N,interaction_r=safe_cor(observed,predicted),
  measured_log2=median(observed),predicted_log2=median(predicted)),by=.(method,split,partner)]
save_table(partner,'partner_train_held')

# Same underlying Tahoe expression, full-bank versus only genes measured in BT20.
# This isolates missing-gene effects; it cannot establish cross-platform equivalence.
bank <- qread('/work/results/program_pilot/61_individual_responses_combination_transfer/data/representation_scores.qs')
molecular <- qread('/work/results/program_pilot/61_individual_responses_combination_transfer/data/matched_molecular.qs')
W <- qread('/work/results/program_pilot/61_individual_responses_combination_transfer/models/original40_weights.qs')
e <- qread('/work/results/program_pilot/72_almanac_full_program_interactions/data/prepared_inputs.qs')
measured <- intersect(rownames(W),rownames(molecular$displacement))
stopifnot(length(measured)>0L,ncol(W)==40L)
full <- partial <- list()
for(cell in unique(e$metadata$recipient)) {
  path <- file.path('/work/results/program_pilot/21_almanac_two_partner_ward_programs',
    paste0('cell_',gsub('[^A-Za-z0-9]+','_',cell)),'models/partner_queries.qs')
  q <- qread(path)
  full[[cell]] <- as.matrix(q$displacement[,rownames(W),drop=FALSE]%*%W)
  partial[[cell]] <- as.matrix(q$displacement[,measured,drop=FALSE]%*%W[measured,,drop=FALSE])
}
F <- do.call(rbind,full)
P <- do.call(rbind,partial)
stopifnot(identical(dim(F),dim(P)),ncol(P)==40L)
coverage <- data.table(program=colnames(F),full_vs_partial_r=vapply(seq_len(ncol(F)),
  function(j) safe_cor(F[,j],P[,j]),numeric(1)))
coverage <- merge(coverage,bank$coverage[,.(program,absolute_weight_fraction,measured_genes,total_genes)],by='program')
save_table(coverage,'program_measurement_consistency')

# Replay BT20-only and joint fitted coefficients in each source's original scaled
# coordinates; compare training compatibility without re-estimating anything.
cross <- function(a,b) {
  ij <- which(upper.tri(matrix(0,40,40),diag=TRUE),arr.ind=TRUE)
  x <- vapply(seq_len(nrow(ij)),function(k) (a[,ij[k,1]]*b[,ij[k,2]]+a[,ij[k,2]]*b[,ij[k,1]])/2,numeric(nrow(a)))
  colnames(x) <- paste(colnames(a)[ij[,1]],colnames(a)[ij[,2]],sep='__')
  x
}
er <- sqrt(colMeans(e$query_scores^2))
EX <- cross(sweep(e$A,2,er,'/'),sweep(e$B,2,er,'/'))
obs <- fread('/work/results/program_pilot/61_individual_responses_combination_transfer/data/matched_molecular_links.tsv.gz')[time_h==3 & molecular_supported==TRUE]
setorder(obs,anchor,partner,nominal_A_uM,nominal_B_uM,biological_replicate)
r <- bank$original40[bank$metadata[time_h==3,profile_id],,drop=FALSE]
r <- sweep(r,2,sqrt(colMeans(r^2)),'/')
BX <- cross(r[obs$profile_A,,drop=FALSE],r[obs$profile_B,,drop=FALSE])
y <- log2(obs$ratio_AB/(obs$ratio_A*obs$ratio_B))
compatibility <- objective_mass <- list()
for(held in sort(unique(obs$partner))) {
  joint <- qread(file.path(root,'models',paste0(held,'_Full.qs')))
  old <- qread(file.path('/work/results/program_pilot/71_bootstrap_screened_program_interactions/models',paste0(held,'_Full.qs')))
  bt <- which(obs$partner!=held)
  ex <- joint$ALMANAC_training_rows
  for(method in c('BT20-only','Joint')) {
    beta <- if(method=='Joint') joint$raw_coefficients else old$raw_coefficients
    stopifnot(identical(names(beta),colnames(EX)),identical(names(beta),colnames(BX)))
    for(source in c('BT20','ALMANAC')) {
      pred <- as.numeric(if(source=='BT20') BX[bt,,drop=FALSE]%*%beta else EX[ex,,drop=FALSE]%*%beta)
      yy <- if(source=='BT20') y[bt] else e$metadata$observed_g[ex]
      compatibility[[paste(held,method,source)]] <- data.table(held,method,source,n=length(yy),
        r=safe_cor(yy,pred),observed_median=median(yy),predicted_median=median(pred),
        observed_p05=quantile(yy,.05),observed_p95=quantile(yy,.95),
        predicted_p05=quantile(pred,.05),predicted_p95=quantile(pred,.95))
    }
  }
  # Internal fitting diagnostic only: test whether near-zero outcomes dominate
  # the actual BT20 raw+contrast objective. Not a performance metric.
  group <- interaction(obs$anchor[bt],obs$partner[bt],obs$biological_replicate[bt],drop=TRUE)
  w <- 1/as.numeric(table(group)[as.character(group)])
  w <- w/sum(w)
  context <- interaction(obs$anchor[bt],obs$nominal_A_uM[bt],obs$biological_replicate[bt],drop=TRUE)
  pred <- as.numeric(BX[bt,,drop=FALSE]%*%joint$raw_coefficients)
  residual <- y[bt]-pred
  centered <- residual-ave(residual,context,FUN=mean)
  mass <- .5*w*(residual^2+centered^2)
  near <- abs(y[bt])<.25
  objective_mass[[held]] <- data.table(held,near_zero_rows=sum(near),rows=length(bt),
    near_zero_weight=sum(w[near]),near_zero_objective_fraction=sum(mass[near])/sum(mass))
}
save_table(rbindlist(compatibility),'source_transfer')
save_table(rbindlist(objective_mass),'objective_mass')

# All-partner training/held view directly exposes transfer-specific failures.
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
d <- means[method=='Singles + program interactions']
d[,split:=factor(split,levels=c('training','held'))]
labels <- d[,.(label=sprintf('r = %.2f',safe_cor(observed,predicted))),by=.(split,partner)]
p <- ggplot(d,aes(observed,predicted))+geom_point(aes(color=anchor),alpha=.4,size=1)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.4)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.2,size=3)+
  facet_wrap(vars(split,partner),ncol=6,scales='free')+theme(aspect.ratio=1,strip.text=element_text(size=8))+
  labs(title='The same partners when seen in training versus entirely excluded',
    subtitle='Training includes five distinct outer fits per condition; held predictions come from the one partner-excluding fit.',
    x='Observed log2 departure from independence',y='Predicted log2 departure',color='First drug')
ggsave(file.path(root,'figures/diagnostic_each_partner_train_held.pdf'),p,width=14,height=6,device=cairo_pdf)
ggsave(file.path(root,'figures/diagnostic_each_partner_train_held.png'),p,width=14,height=6,dpi=140)
notes <- c('Diagnostics only: no model refit or outcome-based selection.',
  'Full/partial program comparison uses identical Tahoe profiles, changing only measured gene availability.',
  'Source replay uses each original model coefficient vector without recalibration; different sources have separate RMS coordinates.',
  'Objective fractions are internal diagnostics of what dominates fitting, not evidence of biological predictive performance.',
  '',capture.output(print(strong)), '',capture.output(print(partner)),
  '',capture.output(print(coverage)), '',capture.output(print(rbindlist(objective_mass))),
  '',capture.output(print(rbindlist(compatibility)[,.(median_r=median(r),
    median_predicted_p05=median(predicted_p05),median_predicted_p95=median(predicted_p95)),by=.(method,source)])))
writeLines(notes,file.path(root,'DIAGNOSTICS.txt'))
cat(paste(notes,collapse='\n'))
