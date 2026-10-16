#!/usr/bin/env Rscript
# Shared PLS1 signatures: measured displacement -> sensitivity variation.
# Run with Rscript --vanilla; use existing project libraries without renv startup.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/90_shared_pls_signature_decoder'
for (d in c('data','models','figures')) dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
v84 <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/data'
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge/data'
programs <- sprintf('P%02d',1:40)
max_components <- 6L
selected <- fread(file.path(v87,'selected_doses.tsv'))
points <- merge(fread(file.path(v84,'observed_program_viability_points.tsv.gz')),selected[,.(drug,dose)],by=c('drug','dose'))
z <- dcast(points,depmap_id+drug+dose+observed_log2_abundance~program,value.var='observed_displacement')
z <- merge(z,fread(file.path(v87,'recipient_folds.tsv'))[,.(depmap_id,fold)],by='depmap_id')
setorder(z,drug,depmap_id)
X <- as.matrix(z[,..programs])
y <- z$observed_log2_abundance
stopifnot(nrow(z)==256L,uniqueN(z$depmap_id)==30L,uniqueN(z$drug)==9L,!anyNA(X),!anyDuplicated(z[,.(depmap_id,drug)]))
write_table <- function(x,n) fwrite(x,file.path(root,'data',paste0(n,'.tsv')),sep='\t')
pearson <- function(a,b) if(length(a)>=4 && sd(a)>1e-12 && sd(b)>1e-12) cor(a,b) else NA_real_
# PLS packages are unavailable. Standard PLS1 NIPALS deflation implemented here:
# w proportional to X'y; t=Xw; deflate X by tp' and y by tq.
# Equal total weight per drug; scale programs using training residual SD only.
fit_pls <- function(rows) {
  groups <- match(z$drug[rows],unique(z$drug[rows]))
  counts <- tabulate(groups)
  xx <- X[rows,,drop=FALSE]
  yy <- y[rows]
  for(g in unique(groups)) {
    ii <- which(groups==g)
    xx[ii,] <- sweep(xx[ii,,drop=FALSE],2,colMeans(xx[ii,,drop=FALSE]))
    yy[ii] <- yy[ii]-mean(yy[ii])
  }
  weights <- 1/counts[groups]
  weights <- weights/mean(weights)
  scale <- sqrt(colSums(xx^2*weights)/sum(weights))
  stopifnot(all(scale>1e-12))
  E <- sweep(xx,2,scale,'/')*sqrt(weights)
  f <- yy*sqrt(weights)
  W <- P <- matrix(0,40,max_components)
  Q <- numeric(max_components)
  B <- matrix(0,40,max_components)
  for(k in seq_len(max_components)) {
    w <- as.numeric(crossprod(E,f))
    stopifnot(sum(w^2)>1e-20)
    w <- w/sqrt(sum(w^2))
    t <- as.numeric(E%*%w)
    p <- as.numeric(crossprod(E,t))/sum(t^2)
    q <- sum(t*f)/sum(t^2)
    W[,k] <- w
    P[,k] <- p
    Q[k] <- q
    E <- E-tcrossprod(t,p)
    f <- f-t*q
    rotation <- W[,1:k,drop=FALSE]%*%solve(crossprod(P[,1:k,drop=FALSE],W[,1:k,drop=FALSE]))
    B[,k] <- as.numeric(rotation%*%Q[1:k])/scale
    stopifnot(max(abs(as.numeric((xx%*%B[,k])*sqrt(weights))-(yy*sqrt(weights)-f)))<1e-7)
  }
  list(beta=B,rotation=rotation,scale=scale,W=W,P=P,Q=Q,training_cells=unique(z$depmap_id[rows]))
}
predictions <- selections <- loadings <- list()
jobs <- c('fold_1','fold_2','fold_3','all_data')
for(j in seq_along(jobs)) {
  job <- jobs[j]
  train <- if(j==4) seq_len(nrow(z)) else which(z$fold!=j)
  test <- if(j==4) train else which(z$fold==j)
  cells <- sort(unique(z$depmap_id[train]))
  set.seed(89012+j)
  inner <- setNames(sample(rep(1:3,length.out=length(cells))),cells)
  predictions_inner <- matrix(NA_real_,length(train),max_components)
  for(k in 1:3) {
    held <- which(inner[z$depmap_id[train]]==k)
    model <- fit_pls(train[-held])
    stopifnot(!any(z$depmap_id[train[held]] %in% model$training_cells))
    predictions_inner[held,] <- X[train[held],,drop=FALSE]%*%model$beta
  }
  residual_y <- y[train]
  groups <- paste(z$drug[train],inner[z$depmap_id[train]])
  for(g in unique(groups)) {
    ii <- which(groups==g)
    residual_y[ii] <- residual_y[ii]-mean(residual_y[ii])
    predictions_inner[ii,] <- sweep(predictions_inner[ii,,drop=FALSE],2,colMeans(predictions_inner[ii,,drop=FALSE]))
  }
  scores <- sapply(1:max_components,function(k) mean(sapply(unique(z$drug),function(d) {
    ii <- which(z$drug[train]==d)
    pearson(residual_y[ii],predictions_inner[ii,k])
  })))
  stopifnot(all(is.finite(scores)))
  chosen <- which.max(scores)
  model <- fit_pls(train)
  model$selected_components <- chosen
  if(j<4) stopifnot(!any(z$depmap_id[test] %in% model$training_cells))
  out <- z[test,.(depmap_id,drug,dose,fold,observed_log2=observed_log2_abundance)]
  out[, `:=`(evaluation=if(j==4)'fitted' else 'held_cell', components=chosen,
    predicted_function_score=as.numeric(X[test,,drop=FALSE]%*%model$beta[,chosen]))]
  predictions[[job]] <- out
  selections[[job]] <- data.table(job,components=1:max_components,inner_mean_r=scores,selected=(1:max_components)==chosen)
  rotation <- model$W[,1:chosen,drop=FALSE]%*%solve(crossprod(model$P[,1:chosen,drop=FALSE],model$W[,1:chosen,drop=FALSE]))
  loadings[[job]] <- data.table(job,program=rep(programs,chosen),signature=rep(1:chosen,each=40),weight=as.numeric(rotation),physical_weight=as.numeric(rotation/model$scale))
  qsave(model,file.path(root,'models',paste0(job,'.qs')))
  write_table(data.table(depmap_id=cells,inner_fold=inner[cells]),paste0(job,'_inner_folds'))
  cat(job,':',chosen,'signatures; inner r',round(scores[chosen],3),'\n')
}
pred <- rbindlist(predictions)
write_table(pred,'predictions_frozen')
write_table(rbindlist(selections),'component_selection')
write_table(rbindlist(loadings),'signature_weights')
pred[, observed_variation:=observed_log2-mean(observed_log2),by=.(drug,evaluation,fold)]
pred[, predicted_variation:=predicted_function_score-mean(predicted_function_score),by=.(drug,evaluation,fold)]
metrics <- pred[,.(n=.N,raw_r=pearson(observed_log2,predicted_function_score),within_fold_r=pearson(observed_variation,predicted_variation)),by=.(drug,dose,evaluation)]
write_table(metrics,'correlations')
print(metrics)
# Descriptive fold centering is only for evaluation/plotting, not prediction.
held <- pred[evaluation=='held_cell']
held[, panel:=paste0(drug,' | ',dose,' µM')]
labels <- metrics[evaluation=='held_cell',.(panel=paste0(drug,' | ',dose,' µM'),label=sprintf('r = %.2f; n = %d',within_fold_r,n))]
p <- ggplot(held,aes(observed_variation,predicted_variation))+
  geom_point(aes(color=factor(fold)),size=1.6)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.3,size=3)+
  facet_wrap(~panel,ncol=3,scales='free')+
  scale_color_manual(values=c('#2378a0','#bd5a28','#72419d'),name='Held-cell fold')+
  theme_bw(base_size=10)+theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=8))+
  labs(title='Shared PLS: measured programs to held-cell sensitivity',x='Observed PRISM log2 abundance variation',y='Predicted functional-score variation',
    caption='Drug/dose/fold means removed for display and reported r. Whole cell lines held out.\nDose selection was retrospective; scores are not calibrated absolute viability.')
ggsave(file.path(root,'figures/held_cell_pls.pdf'),p,width=10,height=10)
ggsave(file.path(root,'figures/held_cell_pls.png'),p,width=10,height=10,dpi=120)
writeLines(c('Measured RNA only; shared linear PLS decoder; 1–6 components; no baseline or drug identity input.',
  'Equal drug weight; training-only within-drug centering and scaling. Whole-cell nested component selection.',
  'No combination/HGF or baseline-to-displacement fit in this decoder diagnostic. Dose selection remains retrospective.',
  capture.output(sessionInfo())),file.path(root,'RUN_DETAILS.txt'))
