#!/usr/bin/env Rscript
## Untouched secondary cohort: GFP-vector T47D samples from GSE153509.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
iteration <- sub('^--iteration=','',args[grepl('^--iteration=',args)])
if(!length(iteration)) iteration <- '04'
iteration <- sprintf('%02d',as.integer(iteration))
root <- '/work/tahoe_prism_v2'; out <- file.path(root,'results/28_fgf2_fulvestrant_bypass',paste0('iteration_',iteration))
fig <- file.path(out,'figures'); frozen <- readRDS(file.path(root,'artifacts/28_fgf2_fulvestrant_bypass',paste0('iteration_',iteration),'frozen_programs.rds'))
memb <- as.data.table(frozen$membership)[order(plotting_order)]; genes <- memb$gene
orient <- setNames(memb$orientation,memb$gene); roles <- c('ER','FGF','Root')
labels <- c(ER='ER branch',FGF='FGF receptor branch',Root='Shared proliferation root')
cols <- c(ER='#7A3E9D',FGF='#D97818',Root='#167C80')
meta <- fread(file.path(root,'results/28_fgf2_fulvestrant_bypass/28a_gse153509_complete_sample_design.tsv'))[construct=='GFP']
raw <- fread(file.path(root,'artifacts/28_fgf2_fulvestrant_bypass/GSE153509_FGFR.Rescue.seq.counts.QC.csv.gz'));setnames(raw,1,'gene');raw<-raw[gene %chin% genes]
sc <- intersect(meta$matrix_column,names(raw));meta<-meta[match(sc,matrix_column)]
counts<-as.matrix(raw[,..sc]);rownames(counts)<-raw$gene
lc<-cpm(calcNormFactors(DGEList(counts=counts)),log=TRUE,prior.count=2);lc<-lc[genes,,drop=FALSE]
z<-t(scale(t(lc)));z[!is.finite(z)]<-0;z<-sweep(z,1,orient[rownames(z)],'*')
scores<-sapply(roles,function(r)colMeans(z[memb[role==r,gene],,drop=FALSE]))
contrasts<-data.table(contrast=c('FGF2 rescue under fulvestrant','Palbociclib blocks FGF2 rescue'),
 treated=c('Fulvestrant + FGF2','Fulvestrant + Palbociclib + FGF2'),reference=c('Fulvestrant','Fulvestrant + FGF2'))
res<-rbindlist(lapply(seq_len(nrow(contrasts)),function(i)rbindlist(lapply(seq_along(roles),function(j){
 a<-scores[meta$arm==contrasts$treated[i],j];b<-scores[meta$arm==contrasts$reference[i],j];tt<-t.test(a,b)
 data.table(contrast=contrasts$contrast[i],role=unname(labels[roles[j]]),effect=mean(a)-mean(b),
  ci_low=tt$conf.int[1],ci_high=tt$conf.int[2],p_value=tt$p.value,n_treated=length(a),n_reference=length(b))
}))))
res[,q_value:=p.adjust(p_value,'BH')];fwrite(res,file.path(out,'28c_gfp_holdout_contrasts.tsv'),sep='\t')
criteria<-res[(contrast=='FGF2 rescue under fulvestrant' & role %chin% labels[c('FGF','Root')])|
 (contrast=='Palbociclib blocks FGF2 rescue' & role==labels['Root'])]
criteria[,expected_sign:=fifelse(contrast=='Palbociclib blocks FGF2 rescue',-1,1)]
criteria[,passed:=sign(effect)==expected_sign&q_value<.05]
fwrite(criteria,file.path(out,'28c_gfp_holdout_primary_criteria.tsv'),sep='\t')
res[,contrast:=factor(contrast,levels=rev(contrasts$contrast))]
p<-ggplot(res,aes(effect,contrast,colour=role))+geom_vline(xintercept=0,colour='#777777',linewidth=.45)+
 geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=.16,linewidth=.75)+
 geom_point(aes(shape=q_value<.05),size=3)+facet_wrap(~role,nrow=1,scales='free_x')+
 scale_colour_manual(values=setNames(cols,unname(labels)))+scale_shape_manual(values=c(`TRUE`=19,`FALSE`=1))+
 labs(title='Untouched GFP-vector T47D validation cohort',x='Difference in frozen program score',y=NULL)+
 theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),axis.title=element_text(face='bold'),
  strip.background=element_rect(fill='#F2F2F2',colour=NA),strip.text=element_text(face='bold'),legend.position='none',aspect.ratio=1)
ggsave(file.path(fig,'28c_06_gfp_holdout_contrasts.pdf'),p,width=8.6,height=3.5,device=cairo_pdf)
ggsave(file.path(fig,'28c_06_gfp_holdout_contrasts.png'),p,width=8.6,height=3.5,dpi=220,bg='white')
