#!/usr/bin/env Rscript
# Joint RNA is a measured-information diagnostic, not an input to held BT20
# predictions. Reuse the original40 projections already calculated in61.
suppressPackageStartupMessages({library(data.table); library(qs); library(limma); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/65_sparse_program_interaction_patterns'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
source_root <- '/work/results/program_pilot/59_mechanistic_evidence_datasets/data'
gene_root <- '/work/results/program_pilot/60_defensible_response_composition/data'
setDTthreads(2L)
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
save_plot <- function(p,name,width,height) {
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=width,height=height,dpi=150)
}
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
input_paths <- c(file.path(prior,c('original40_Diaz_scores.tsv.gz','original40_program_annotations.tsv')),
  file.path(source_root,'diaz_functional_factorials.tsv.gz'),
  file.path(gene_root,paste0('representation_',c('MCF7','LNCaP'),'_gene_interactions.tsv.gz')))
save_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),
  'joint_RNA_source_manifest')
scores <- fread(file.path(prior,'original40_Diaz_scores.tsv.gz'))[measurement=='full' & time_h<=24]
annotations <- fread(file.path(prior,'original40_program_annotations.tsv'))
function_data <- fread(file.path(source_root,'diaz_functional_factorials.tsv.gz'))[
  assay=='InCell' & arm_AB %chin% c('TM','TW','MW') & time_h<=24]
function_data[,`:=`(independent_pct=A_pct*B_pct/100,
  functional_departure_log2=log2(AB_pct/(A_pct*B_pct/100)),
  extra_inhibition_pp=A_pct*B_pct/100-AB_pct)]
stopifnot(nrow(function_data)==36L,all(function_data$AB_pct>0))
programs <- sprintf('P%02d',1:40)
effects <- gene_counts <- list()
for(cell in c('MCF7','LNCaP')) {
  s <- scores[cell_line==cell]
  md <- unique(s[,.(sample_id,arm,time_h)])
  wide <- dcast(s,program~sample_id,value.var='score')
  columns <- md$sample_id
  mat <- as.matrix(wide[match(programs,program),..columns])
  rownames(mat) <- programs
  stopifnot(identical(colnames(mat),md$sample_id),all(is.finite(mat)))
  group <- factor(paste(md$arm,md$time_h,sep='_'))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  cc <- list()
  cm <- list()
  for(tt in c(0,3,6,9,12,24)) for(pair in c('TM','TW','MW')) {
    aa <- strsplit(pair,'')[[1]]
    terms <- paste(c(pair,aa,'DMSO'),tt,sep='_')
    stopifnot(all(terms %in% colnames(design)))
    v <- setNames(numeric(ncol(design)),colnames(design))
    v[terms] <- c(1,-1,-1,1)
    id <- paste(pair,tt,sep='_')
    cc[[id]] <- v
    cm[[id]] <- data.table(contrast_id=id,arm_AB=pair,time_h=tt)
  }
  fit <- eBayes(contrasts.fit(lmFit(mat,design),do.call(cbind,cc)),trend=TRUE,robust=TRUE)
  # RMS uses individual molecular responses only, without functional labels.
  single_means <- s[arm %chin% c('T','M','W') & time_h>0,.(value=mean(score)),by=.(program,arm,time_h)]
  scale <- single_means[,.(single_rms=sqrt(mean(value^2))),by=program]
  stopifnot(all(scale$single_rms>0))
  effect <- rbindlist(lapply(seq_len(ncol(fit$coefficients)),function(j) data.table(
    program=programs,contrast_id=colnames(fit$coefficients)[j],effect=fit$coefficients[,j],
    standard_error=fit$stdev.unscaled[,j]*sqrt(fit$s2.post),p_value=fit$p.value[,j],
    q_value=p.adjust(fit$p.value[,j],'BH'))))
  effect <- merge(merge(effect,rbindlist(cm),by='contrast_id'),scale,by='program')
  effect[,`:=`(cell_line=cell,scaled_departure=effect/single_rms)]
  effects[[cell]] <- effect
  genes <- fread(file.path(gene_root,paste0('representation_',cell,'_gene_interactions.tsv.gz')))[
    contrast_type=='pair_interaction']
  gene_counts[[cell]] <- genes[,.(clear_genes=sum(clear_interaction)),by=.(cell_line,arm,time_h)]
}
effect <- rbindlist(effects)
effect <- merge(effect,function_data[,.(cell_line,arm_AB,time_h,AB_pct,independent_pct,
  extra_inhibition_pp,functional_departure_log2)],by=c('cell_line','arm_AB','time_h'))
effect <- merge(effect,annotations[,.(program,display_name)],by='program')
save_table(effect,'joint_RNA_original40_program_departures')
gene_counts <- rbindlist(gene_counts)
setnames(gene_counts,'arm','arm_AB')
summary <- effect[,.(molecular_departure_rms=sqrt(mean(scaled_departure^2)),
  significant_programs=sum(q_value<=.05),large_significant_programs=sum(q_value<=.05 & abs(scaled_departure)>=.5)),
  by=.(cell_line,arm_AB,time_h,AB_pct,independent_pct,extra_inhibition_pp,functional_departure_log2)]
summary <- merge(summary,gene_counts,by=c('cell_line','arm_AB','time_h'))
save_table(summary,'joint_RNA_functional_diagnostic')
pair_names <- c(TM='Tamoxifen + mefloquine',TW='Tamoxifen + withaferin A',MW='Mefloquine + withaferin A')

# Same ordering across the measured counts, gene departures and program map.
for(cell in c('MCF7','LNCaP')) {
  d <- summary[cell_line==cell]
  curves <- rbind(d[,.(arm_AB,time_h,value=AB_pct,measurement='Measured combination')],
    d[,.(arm_AB,time_h,value=independent_pct,measurement='Independent individual effects')])
  curves[,pair:=pair_names[arm_AB]]
  p <- ggplot(curves,aes(time_h,value,color=measurement))+geom_line(linewidth=.8)+geom_point(size=1.8)+
    facet_wrap(~pair,nrow=1)+scale_color_manual(values=c('Measured combination'='black','Independent individual effects'='#777777'))+
    labs(title=paste(cell,': functional departure from the individual effects'),
      subtitle='Measured nuclei counts; exact same drugs, doses and sampling times as the RNA comparisons',
      x='Time (h)',y='Nuclei (% vehicle)',color=NULL,
      caption='Tamoxifen20 µM; mefloquine10 µM; withaferin A5 µM. Function and RNA are separate wells.')
  save_plot(p,paste0('joint_RNA_',cell,'_functional_time_courses'),10.5,4)
  ee <- copy(effect[cell_line==cell])
  ee[,column_id:=paste(arm_AB,sprintf('%02dh',time_h),sep=' | ')]
  wide <- dcast(ee,program~column_id,value.var='scaled_departure')
  mat <- as.matrix(wide[,-1])
  rownames(mat) <- wide$program
  order <- hclust(dist(mat),method='ward.D2')$order
  mat <- mat[order,,drop=FALSE]
  md <- unique(ee[,.(column_id,arm_AB,time_h,extra_inhibition_pp)])
  md <- md[match(colnames(mat),column_id)]
  stopifnot(identical(colnames(mat),md$column_id))
  label <- substr(gsub('GOBP ','',annotations$display_name[match(rownames(mat),annotations$program)],fixed=TRUE),1,48)
  limit <- max(abs(mat))
  hm <- draw_heatmap(mat,col_anno_df=as.data.frame(md),col_anno_cols='arm_AB',
    ann_top_args=list('Extra inhibition (pp)'=ComplexHeatmap::anno_barplot(md$extra_inhibition_pp,
      height=grid::unit(1.3,'cm'),gp=grid::gpar(fill='#555555'),border=FALSE),
      annotation_name_gp=grid::gpar(fontsize=8)),
    scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    hm_name='Joint RNA\ndeparture',show_row_names=TRUE,show_column_names=TRUE,row_labels=label,
    row_names_side='left',row_names_max_width=grid::unit(6,'cm'),
    row_names_gp=grid::gpar(fontsize=6.6),column_names_gp=grid::gpar(fontsize=7),column_names_rot=70,
    title=paste(cell,'| original40 | combination minus both singles plus vehicle'))
  cairo_pdf(file.path(root,'figures',paste0('joint_RNA_',cell,'_program_departures.pdf')),width=11,height=8.5)
  ComplexHeatmap::draw(hm,heatmap_legend_side='bottom',annotation_legend_side='bottom')
  dev.off()
  png(file.path(root,'figures',paste0('joint_RNA_',cell,'_program_departures.png')),width=1650,height=1275,res=150)
  ComplexHeatmap::draw(hm,heatmap_legend_side='bottom',annotation_legend_side='bottom')
  dev.off()
}

# Plot all trajectories, not just the most responsive combination. This is a
# descriptive association: repeated times share biology and are not independent.
d <- summary[time_h>0]
d[,pair:=pair_names[arm_AB]]
labels <- d[,.(label=sprintf('r = %.2f; %d pair-times',cor(molecular_departure_rms,extra_inhibition_pp),.N)),by=cell_line]
p <- ggplot(d,aes(molecular_departure_rms,extra_inhibition_pp))+
  geom_path(aes(group=pair,color=pair),alpha=.35)+geom_point(aes(color=pair,size=time_h),alpha=.8)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.55)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.2,size=3.5)+
  facet_wrap(~cell_line,nrow=1,scales='free')+theme(aspect.ratio=1)+
  scale_size_continuous(breaks=c(3,6,12,24))+
  labs(title='Do measured joint molecular departures accompany extra functional inhibition?',
    subtitle='Separate Diaz diagnostic; observed combination RNA is never supplied to BT20 predictions',
    x='RMS program departure (scaled by single-response magnitude)',y='Extra inhibition beyond individual effects (pp)',
    color='Combination',size='Time (h)',caption='Descriptive time courses; association does not establish causation or prediction from individual responses.')
save_plot(p,'joint_RNA_departure_vs_function',10,5.3)
top <- effect[time_h>0 & q_value<=.05][order(-abs(scaled_departure)),head(.SD,5),by=.(cell_line,arm_AB,time_h)]
save_table(top,'joint_RNA_largest_program_departures')
writeLines(c('JOINT MOLECULAR INFORMATION DIAGNOSTIC',
  'BT20 has individual RNA and combination function, but no measured combination RNA.',
  'Diaz provides the required factorial molecular measurements in MCF7 and LNCaP.',
  'Original40 frozen full-RNA scores reused exactly from61; no program selection on function.',
  'Molecular departure is AB-A-B+vehicle on log2(1+TMM CP10K), projected with frozen signed weights.',
  'Plot scaling uses RMS individual responses. Limma uncertainty uses deposited wells; biological independence unresolved.',
  'Functional departure is individual count-product minus measured joint count, in percentage points.',
  'All three pairs at all six matched times retained. Time-course correlations are descriptive, not independent tests.',
  'This checks whether joint expression contains additional observed information. It does not establish that',
  'BT20 interaction patterns predict those changes, nor that these molecular changes cause loss of viable cells.',
  '',capture.output(print(summary[time_h %in% c(3,24)]))),file.path(root,'JOINT_RNA_DIAGNOSTIC.txt'))
print(summary[time_h %in% c(3,24)])
