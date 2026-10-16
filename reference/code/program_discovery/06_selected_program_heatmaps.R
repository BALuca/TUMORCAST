suppressPackageStartupMessages(library(data.table))
source('/work/src/R/lib_heatmap.R')
setDTthreads(2)
ROOT='almanac_conditionals/state_conditioned_bridge/program_first';OUT=file.path(ROOT,'sensitive_resistant/three_panels');selection=fread(file.path(ROOT,'results/program_selection.tsv'));if(!any(selection$selected)){writeLines('No program passed the predeclared association threshold; selected-program heatmaps are not defined.',file.path(OUT,'NO_SELECTED_PROGRAMS.txt'));quit(save='no',status=0)};p=readRDS(file.path(ROOT,'prepared.rds'));stopifnot(!'ACH-000248'%in%p$m$depmap_id);P=colnames(p$W);K=length(P);stopifnot(K>=1)
cluster_rows_safe=function(x)if(nrow(x)>1)hclust(dist(x),method='ward.D2') else FALSE
order_safe=function(x)if(identical(x,FALSE))1L else x$order
file='outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv';hd=names(fread(file,nrows=0));sy=sub(' .*','',hd[-1])
# Read all expression genes, retain the 45 Tahoe recipients, and average duplicate-symbol columns.
e=fread(file);e=e[e[[1]]%in%unique(p$m$depmap_id)];unique_names=make.unique(sy);setnames(e,c('depmap_id',unique_names))
for(g in unique(sy[duplicated(sy)])){cols=unique_names[sy==g];e[,(cols[1]):=rowMeans(as.matrix(.SD)),.SDcols=cols]}
if(anyDuplicated(sy))e[,(unique_names[duplicated(sy)]):=NULL]
candidate_genes=sort(unique(sy));stopifnot(length(candidate_genes)==53960,nrow(e)==45,!anyDuplicated(names(e)))
# Baseline-only eligibility filters, independent of sensitivity values.
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
annotation=as.data.table(as.data.frame(ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,columns=c('gene_id','gene_name','gene_biotype'),return.type='DataFrame')))
ensid=sub('.*\\((ENSG[0-9]+).*','\\1',hd[-1])
mapping=data.table(gene=sy,ensembl_id=ensid);mapping[,biotype:=annotation$gene_biotype[match(ensembl_id,annotation$gene_id)]]
biotypes=mapping[,.(protein_coding=all(!is.na(biotype)&biotype=='protein_coding'),ensembl_ids=paste(ensembl_id,collapse=';'),biotypes=paste(unique(biotype),collapse=';')),by=gene]
TPM=pmax(2^t(as.matrix(e[,candidate_genes,with=FALSE]))-1,0)
gate=data.table(gene=candidate_genes,median_TPM=matrixStats::rowMedians(TPM),fraction_ge_1_TPM=rowMeans(TPM>=1),n_lines=45)
gate=merge(gate,biotypes,by='gene');gate[,eligible:=protein_coding&median_TPM>=5&fraction_ge_1_TPM>=.8];gate[,annotation_source:='EnsDb.Hsapiens.v86; exact Ensembl ID match'];gate[,expression_rule:='Median >=5 TPM and >=1 TPM in >=80% of 45 Tahoe lines']
fwrite(gate,file.path(OUT,'gene_expression_eligibility.tsv'),sep='\t');candidate_genes=gate[eligible==TRUE,gene];stopifnot(length(candidate_genes)>=10)
message('Eligible well-expressed protein-coding genes: ',length(candidate_genes),' of ',nrow(gate))
rm(TPM)


bc=c('#2568A0','#FFFFFF','#BF3434');dc=c('#69459C','#FFFFFF','#D77519')
for(slug in c('lapatinib','crizotinib')){
 dr=if(slug=='lapatinib')'Lapatinib ditosylate' else 'crizotinib';lab=if(slug=='lapatinib')'Lapatinib' else 'Crizotinib'
 d=fread(file.path(OUT,paste0(slug,'_rows.tsv')));d[,group_order:=match(group,c('Sensitive','Resistant','Unclassified'))];setorder(d,group_order,abundance_percent,line,na.last=TRUE);ids=d$depmap_id
 stopifnot(nrow(d)==45,all(d[group=='Sensitive',abundance_percent]<=50),all(d[group=='Resistant',abundance_percent]>50),all(is.na(d[group=='Unclassified',abundance_percent])))
 for(gr in c('Sensitive','Resistant'))stopifnot(all(diff(d[group==gr,abundance_percent])>=0))
 fwrite(d[,.(depmap_id,line,group,abundance_percent,column_order=seq_len(.N))],file.path(OUT,paste0(slug,'_rows.tsv')),sep='\t')
 A=zscale_rows(t(as.matrix(e[match(ids,e$depmap_id),c('ERBB2','ERBB3','EGFR','MET','HGF'),with=FALSE])));rownames(A)=c('HER2 / ERBB2','ERBB3','EGFR','MET','HGF')
 B=zscale_rows(t(as.matrix(p$controls[match(ids,p$controls$depmap_id),..P])))
 D=do.call(rbind,lapply(c(.05,.5,5),function(ds){x=p$m[drug==dr&dose==ds];t(as.matrix(x[match(ids,x$depmap_id),..P]))}));rownames(D)=rep(P,3)
 raw_D=D;D=zscale_rows(raw_D)
 stopifnot(max(abs(rowMeans(D)))<1e-10,max(abs(apply(D,1,sd)-1))<1e-10)
 stats=data.table(program=rep(P,3),dose=rep(c(.05,.5,5),each=K),raw_mean=rowMeans(raw_D),raw_sd=apply(raw_D,1,sd),n=45)
 fwrite(stats,file.path(OUT,paste0(slug,'_displacement_scaling.tsv')),sep='\t')
 raw_export=as.data.table(raw_D);setnames(raw_export,ids);fwrite(cbind(stats[,.(program,dose)],raw_export),file.path(OUT,paste0(slug,'_displacement_raw.tsv')),sep='\t')
 z_export=as.data.table(D);setnames(z_export,ids);fwrite(cbind(stats[,.(program,dose)],z_export),file.path(OUT,paste0(slug,'_displacement_z.tsv')),sep='\t')
 for(n in c('A','B','D')){x=get(n);colnames(x)=ids;assign(n,x)}
 stopifnot(ncol(A)==45,identical(colnames(A),colnames(B)),identical(colnames(B),colnames(D)),all(is.finite(c(A,B,D))))
 # Exploratory genome-wide ranking; sensitivity is 100 minus measured percent abundance.
 GX=t(as.matrix(e[match(ids,e$depmap_id),candidate_genes,with=FALSE]));colnames(GX)=ids
 sensitivity=100-d$abundance_percent;known=is.finite(sensitivity);stopifnot(sum(known)==33,all(d$group[!known]=='Unclassified'))
 X=GX[,known,drop=FALSE];y=sensitivity[known];centered=sweep(X,1,rowMeans(X),'-');yc=y-mean(y);den=sqrt(rowSums(centered^2)*sum(yc^2));r=as.numeric(centered%*%yc)/den;r[!is.finite(r)]=NA_real_;r=pmax(-1,pmin(1,r));pv=2*pt(abs(r)*sqrt((sum(known)-2)/pmax(1-r^2,.Machine$double.eps)),df=sum(known)-2,lower.tail=FALSE)
 corr=data.table(gene=rownames(GX),r=r,p=pv,n=sum(known),baseline_sd=apply(X,1,sd));corr[,q:=p.adjust(p,'BH')];corr=corr[order(-abs(r),gene,na.last=TRUE)];corr[,rank:=seq_len(.N)];corr[,selected:=rank<=10];corr[,sensitivity_definition:='100 - measured abundance percent at 2.5 uM'];corr[,status:='exploratory well-expressed protein-coding genome-wide selection; not validated prediction'];corr[,test_status:=ifelse(is.finite(r),'tested','constant or nonfinite expression')]
 stopifnot(nrow(corr)==length(candidate_genes),all(is.finite(corr[1:10,r])))
 # Check the vectorized calculation against cor.test for all displayed genes.
 for(g in corr[1:10,gene]){ct=cor.test(GX[g,known],y);stopifnot(abs(unname(ct$estimate)-corr[gene==g,r])<1e-10,abs(ct$p.value-corr[gene==g,p])<1e-10)}
 message(slug,': ',sum(is.finite(corr$r)),' variable genes tested from ',nrow(corr),' unique symbols')

 fwrite(corr,file.path(OUT,paste0(slug,'_gene_sensitivity_correlations.tsv')),sep='\t')
 topgenes=corr[1:10,gene];G=zscale_rows(GX[topgenes,,drop=FALSE]);hcG=hclust(dist(G),method='ward.D2');gene_labels=sprintf('%s  (r=%+.2f)',topgenes,corr[match(topgenes,gene),r]);stopifnot(identical(colnames(G),ids))
 fwrite(as.data.table(GX[topgenes,,drop=FALSE],keep.rownames='gene'),file.path(OUT,paste0(slug,'_top10_genes_raw.tsv')),sep='\t');fwrite(as.data.table(G,keep.rownames='gene'),file.path(OUT,paste0(slug,'_top10_genes_z.tsv')),sep='\t')
 # Cluster standardized row profiles, without using sensitivity labels.
 hcA=hclust(dist(A),method='ward.D2');hcB=cluster_rows_safe(B)
 jointD=do.call(cbind,lapply(0:2,function(k)D[k*K+seq_len(K),,drop=FALSE]));rownames(jointD)=P
 hcD=cluster_rows_safe(jointD)
 average_raw=Reduce('+',lapply(0:2,function(k)raw_D[k*K+seq_len(K),,drop=FALSE]))/3;rownames(average_raw)=P;colnames(average_raw)=ids
 E=zscale_rows(average_raw);hcE=cluster_rows_safe(E);stopifnot(max(abs(rowMeans(E)))<1e-10,max(abs(apply(E,1,sd)-1))<1e-10,identical(colnames(E),ids))
 fwrite(as.data.table(average_raw,keep.rownames='program'),file.path(OUT,paste0(slug,'_dose_average_raw.tsv')),sep='\t');fwrite(as.data.table(E,keep.rownames='program'),file.path(OUT,paste0(slug,'_dose_average_z.tsv')),sep='\t')
 orders=rbindlist(list(data.table(panel='A baseline features',feature=rownames(A)[hcA$order]),data.table(panel='B baseline programs',feature=rownames(B)[order_safe(hcB)]),data.table(panel='C displacement (all doses)',feature=P[order_safe(hcD)]),data.table(panel='D mean across doses',feature=P[order_safe(hcE)]),data.table(panel='E top10 baseline genes',feature=topgenes[hcG$order])));orders[,position:=seq_len(.N),by=panel]
 fwrite(orders,file.path(OUT,paste0(slug,'_row_clustering.tsv')),sep='\t');saveRDS(list(A=hcA,B=hcB,C=hcD,D=hcE,E=hcG,method='ward.D2',distance='Euclidean on standardized row profiles',C_input='Selected programs x 135 line-dose values'),file.path(OUT,paste0(slug,'_row_clustering.rds')))
 # Row-wise association with the same measured continuous sensitivity endpoint.
 rowstats=function(X,panel,doses=rep(NA_real_,nrow(X))){rbindlist(lapply(seq_len(nrow(X)),function(j){ct=cor.test(X[j,known],y,method='pearson');data.table(panel,feature=rownames(X)[j],dose=doses[j],n=sum(known),r=unname(ct$estimate),p=ct$p.value)}))}
 assoc=rbindlist(list(rowstats(A,'A'),rowstats(B,'B'),rowstats(raw_D,'C',rep(c(.05,.5,5),each=K)),rowstats(average_raw,'D')))
 stopifnot(nrow(assoc)==5+5*K);assoc[,q:=p.adjust(p,'BH')];assoc[,fdr_family:=sprintf('%d feature/program/dose rows; program selection not adjusted',5+5*K)]
 sg=corr[match(topgenes,gene),.(panel='E',feature=gene,dose=NA_real_,n,r,p,q,fdr_family='9063 eligible protein-coding genes per drug')];assoc=rbind(assoc,sg);assoc[,neglog10p:=-log10(pmax(p,.Machine$double.xmin))];assoc[,significant:=q<.05];assoc[,drug:=dr]
 fwrite(assoc,file.path(OUT,paste0(slug,'_all_row_sensitivity_statistics.tsv')),sep='\t')
 sa=assoc[panel=='A'];sb=assoc[panel=='B'];sc=assoc[panel=='C'];sdose=assoc[panel=='D'];se=assoc[panel=='E']
 stopifnot(identical(sa$feature,rownames(A)),identical(sb$feature,rownames(B)),identical(sc$feature,rownames(D)),identical(sdose$feature,rownames(E)),identical(se$feature,rownames(G)))
 pmaxaxis=max(2,ceiling(max(assoc$neglog10p)))
 bars=function(st,panel){
 names_r=paste0('r_',panel);names_p=paste0('p_',panel)
 args=list(anno_barplot(st$r,which='row',baseline=0,ylim=c(-1,1),extend=0,bar_width=.75,gp=gpar(fill=ifelse(st$r>=0,'#BF3434','#2568A0'),col=NA),border=TRUE,width=unit(23,'mm'),axis_param=list(side='bottom',at=c(-1,0,1),gp=gpar(fontsize=8))),anno_barplot(st$neglog10p,which='row',baseline=0,ylim=c(0,pmaxaxis),extend=0,bar_width=.75,gp=gpar(fill=ifelse(st$significant,'#69459C','#BDBDBD'),col=NA),border=TRUE,width=unit(23,'mm'),axis_param=list(side='bottom',at=unique(c(0,2,pmaxaxis)),gp=gpar(fontsize=8))))
 names(args)=c(names_r,names_p)
 do.call(rowAnnotation,c(args,list(annotation_label=c('Pearson r','-log10(p)'),annotation_name_side='top',annotation_name_rot=0,annotation_name_gp=gpar(fontsize=9,fontface='bold'),gap=unit(4,'mm'))))
 }
 starlabels=function(st)paste0(st$feature,ifelse(st$significant,' *',''))
 lev=c('Sensitive','Resistant','Unclassified');split=factor(paste0(d$group,' (n=',as.integer(table(d$group)[d$group]),')'),levels=paste0(lev,' (n=',as.integer(table(factor(d$group,levels=lev))),')'));lim=3
 common=list(scale_rows=FALSE,row_dend_reorder=FALSE,show_row_dend=TRUE,row_dend_width=unit(12,'mm'),cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=split,column_gap=unit(3,'mm'),show_row_names=TRUE,row_names_side='left',row_names_gp=gpar(fontsize=10),column_labels=d$line,column_names_gp=gpar(fontsize=9),column_names_rot=90,border=TRUE,width=unit(267,'mm'),show_heatmap_legend=FALSE)
 a=do.call(draw_heatmap,c(list(mat=A,right_annotation=bars(sa,'A'),row_labels=starlabels(sa),cluster_rows=hcA,hm_name='A',show_column_names=FALSE,title=levels(split),column_title_gp=gpar(fontsize=11,fontface='bold'),value_range=c(-3,0,3),value_colors=bc,height=unit(25,'mm')),common))
 b=do.call(draw_heatmap,c(list(mat=B,right_annotation=bars(sb,'B'),row_labels=starlabels(sb),cluster_rows=hcB,hm_name='B',show_column_names=FALSE,title=NULL,value_range=c(-3,0,3),value_colors=bc,height=unit(4.7*K,'mm')),common))
 cc=do.call(draw_heatmap,c(list(mat=D,right_annotation=bars(sc,'C'),row_labels=starlabels(sc),cluster_rows=if(K>1)function(x)hcD else FALSE,hm_name='C',show_column_names=FALSE,title=NULL,value_range=c(-lim,0,lim),value_colors=dc,height=unit(12.9*K,'mm'),row_split=factor(rep(c('0.05 uM','0.5 uM','5 uM'),each=K),levels=c('0.05 uM','0.5 uM','5 uM')),cluster_row_slices=FALSE,row_gap=unit(2,'mm'),row_title_gp=gpar(fontsize=10,fontface='bold')),common))
 ee=do.call(draw_heatmap,c(list(mat=E,right_annotation=bars(sdose,'D'),row_labels=starlabels(sdose),cluster_rows=hcE,hm_name='D',show_column_names=FALSE,title=NULL,value_range=c(-lim,0,lim),value_colors=dc,height=unit(4.7*K,'mm')),common))
 gg=do.call(draw_heatmap,c(list(mat=G,right_annotation=bars(se,'E'),cluster_rows=hcG,hm_name='E',row_labels=starlabels(se),show_column_names=TRUE,title=NULL,value_range=c(-3,0,3),value_colors=bc,height=unit(47,'mm')),common))
 ymax=ceiling(max(d$abundance_percent,na.rm=TRUE)/25)*25
 barcols=unname(c(Sensitive='#B44745',Resistant='#426A89',Unclassified='#BBBBBB')[d$group])
 ha=HeatmapAnnotation(Abundance=anno_barplot(d$abundance_percent,ylim=c(0,ymax),extend=0,gp=gpar(fill=barcols,col=NA),border=TRUE,bar_width=.8,height=unit(20,'mm'),axis_param=list(side='left',at=unique(c(0,50,100,ymax)),labels=unique(c(0,50,100,ymax)),gp=gpar(fontsize=8))),annotation_label='Remaining cells (%)\n2.5 uM; lower = more sensitive',annotation_name_side='left',annotation_name_rot=0,annotation_name_gp=gpar(fontsize=9))
 h=ha%v%a%v%b%v%cc%v%ee%v%gg
 legends=list(Legend(title='A + B + E: baseline level relative to other lines',col_fun=colorRamp2(c(-3,0,3),bc),at=c(-3,0,3),labels=c('Lower (-3 SD)','Mean (0)','Higher (+3 SD)'),direction='horizontal',legend_width=unit(55,'mm')),Legend(title='C + D: displacement relative to other lines',col_fun=colorRamp2(c(-lim,0,lim),dc),at=c(-lim,0,lim),labels=c('Lower (-3 SD)','Mean (0)','Higher (+3 SD)'),direction='horizontal',legend_width=unit(65,'mm')))
 paint=function(){drawn=draw(h,ht_gap=unit(c(13,11,11,11,11),'mm'),heatmap_legend_list=list(packLegend(list=legends,direction='horizontal',gap=unit(15,'mm'))),heatmap_legend_side='bottom',padding=unit(c(25,3,24,3),'mm'))
 actual=row_order(drawn,name='C');stopifnot(length(actual)==3);for(k in 1:3)stopifnot(identical(rownames(D)[actual[[k]]],P[order_safe(hcD)]))
 for(j in 1:2)decorate_annotation('Abundance',slice=j,{grid.lines(x=c(0,1),y=rep(50/ymax,2),gp=gpar(lty=2,col='#222222',lwd=1))})
 decorate_annotation('Abundance',slice=3,{grid.text('No matched measurement',gp=gpar(fontsize=9,col='#555555'))})
 grid.text(paste0(lab,' | Programs discovered first | All 45 lines; AU565 excluded'),x=.5,y=unit(1,'npc')-unit(7,'mm'),gp=gpar(fontsize=17,fontface='bold'))
 decorate_column_title('A',slice=1,{grid.text('A  Baseline features',x=unit(0,'npc'),y=unit(8,'mm'),just=c('left','bottom'),gp=gpar(fontsize=13,fontface='bold'))})
 decorate_heatmap_body('B',column_slice=1,{grid.text('B  Baseline programs | selected from the full-gene program bank',x=0,y=unit(1,'npc')+unit(3,'mm'),just=c('left','bottom'),gp=gpar(fontsize=13,fontface='bold'))})
 decorate_heatmap_body('C',row_slice=1,column_slice=1,{grid.text('C  Measured displacement | standardized within each program and dose',x=0,y=unit(1,'npc')+unit(3,'mm'),just=c('left','bottom'),gp=gpar(fontsize=13,fontface='bold'))})
 decorate_heatmap_body('D',column_slice=1,{grid.text('D  Mean displacement across doses | average raw changes, then standardize',x=0,y=unit(1,'npc')+unit(3,'mm'),just=c('left','bottom'),gp=gpar(fontsize=13,fontface='bold'))})
 decorate_heatmap_body('E',column_slice=1,{grid.text('E  Top 10 well-expressed protein-coding genes | sensitivity correlation (n=33)',x=0,y=unit(1,'npc')+unit(3,'mm'),just=c('left','bottom'),gp=gpar(fontsize=13,fontface='bold'))})
 grid.text('Discovery associations: program selection used these outcomes. Side bars: r and -log10(p); purple/*: nominal BH q<0.05, not selection-adjusted.',x=.5,y=unit(20,'mm'),gp=gpar(fontsize=9))
 grid.text('Sensitivity = 100 - remaining cells (%). B-D show selected bank programs; E remains baseline genes. Protein-coding, median >=5 TPM, >=1 TPM in >=80% of lines.',x=.5,y=unit(15,'mm'),gp=gpar(fontsize=9))
 grid.text('Each row standardized across all 45 lines; colors saturate at +/-3 SD. White in C/D means mean displacement, not zero change.',x=.5,y=unit(10,'mm'),gp=gpar(fontsize=9))
 grid.text('Sensitive: <=50% abundance at 2.5 uM; resistant: >50%; unclassified: no matched PRISM outcome.',x=.5,y=unit(5,'mm'),gp=gpar(fontsize=10))}
 pdf(file.path(OUT,paste0(slug,'_annotated_sensitivity.pdf')),width=16,height=9.72+0.878*K,useDingbats=FALSE);paint();dev.off()
 png(file.path(OUT,paste0(slug,'_annotated_sensitivity.png')),width=16,height=9.72+0.878*K,units='in',res=160);paint();dev.off()
}
