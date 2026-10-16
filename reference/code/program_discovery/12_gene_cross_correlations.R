suppressPackageStartupMessages(library(data.table));source('/work/src/R/lib_heatmap.R');setDTthreads(2)
BASE='almanac_conditionals/state_conditioned_bridge/program_first';OUT=file.path(BASE,'figures/cross_correlations');dir.create(OUT,recursive=T,showWarnings=F)
genes=fread(file.path(BASE,'results/selected_gene_universe.tsv'))$gene;m=fread(file.path(BASE,'results/frozen_program_weights.tsv'));programs=sort(unique(m$program));stopifnot(setequal(genes,m$gene),!anyDuplicated(m$gene));N=length(genes)
read_npy=function(path){con=file(path,'rb');on.exit(close(con));stopifnot(rawToChar(readBin(con,'raw',6))=='\x93NUMPY');v=readBin(con,'integer',2,size=1,signed=F);n=readBin(con,'integer',1,size=if(v[1]==1)2 else 4,signed=F,endian='little');h=rawToChar(readBin(con,'raw',n));dims=as.integer(strsplit(sub(".*'shape': \\(([^)]*)\\).*",'\\1',h),',')[[1]]);dims=dims[!is.na(dims)];stopifnot(identical(dims,c(N,N)));x=matrix(readBin(con,'numeric',n=prod(dims),size=if(grepl('f8',h))8 else 4,endian='little'),nrow=N,byrow=!grepl('True',h));dimnames(x)=list(genes,genes);x}
C=read_npy(file.path(BASE,'artifacts/discovery_joint_correlation.npy'))
# Standard pipeline ordering: frozen blocks, hierarchical clustering within each
# block using the sign-oriented joint co-response correlation.
ord=unlist(lapply(programs,function(p){g=m[program==p,gene];s=m[match(g,gene),orientation];x=C[g,g,drop=F]*tcrossprod(s);if(length(g)>1)g[hclust(as.dist(pmax(1-x,0)),method='average')$order] else g}),use.names=F)
stopifnot(length(ord)==N,!anyDuplicated(ord),setequal(ord,genes));mi=m[match(ord,gene)];stopifnot(identical(mi$gene,ord));pg=factor(mi$program,levels=programs);sizes=as.numeric(table(pg));edges=cumsum(sizes)/N;edges=edges[-length(edges)];mid=cumsum(sizes)-floor(sizes/2);colnames(C)=rownames(C)=genes
pcols=setNames(hcl.colors(length(programs),'Dark 3'),programs);key=m[order(program,-absolute_weight),.(label=paste0(program,'  ',paste(head(gene,2),collapse=', '))),by=program];labels=key$label[match(programs,key$program)]
fwrite(data.table(position=seq_len(N),gene=ord,program=as.character(pg),orientation=mi$orientation),file.path(BASE,'results/cross_correlation_gene_order.tsv'),sep='\t')
# Display aggregation only: contiguous bins of <=8 genes within each program.
# Every matrix entry contributes; original correlations and gene order remain intact.
bins=integer(N);offset=0L
for(pr in programs){ix=which(pg==pr);local=ceiling(seq_along(ix)/8);bins[ix]=offset+local;offset=offset+max(local)}
bs=as.numeric(tabulate(bins));bpg=factor(as.character(pg[match(seq_along(bs),bins)]),levels=programs);display_ids=paste0('bin_',seq_along(bs));sizes=as.numeric(table(bpg));edges=head(cumsum(sizes)/sum(sizes),-1);mid=cumsum(sizes)-floor(sizes/2)
fwrite(data.table(position=seq_len(N),gene=ord,program=as.character(pg),display_bin=bins),file.path(BASE,'results/cross_correlation_display_bins.tsv'),sep='\t')
lim=.35;colours=c('#244B6B','#9DB8C8','#F7F5F0','#D99A83','#9E3D43')
render=function(mat,slug,title,oriented=T,png_preview=F){
 x=mat[ord,ord,drop=F];if(oriented)x=sweep(sweep(x,1,mi$orientation,'*'),2,mi$orientation,'*')
 stopifnot(identical(rownames(x),ord),identical(colnames(x),ord),all(is.finite(x)))
 marks=rowAnnotation(Programs=anno_mark(at=mid,labels=labels,labels_gp=gpar(fontsize=6.5),link_width=unit(2,'mm'),padding=unit(.4,'mm')),show_annotation_name=F,width=unit(42,'mm'))
 x=rowsum(x,bins,reorder=FALSE)/bs;x=t(rowsum(t(x),bins,reorder=FALSE));x=sweep(x,2,bs,'/');dimnames(x)=list(display_ids,display_ids);stopifnot(nrow(x)==length(bs),sum(bs)==N,max(bs)<=8)
 ann=data.frame(Program=bpg,row.names=display_ids)
 hm=draw_heatmap(x,hm_name=if(oriented)'Oriented correlation' else 'Raw correlation',scale_rows=F,value_range=c(-lim,-lim/2,0,lim/2,lim),value_colors=colours,row_anno_df=ann,row_anno_cols='Program',col_anno_df=ann,col_anno_cols='Program',palette_overrides=list(Program=pcols),ann_row_args=list(show_annotation_name=F,show_legend=F,simple_anno_size=unit(2,'mm')),ann_top_args=list(show_legend=F,simple_anno_size=unit(2,'mm')),right_annotation=marks,cluster_rows=F,cluster_columns=F,show_row_names=F,show_column_names=F,border=T,width=unit(165,'mm'),height=unit(165,'mm'),raster_resize_mat=FALSE,heatmap_legend_param=list(direction='horizontal',title_position='leftcenter',legend_width=unit(55,'mm'),at=c(-lim,0,lim),labels=c('<= -0.35','0','>= +0.35')))
 paint=function(){draw(hm,heatmap_legend_side='bottom',padding=unit(c(20,3,20,3),'mm'));decorate_heatmap_body(if(oriented)'Oriented correlation' else 'Raw correlation',{for(e in edges){grid.lines(x=c(e,e),y=c(0,1),gp=gpar(col='#FFFFFFB0',lwd=.35));grid.lines(x=c(0,1),y=c(1-e,1-e),gp=gpar(col='#FFFFFFB0',lwd=.35))}});upViewport(0);grid.text(title,x=.5,y=unit(1,'npc')-unit(4,'mm'),gp=gpar(fontsize=13,fontface='bold'));grid.text(sprintf('All %s genes | %d programs | identical gene order across sources',format(N,big.mark=','),length(programs)),x=.5,y=unit(1,'npc')-unit(10,'mm'),gp=gpar(fontsize=9));grid.text(if(oriented)'Standard pipeline view: gene signs aligned to their program direction. Red = aligned; blue = opposed.' else 'Original correlation signs, without program orientation. Red = positive; blue = negative.',x=.5,y=unit(9,'mm'),gp=gpar(fontsize=8));grid.text('White lines mark programs. No row scaling; colors saturate at +/-0.35. Display bins average <=8 adjacent genes; all genes included.',x=.5,y=unit(4,'mm'),gp=gpar(fontsize=7.5))}
 pdf(file.path(OUT,paste0(slug,'.pdf')),width=9,height=8.8,useDingbats=F);paint();dev.off()
 if(png_preview){png(file.path(OUT,paste0(slug,'.png')),width=9,height=8.8,units='in',res=180);paint();dev.off()}
 cat('Rendered',slug,'\n');invisible(NULL)
}
render(C,'01_joint_gene_correlations','Gene cross-correlations | Tahoe + CytoSig',png_preview=T)
render(C,'04_joint_gene_correlations_raw','Gene cross-correlations | Unmodified signs',oriented=F)
rm(C);gc()
C=read_npy(file.path(BASE,'artifacts/discovery_tahoe_correlation.npy'));render(C,'02_tahoe_gene_correlations','Gene cross-correlations | Tahoe');rm(C);gc()
C=read_npy(file.path(BASE,'artifacts/discovery_cytosig_correlation.npy'));render(C,'03_cytosig_gene_correlations','Gene cross-correlations | CytoSig');rm(C);gc()
status=system2('/work/.venv/bin/python',file.path(BASE,'scripts/13_bundle_cross_correlations.py'));stopifnot(status==0)
