#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(digest);library(limma);
  library(AnnotationDbi);library(org.Hs.eg.db)})
root <- '/work/results/program_pilot/130_HGF_source_platform_audit/linear'
dir.create(file.path(root,'data/source'),recursive=TRUE,showWarnings=FALSE)
native <- list()
for(study in c('GSE26732','GSE26735')) {
  name<-paste0(study,'_series_matrix.txt.gz')
  path<-file.path(root,'data/source',name)
  if(!file.exists(path))download.file(paste0('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/',study,'/matrix/',name),path,quiet=TRUE,mode='wb')
  connection<-gzfile(path,'rt');lines<-readLines(connection,warn=FALSE);close(connection)
  metadata<-lines[grepl('^!Sample_(title|geo_accession|characteristics_ch1|treatment_protocol_ch1|data_processing|platform_id)',lines)]
  writeLines(metadata,file.path(root,'data',paste0(study,'_sample_metadata.txt')))
  cat(study,'\n',paste(metadata,collapse='\n'),'\n')
  begin<-grep('^!series_matrix_table_begin',lines)
  end<-grep('^!series_matrix_table_end',lines)
  table<-fread(text=paste(lines[(begin+1):(end-1)],collapse='\n'))
  native[[study]]<-table
}
definition<-qread('/work/results/program_pilot/121_HGF_consensus_coordinate_bridge/linear/data/HGF_consensus_coordinate_definition.qs')
meta<-definition$metadata[study %in% c('GSE26732','GSE26735')]
cy<-fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',meta$profile_id))
setnames(cy,1,'gene')
columns<-meta$profile_id
Z<-as.matrix(cy[,..columns])
fwrite(data.table(genes=sum(complete.cases(Z)),between_platform_gene_logFC_r=cor(Z[,1],Z[,2],use='complete.obs'),
  Affymetrix_gene_logFC_SD=sd(Z[,which(grepl('GSE26732',columns))],na.rm=TRUE),
  Illumina_gene_logFC_SD=sd(Z[,which(grepl('GSE26735',columns))],na.rm=TRUE)),
  file.path(root,'data/cytosig_platform_agreement.tsv'),sep='\t')

# Reconstruct exactly declared native24h-HGF versus untreated-WT contrasts.
# Affymetrix is already RMA log2. Illumina release is linear, includes negative
# background-corrected values: compare the transparent floor-log convention
# and standard normexp+quantile preprocessing as a molecular sensitivity check.
responses<-replicate_audit<-list()
for(study in names(native)) {
  gpl<-if(study=='GSE26732')'GPL96' else 'GPL2725'
  gpl_path<-file.path(root,'data/source',paste0(gpl,'.annot.gz'))
  if(!file.exists(gpl_path))download.file(paste0('https://ftp.ncbi.nlm.nih.gov/geo/platforms/',
    if(gpl=='GPL96')'GPLnnn' else 'GPL2nnn','/',gpl,'/annot/',gpl,'.annot.gz'),gpl_path,quiet=TRUE,mode='wb')
  lines<-readLines(gzfile(gpl_path),warn=FALSE)
  begin<-grep('^!platform_table_begin',lines,ignore.case=TRUE);end<-grep('^!platform_table_end',lines,ignore.case=TRUE)
  annotation<-fread(text=paste(lines[(begin+1):(end-1)],collapse='\n'))
  entrez<-sub(' .*','',sub('///.*','',annotation[['Gene ID']]))
  valid<-unique(entrez[grepl('^[0-9]+$',entrez)])
  symbols<-unname(mapIds(org.Hs.eg.db,keys=valid,keytype='ENTREZID',column='SYMBOL',multiVals='first'))
  annotation[,gene:=symbols[match(entrez,valid)]]
  table<-native[[study]]
  E<-as.matrix(table[,-1])
  profiles<-if(study=='GSE26732')list(native_RMA_log2=E) else {
    corrected<-backgroundCorrect(new('EListRaw',list(E=E)),method='normexp',offset=16)$E
    list(released_signal_floor_log2=log2(pmax(E,0)+1),
      normexp_quantile_log2=normalizeBetweenArrays(log2(corrected),method='quantile'))
  }
  for(method in names(profiles)) {
    X<-profiles[[method]]
    gene<-annotation$gene[match(table$ID_REF,annotation$ID)]
    df<-data.table(gene,logFC=rowMeans(X[,9:10])-rowMeans(X[,1:2]),
      WT_difference=X[,1]-X[,2],HGF_difference=X[,9]-X[,10],
      replicate1_logFC=X[,9]-X[,1],replicate2_logFC=X[,10]-X[,2],
      control_mean=rowMeans(X[,1:2]))[!is.na(gene)&gene!='']
    df<-df[,lapply(.SD,mean),by=gene]
    index<-match(df$gene,cy$gene)
    study_name<-study
    cy_column<-meta[study==study_name,profile_id]
    stopifnot(length(cy_column)==1L)
    cy_fc<-cy[[cy_column]][index]
    responses[[paste(study,method)]]<-df[,`:=`(study=study,preprocessing=method,cached_CytoSig_logFC=cy_fc)]
    # Whole-gene moments are descriptives, not a gene-selection filter.
    replicate_audit[[paste(study,method)]]<-data.table(study,preprocessing=method,genes=nrow(df),
      response_SD=sd(df$logFC),WT_difference_SD=sd(df$WT_difference),HGF_difference_SD=sd(df$HGF_difference),
      expected_mean_contrast_noise_SD=sqrt((var(df$WT_difference)+var(df$HGF_difference))/4),
      replicate_response_r=cor(df$replicate1_logFC,df$replicate2_logFC),
      cached_CytoSig_r=cor(df$logFC,cy_fc,use='complete.obs'))
  }
}
fwrite(rbindlist(responses),file.path(root,'data/reconstructed_platform_gene_responses.tsv'),sep='\t')
fwrite(rbindlist(replicate_audit),file.path(root,'data/platform_replicate_reliability.tsv'),sep='\t')
print(rbindlist(replicate_audit))
