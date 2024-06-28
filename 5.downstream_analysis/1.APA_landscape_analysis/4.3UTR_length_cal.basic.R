#!/bin/bash
library(data.table)
library(dplyr)
library(stringr)

UTR_length_cal<-function(samp){
  if(0){   # run once
    # read in APA annotation
    anno_loc = read.table('All.pasites.genome_loc.tsv',head=F)
    anno = read.table('All.pasites.annotation.tsv',head=F)
    anno_loc$idx = paste0(anno_loc$V1,":",anno_loc$V2)
    anno$idx = paste0(anno$V6,":",anno$V2)
    anno = merge(anno_loc,anno,by="idx",suffix=c(".gppy",".scape"))
    # filter for 3'UTR 
    anno = anno[anno$V10 %in% c('3UTRs','3UTRs_1kb','3UTRs_2kb','LastExon1Kb') & anno$V4.gppy %in% c('downstream','exon'),]
    
    tx_info = read.table('hg38.10x.genes.info.txt',head=T,sep="\t")
    # consider coding-genes only
    tx_info = tx_info[tx_info$gene_type == "protein_coding",]
    tx_info.cut = tx_info[,c('tx_name','utr5_len','cds_len','utr3_len')]
    
    anno.info = merge(anno,tx_info.cut,by.x='V1.gppy',by.y='tx_name')
    anno.info.utr = anno.info[anno.info$V4.gppy == "exon",]
    anno.info.dstr = anno.info[anno.info$V4.gppy == "downstream",]
    anno.info.utr$pA_UTR3_len = anno.info.utr$V3.gppy - anno.info.utr$utr5_len - anno.info.utr$cds_len
    anno.info.dstr$pA_UTR3_len = anno.info.dstr$utr3_len + anno.info.dstr$V3.gppy
    anno.info = rbind(anno.info.utr,anno.info.dstr)
    anno.info.cut = anno.info[,c('V1.gppy','idx','V3.gppy','V4.gppy','V12','V8','utr5_len','cds_len','utr3_len','pA_UTR3_len')]
    colnames(anno.info.cut)[1:6] = c('TxName','Tx.pa_idx','Tx.pos','pa_Type','Genome.pa_idx','GeneName')
    rownames(anno.info.cut) = anno.info.cut$Genome.pa_idx
    write.table(anno.info.cut,'All.pasites.annotation.pa_UTR_length.tsv',quote=F,sep="\t",row.names=F)
  }
  anno.info.cut = read.table('All.pasites.annotation.pa_UTR_length.tsv',head=T)
  rownames(anno.info.cut) = anno.info.cut$Genome.pa_idx
  pa = read.table(paste0(samp,'.APA/pasite.csv.gz'),head=T,sep=",")
  pa.name = data.frame(str_split_fixed(pa[,1],':',4))
  rownames(pa) = paste0(pa.name[,1],':',pa.name[,2],':',pa.name[,4])
  pa = pa[,-1]
  pa.name$idx = paste0(pa.name[,1],':',pa.name[,2],':',pa.name[,4])
  pa.name = pa.name[pa.name$idx %in% anno.info.cut$Genome.pa_idx,]
  anno.info.cut.samp = anno.info.cut[as.character(pa.name$idx),]
  pa = pa[as.character(anno.info.cut.samp$Genome.pa_idx),] 
  lengths = anno.info.cut.samp$pA_UTR3_len
  # calculate weighted length for each gene in each cell
  pa2 = apply(pa, 2, function(x) x*lengths)  
  pa.bygene.total = aggregate(pa2, by=list(as.character(anno.info.cut.samp$TxName)), FUN=sum)
  pa.bygene.npa = aggregate(pa, by=list(as.character(anno.info.cut.samp$TxName)), FUN=sum)	
  all(rownames(pa.bygene.total) == rownames(pa.bygene.npa)) # check order
  all(colnames(pa.bygene.total) == colnames(pa.bygene.npa))
  pa.bygene.total.mat = as.matrix(pa.bygene.total[,-1])
  pa.bygene.npa.mat = as.matrix(pa.bygene.npa[,-1])
  pa.bygene.weighted.len = pa.bygene.total.mat/pa.bygene.npa.mat
  rownames(pa.bygene.weighted.len) = pa.bygene.total[,1]
  saveRDS(pa.bygene.weighted.len,paste0(samp,'.pa_UTR3_avg_length.bygene.rds'))
}


mean_rm_na<- function(x){
  s<- mean(na.omit(x))
}

### calculate mean 3 'UTR length 
mean_cal_cell<- function(id){
  df<- readRDS(paste0( id,'.pa_UTR3_avg_length.bygene.rds'))%>%
    as.data.frame()
  ## across gene
  df.sum<- data.frame(index=paste(id, str_split_i(names(df), '\\.',1),sep=':'),
                      length_cell=apply(df,2, mean_rm_na) )
  df.sum<- merge(df.info1%>% select(cellid, major.celltype), df.sum,by.x='cellid', by.y='index') 
  
  df.sum1<- data.frame()
  for(ct in unique(df.sum$major.celltype)){
    x<- data.frame(ID=id, major.celltype= ct, length_cell=mean(df.sum[which(df.sum$major.celltype==ct),]$length_cell) )
    df.sum1<- rbind(df.sum1,x)
  }
  
  ## across cell
  df.t<- t(df) %>% as.data.frame()%>% mutate(index= paste(id,str_split_i(names(df),'\\.',1),sep=':'))
  df.t<- merge(df.info%>% select(cellid, major.celltype), df.t, by.x='cellid', by.y='index')
  
  df.sum_g<- data.frame()
  for(ct in unique(df.t$major.celltype)){
    df.ct<- subset(df.t, major.celltype== ct)
    df.ct.sum<- apply(df.ct[,-c(1:2)], 2, mean_rm_na) 
    df.ct.sum<- data.frame(ID=id, major.celltype= ct,length_gene=mean(na.omit(df.ct.sum)))
    df.sum_g<-rbind(df.sum_g, df.ct.sum)
  }
  
  df.sum_all<- merge(df.sum1, df.sum_g)
  if(!(file.exists('10.3UTR_length/01.celltypes_3UTR_length.res'))){
    write.table(df.sum_all, file='10.3UTR_length/01.celltypes_3UTR_length.res',sep='\t',
                row.names=F, quote=F)
  }else{
    write.table(df.sum_all, file='10.3UTR_length/01.celltypes_3UTR_length.res',sep='\t',
                row.names=F, quote=F, col.names = F, append=T)
  }
  
}

args<- commandArgs(T)
id=args[1]
UTR_length_cal(id)
mean_cal_cell(id)
