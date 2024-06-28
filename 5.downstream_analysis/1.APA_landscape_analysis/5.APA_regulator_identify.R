library(data.table)
library(dplyr)
library(ggplot2)


cor_cal<- function(df, gr){
  ### RBP and TF datasets were downloaded from ChIP-seq and CLIP-seq of ENCODE
  meta.tf<- fread('ENCODE/metadata.tsv', header=T) %>%
    select(Assay,`Experiment target`) %>%
    mutate(type=ifelse(Assay=='TF ChIP-seq', 'TF', 'RBP'),
           protein=str_split_i(`Experiment target`, '-',1))
  meta.rbp<- fread('RNA_binding_anntotation/human.txt.gz', header=F) %>%
    select(V6)
  
  df.meta<- merge(df, meta, by='cellid')
  df.meta<- merge(df.meta, df.id) %>%
    group_by(projid,major.celltype) %>%
    summarise(mean.length=mean(mean.len))
  
  Merge.cor<- data.frame()
  Merge.e<- data.frame()
  for(ct in cts.list){
    ##read expression matrix
    df.e<- fread(paste0(ct,'_expression_QC.voom.sort.bed.gz'),header=T) %>%
      filter(`#Chr` %in% seq(22))%>%
      t()%>%
      as.data.frame()
    names(df.e)<- df.e[row.names(df.e)=='Symbol',]
    df.e<- df.e[-(1:4),] 
    df.e$projid<- rownames(df.e)
    df.e1<- data.frame(celltype=ct, df.e)
    if(nrow(Merge.e)==0){
      Merge.e<- df.e1
    }else{
      Merge.e<- rbind(subset(Merge.e,select=intersect(names(Merge.e), names(df.e1))), subset(df.e1, select=intersect(names(Merge.e), names(df.e1))))
    }
    df.meta.e<- merge(df.meta[which(df.meta$major.celltype==ct),], df.e, by='projid')
    df.meta.e[,-c(1:2)]<- apply(df.meta.e[,-c(1:2)], 2, as.numeric)
    for(i in seq(4, ncol(df.meta.e))){
      a<- cor.test(df.meta.e[,4],df.meta.e[,i])
      df.cor<- data.frame(celltype=ct, gene=names(df.meta.e)[i], cor=a$estimate, pvalue=a$p.value)
      Merge.cor<- rbind(Merge.cor,df.cor)
    }
  }

  Merge.signif<- merge.cor %>%
    mutate(p.adj= p.adjust(pvalue, method='bonferroni'))
    filter(p.adj < 0.01) %>%
      mutate(type=ifelse(!(gene %in% c(unique(meta.tf$protein),meta.rbp$V6)), 'other',' ' ),
      type= ifelse(gene %in% meta.tf[which(type=='TF'),]$protein & type !='other', 'TF','RBP'),
      type= ifelse(gene %in% intersect(meta.tf[which(type=='TF'),]$protein, c(meta.tf[which(type=='RBP'),]$protein, meta.rbp$V6)),  'TF&RBP',type),
      group=ifelse(gene %in% g.report, 'reported','novel'))
  saveRDS(Merge.signif, file=paste0('all_cts.length_expr.signif.cor.rds'))
}

##read meta information of snRNA-seq data
meta = read.table('consensus_annot_snPFC.joint_annotation_metadata3.tsv.gz',head=T,sep="\t")
df.id<- fread('AD430_fastq_projid_mapping.tsv', header=T)%>%
  select('projid','library_id')
df<- readRDS('cell.mean.len.perc.cut0.05.rds') 
cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')
cts.cols = c("#33A02C", "#B2DF8A","#E31A1C", "#FDBF6F", "#B15928", "#CAB2D6")
##known regulators were reported by Sibylle Mitschka and Christine Mayr Nat Rev Mol Cell Biol. 2021)
g.report<- c('ELAVL1','CDK12','HNRNPC','SRSF3','SRSF7','CPEB4','SCAF4','SCAF8','DICER1','CPEB1',
             'MBNL1', 'MBNL2','RBM3', 'CIRBP', 'NOVA', 'NFKB','CPSF6','CPSF5','CPEB2', paste0('CPSF',seq(7)),'WDR33',
             'FIP1','NUDT21','FIP1L1',paste0('CSTF',seq(3)),'CSTF2T','SYMPK','PAF1','PABPN1','RBBP6','SCAF4','SCAF8')
df$cellid<- rownames(df)
cor_cal(df, 'cut0.05')


