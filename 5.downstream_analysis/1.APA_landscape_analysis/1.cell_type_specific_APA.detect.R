#!/bin/R
library(data.table)
library(ggplot2)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(dplyr)

####extract the APA events more than in four cell-types
broadly_APA_extract<- function(){
  ##summarize cell type-sharing across cell types
  df.g<- as.data.frame()
  for(ct in cts.list){
    mat1<-fread(paste0('../../07',ct,"_imputation_MI_1.csv"),sep=",",header=T, check.names=F)
    gene<-data.frame(celltype=ct,gene=mat1$Gene)
    df.g<- rbind(gene, df.g)
  }
  gene.frq<- as.matrix(table(df.g$gene))
  gene.useful<-names(gene.fre[gene.fre[,1]>=4, ])
  gene.sum.useful<- subset(df.g, gene%in%gene.useful )
  
  return(gene.sum.useful)
}


## PDUI comparison between celltypes 
comparison_APAs<- function(gene.sum.useful,ct){
  gene.list=unlist(subset(gene.sum.useful, celltype==ct, select=gene))
  mat1<-fread(paste0('../../07',ct,"_imputation_MI_1.csv"),sep=",",header=T, check.names=F)
  write.table(data.frame("Genes","P_wilcox.test","P_t.test","log2FC(test/other)"),file=paste0('02.',ct,'_other_celltype.sum'),
              sep="\t",col.names=F, row.names=F,quote=F)
  for (gene in gene.list){
    ct1.list=unlist(subset(gene.sum.useful, V1!=ct & V2==gene,select=V1))
    gene.mat.test<- as.numeric(unlist(subset(mat1, Gene == gene)[,-c(1:5)]))
    for (ct1 in ct1.list){
      mat<- fread(paste0('../../07',ct1,"_imputation_MI_1.csv"),sep=",",header=T, check.names=F)
      write.table(t(subset(mat, Gene == gene)[,-c(1:5)]), file="tmp.txt",append=T,
                  col.names=F, row.names=F)
    }
    mat.other<- fread("tmp.txt",header=F)
    gene.mat.other<- as.numeric(unlist(mat.other$V1))
    sum.test<- cbind(gene, wilcox.test(gene.mat.test, gene.mat.other)$p.value,
                     t.test(gene.mat.test, gene.mat.other)$p.value,
                     log2(mean(gene.mat.test)/mean(gene.mat.other)))
    write.table(sum.test,file=paste0('02.',ct,'_other_celltype.sum'),
                sep="\t",append=T, col.names=F, row.names=F,quote=F)
    file.remove("tmp.txt")
  }
  
  sum_raw<- fread(paste0('02.',ct,'_other_celltype.sum'),header=T)
  P_wilcox.test.fdr<- p.adjust(as.numeric(unlist(sum_raw$P_wilcox.test)),"BH")
  P_wlicox.test.bonf<- p.adjust(as.numeric(unlist(sum_raw$P_wilcox.test)),"bonferroni")
  P_t.test.fdr<- p.adjust(as.numeric(unlist(sum_raw$P_t.test)),"BH")
  P_t.test.bonf<- p.adjust(as.numeric(unlist(sum_raw$P_t.test)),"bonferroni")
  p.sum<- data.frame(sum_raw, P_wilcox.test.fdr,  P_wlicox.test.bonf, P_t.test.fdr,P_t.test.bonf )%>%
    select('Genes','P_wilcox.test.fdr','log2FC.test.other.')%>%
    filter(P_wilcox.test.fdr<=0.05)
  p.sig.short<- data.frame(ct,"Shortening",
                           unique(str_split_i(unlist(p.sum.sig[which(p.sum.sig$log2FC.test.other.<0),1]), '\\|',2)))
  p.sig.long<- data.frame(ct,"Lengthening",
                          unique(str_split_i(unlist(p.sum.sig[which(p.sum.sig$log2FC.test.other.>0),1]), '\\|',2)))
  write.table(p.sig.short, file='06.all_cell_types.sig_gene.csv',
              sep=",",append=T, col.names=F, row.names=F, quote=F)
  write.table(p.sig.long, file='06.all_cell_types.sig_gene.csv',
              sep=",",append=T, col.names=F, row.names=F, quote=F)
  p.gene<- cbind(ct, unique(str_split_i(unlist(p.sum[,1]),'\\|',2)))
  write.table(p.gene, file='07.all_cell_types.all_gene.csv',
              sep=",",append=T, col.names=F,row.names=F,quote=F)
}

###cell-type-specific APA identification
specific_extract<- function(pdui.gen, tp){
  ####filter the gene in more than 2 cell types specific
  pdui.num<- subset(pdui.gen, type==tp) %>% count(transcript)
  pdui.filter<- subset(pdui.gen, type==tp& transcript %in% unlist(pdui.num[pdui.num$n==1,]$transcript))
  ### define the shared APA in Exc and Inh as neuron-sepecific
  for(gen in unique(subset(pdui.num, n==2)$transcript)){
    
    if(subset(pdui.gen, transcript ==gen &type ==tp)$celltype[1] %in% c('Exc','Inh') &subset(pdui.gen, transcript ==gen&type ==tp)$celltype[2] %in% c('Exc','Inh')){
      pdui.filter1<- subset(pdui.gen, transcript ==gen&type ==tp)
      pdui.filter1$celltype <-'Neuron'
      pdui.filter<- rbind(pdui.filter, pdui.filter1)
    }
  }pdui.gen
  pdui.filter<-as.data.frame(unique(pdui.filter[order(factor(pdui.filter$celltype, levels=cts.list)), ] ) )
  row.names(pdui.filter)<- paste0(pdui.filter$transcript,'_', pdui.filter$celltype)
  if(!(file.exists('all_cts.specific.res'))){
    write.table(pdui.filter, file='all_cts.specific.res', row.names=F, sep='\t',quote=F)
  }else{
    write.table(pdui.filter, file='all_cts.specific.res', row.names=F,col.names=F, sep='\t',quote=F, append=T)
  }
  
}

##########Main function###################
cts.list = c('Ast','Exc','Inh','Mic','Oli','Opc')
gene.sum.useful<- broadly_APA_extract()
for(ct in cts.list){
  comparison_APAs(gene.sum.useful,ct)
}
pdui.gen<- fread('../specific_gene_plot/01.celltypes_specific.pdui_mean.res',header=F,
                 col.names= c("transcript",'celltype','type',cts.list[-1])) %>%
  na.omit()%>%
  as.data.frame()
for(tp in c('shortening','lengthening')){
  specific_extract(pdui.gen, tp)
}



