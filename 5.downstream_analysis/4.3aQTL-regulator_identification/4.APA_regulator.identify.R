library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

cts.list=c("Exc", "Inh", "Ast", "Oli", "Opc", "Mic")
cts.cols=c("#33A02C", "#B2DF8A", "#E31A1C", "#FDBF6F", "#B15928", "#CAB2D6")
### extract the regulators significantly enriched by apaQTL
# df<- fread('04.all_cts.cells.filter.est.sum', header=T) %>%
#   mutate(fdr=p.adjust(pval, method='fdr')) %>%
#   filter(logOR>0, fdr<0.05)
# write.csv(df, file='07.all_cts.enrich_significant.csv')

####exatact the regulator expression and target gene expression
target_extract<- function(ct, cellline,type){
  df1<- df%>% 
    filter(celltype==ct, cell==cellline)
  df.b<- fread(paste0('01.', ct,'.snpanno.', cellline,'.filter100.txt.gz'), header=T) 
  df.b<- subset(df.b, select=c('SNP', names(df.b)[str_split_i(names(df.b),'\\_',1) %in% 
                                                                paste(str_split_i(df1$protein, '\\.',1),str_split_i(df1$protein, '\\.',2),sep='.')]))
  if(type == 'apaQTL'){
    df.q<- fread(paste0('/data/slurm/nanjh/scAPA_project/06.3aQTL/3aQTL/',ct,'.apaQTL.signif.txt'), header=T)
  }else if(type =='leadQTL'){
    df.q<- fread(paste0('/data/slurm/nanjh/scAPA_project/06.3aQTL/lead3aQTL/',ct,'.permute.leadVariant.txt'), header=T)
  }
  
  df.pdui<- fread(paste0('/data/slurm/xiongxs/Proj_APA/2.cisQTL/2.voom/',ct,'_imputation_MI_3_naomit.parse.voom.sort.bed.gz'), header=T)
  df.e<- fread(paste0('/data/slurm/nanjh/scAPA_project/04.expression_matrix/re-QC/02.',ct,'_expression_QC.voom.sort.bed.gz'), header=T)
  for(pro1 in intersect(unique(str_split_i(names(df.b)[-1],'\\.',1)),str_split_i(df.pdui$ID,'\\|',2))){
    print(paste0(pro1,' analysis start...'))
    pro.l<- names(df.b)[str_split_i(names(df.b),'\\.',1)%in% pro1]
    for(pro in pro.l){
      print(paste0(pro,' analysis start...'))
      ## extact the SNP binding at regulators
      x<- subset(df.b, select=c('SNP',pro)) %>%
        rename(regulator=pro) %>%
        filter(regulator==1)
      ## extact the gene binding at regulators
      y<- subset(df.q, variant_id %in% x$SNP)
      ## extract the target gene PDUI
      df.p.y<- subset(df.pdui,ID %in% y$phenotype_id) %>%
        select(-`#Chr`, -start,-end) %>%
        t()%>%
        as.data.frame() 
      names(df.p.y)<- unlist(df.p.y[1,])
      df.p.y$projid<- row.names(df.p.y)
      df.p.y<- df.p.y[-1,]
      ##extract the regulator expression
      df.e.p<- subset(df.e, Symbol==str_split_i(pro, '\\.',1)) %>%
        select(-`#Chr`, -start,-end) %>%
        t()%>%
        as.data.frame()
      names(df.e.p)<- unlist(df.e.p[1,])
      df.e.p$projid<- row.names(df.e.p)
      df.e.p<- df.e.p[-1,]
      ###merge the pdui and expression
      z<- merge(df.e.p, df.p.y, by='projid')
      
      saveRDS(z, file=paste0('local_regulator/pdui_expr/',pro,'.',ct,'.pdui_expr.RDS'))
      ### aggrate the mean pdui of all target genes and regulator expression
      z[,-1]<- apply(z[,-1],2, as.numeric)
      z.agg<- data.frame(z[1:2], PDUI=apply(z[,-c(1:2)],1,mean))
      s<- cor.test(z.agg[,2],z.agg[,3])
      res<- data.frame(celltype=ct,cellline=cellline,regulator=pro,cor=s$estimate, signif=s$p.value)
      write.table(res, file='local_regulator/01.all_cts.regulator.cor.res', sep='\t', append=T, col.names = F, row.names = F, quote=F)
      if(s$p.value<0.01){
        ggplot(data=z.agg, aes(x=z[,2], y=z[,3]))+
          geom_point(col=cts.cols[which(cts.list==ct)])+
          geom_smooth(method='lm',se=T, color='black',linetype='dashed')+
          stat_cor()+
          theme_classic()+
          labs(x=names(z.agg)[2],y=names(z.agg)[3])
        
        ggsave(file=paste0('local_regulator/cor_plot/',ct,'.',pro,".",s$estimate,'.cor.pdf'), height=3,width=3)
      }
    }
    print(paste0(pro,' analysis done...'))
  }
}

if(!(file.exists("local_regulator/pdui_expr"))){
  dir.create("local_regulator/pdui_expr")
}
if(!(file.exists("local_regulator/cor_plot"))){
  dir.create("local_regulator/cor_plot")
}

args<- commandArgs(TRUE)
ct=args[1]
cellline=args[2]
type=args[3]
df<- fread('07.all_cts.enrich_significant.csv',header=T)
target_extract(ct, cellline, type)
