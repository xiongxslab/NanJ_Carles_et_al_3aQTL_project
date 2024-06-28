library("rtracklayer")
library(dplyr)
library(ggplot2)
library(stringr)
library(gridExtra)
cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')
df.gAPA<- read.csv('../../../06.3aQTL/3aQTL/all_cts.apaQTL.signif.ps.csv', header=T) %>%
  mutate(gene=str_split_i(phenotype_id, '\\|', 2))%>%
  dplyr::select(celltype, phenotype_id, gene) %>%
  unique()

df.info<- data.frame()
for(ct in cts.list){
  x<- read.csv(paste0('../../../05APAmatrix/02.genetic_matrix/07',ct,'_imputation_MI_3.csv'),header=T, check.names = F)[,-1] %>% ## read 3' UTR and predicted proximal pA site positions
    dplyr:: mutate(strand=str_split_i(Gene, '\\|', 4),
                   Distal_APA= ifelse(strand=='+', str_split_i(Loci, '\\-',2), str_split_i(str_split_i(Loci, '\\-',1),'\\:',2)),
                   UTR_start= ifelse(strand=='+',  str_split_i(str_split_i(Loci, '\\-',1) ,'\\:',2),str_split_i(Loci, '\\-',2)), celltype=ct) %>%
    dplyr::select(celltype,Gene, strand, Predicted_Proximal_APA, Distal_APA, UTR_start)
  df.info<- rbind(df.info,x)
}

df.info<- df.info %>%
  mutate(If_gAPA= ifelse(paste(celltype, Gene) %in% paste(df.gAPA$celltype, df.gAPA$phenotype_id), 'Yes', 'No'),
         distance_pA=abs(as.numeric(Distal_APA) - Predicted_Proximal_APA),
         distance_UTR=abs(as.numeric(UTR_start) - as.numeric(Distal_APA)))
df.info$geneName<- str_split_i(df.info$Gene,'\\|',2)

df.signif.merge<- data.frame()
for(n in c('distance_pA','distance_UTR')){
  for(ct in cts.list){
  x.sig<- wilcox.test(get(n) ~ If_gAPA, data=x)
  x.signif <- data.frame(celltype=ct, groups=n, pvalue=x.sig$p.value)
  df.signif.merge<- rbind(df.signif.merge, x.signif)
  }
  ggplot(data=df.info, aes(x=celltype, y=get(n), fill=If_gAPA))+
    geom_boxplot( outlier.shape =NA)+
    theme_classic()+
    scale_fill_manual(values=c('#fdbf6f','#fb9a99'))+
    labs(y=n)+scale_y_log10()
  ggsave(file=paste0('04.all_cts.',n, '.distribution.pdf'), height=3, width=3.5)
}

write.csv(df.signif.merge, file=paste0('04.all_cts.',n,'.signif_test.csv'))