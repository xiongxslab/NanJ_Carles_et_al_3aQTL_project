library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')
cts.cols = c("#33A02C", "#B2DF8A", "#E31A1C", "#FDBF6F", "#B15928", "#CAB2D6")

est <- fread("03.all_cts.cells.filter.all.est",  header = FALSE, col.names = c("celltype","cell","protein", "logOR", "lci", "uci"))
est <- est %>% mutate(se = (uci - lci)/(2*1.96), zval = logOR/se, pval = 2*(1 - pnorm(abs(zval))),
                      padj = p.adjust(pval, method = "bonferroni", n = length(pval)))
est$celltype<- factor(est$celltype, levels=cts.list)
ggplot(data=est, aes(x=logOR, y=-log10(padj), col=celltype))+
  geom_point()+
  scale_color_manual(values=cts.cols)+
  theme_classic()
ggsave(file='01.all_cts.logOR_bonf.pdf', height=5, width=4)

write.table(est, file = "./04.all_cts.cells.filter.est.sum", quote = FALSE, sep = "\t", row.names = FALSE)

est<- est%>%
  mutate( type=str_split_i(protein, '\\.',2),protein=str_split_i(protein, '\\.',1))

### update the cor. between expression and 3'UTR length

g.r<- c('ELAVL1','CDK12','HNRNPC','SRSF3','SRSF7','CPEB4','SCAF4','SCAF8','DICER1','CPEB1',
        'MBNL1', 'MBNL2','RBM3', 'CIRBP', 'NOVA', 'NFKB','CPSF6','CPSF5','CPEB2', paste0('CPSF',seq(7)),'WDR33',
        'FIP1','NUDT21','FIP1L1',paste0('CSTF',seq(3)),'CSTF2T','SYMPK','PAF1','PABPN1','RBBP6','SCAF4','SCAF8') ##(Mitschka S, Mayr C. 2022 Nature Reviews Molecular Cell Biology)
df.merge<-data.frame()
for (ct in cts.list){
  df<- fread(paste0('01.',ct,'.expr_len.cor.res'), header=T) %>%
    mutate(bonf=p.adjust(pval))
  df.merge<- rbind(df.merge, df)
}

est1<- merge(est, df.merge, by.x=c('celltype','protein'), by.y=c('celltype','gene')) %>%
  na.omit()
est1$celltype<- factor(est1$celltype, levels=cts.list)

est1$group<- '+.n.'
est1[which(est1$logOR<0 & est1$padj<0.05),]$group<-'-.s.'
est1[which(est1$logOR>0 & est1$padj<0.05),]$group<-'+.s.'
est1[which(est1$logOR>0 & est1$padj>0.05),]$group<-'+.n.'
est1$group<- factor(est1$group, levels=c('+.s.','-.s.','+.n.','-.n.'))

##update whether reported
est1<- est1%>%
  mutate(ref=ifelse(protein %in% g.r, 'reported','novel'))
write.table(est1, file='05.all_cts.torus.cor.res',append=F, sep='\t', quote=F, row.names = F)


ggplot(data=est1, aes(x=logOR, y=abs(cor), col=celltype))+
  geom_point()+
  scale_color_manual(values=cts.cols)+
  stat_cor()+
  geom_vline(xintercept = 1, linetype='dashed')+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_classic()

df.sum<- est1%>% filter(logOR>0) %>%
  select(celltype,group) %>%
  table() %>%
  as.data.frame() %>%
  filter( group %in% c('+.s.','+.n.'))

ggplot(data=est1%>% filter(logOR>0), aes(x=celltype, y=-log10(bonf), fill=group))+
  geom_boxplot()+
  geom_text(data=df.sum,aes(x=celltype, color=group, y=0, label=Freq))+
  #scale_fill_manual(values=cts.cols)+
  # geom_vline(xintercept = 1, linetype='dashed')+
  #  geom_hline(yintercept = 0,linetype='dashed')+
  theme_classic()
ggsave(file='05.all.cts.torus.bonf.pdf', height=5, width=6)

