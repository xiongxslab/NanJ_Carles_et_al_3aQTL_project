#!/bin/R
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)

cts.list = c("Exc","Inh","Ast","Oli","Opc","Mic")
cts.cols = c("#33A02C", "#B2DF8A", "#E31A1C", "#FDBF6F", "#B15928", "#CAB2D6")

df.longest<- read.table('GencodeV32.20230414.hg38_3UTR_annotation.longest.bed', header=F,
                        col.names = c('CHR','start','end','phenotype_id','strand','gene'))

df.bulk <-read.table('01.Brain_Cortex.3aQTL.hg38.txt', header=T) ##bulk-3'aQTL after liftover 
snp.map<- fread(map, header=F,
                col.names=c('variant_id','CHR','position',"Reference_allele","Effect_allele"))
df.bulk<- merge(snp.map, df.bulk, by=c("CHR","position")) %>%
  mutate(beta= ifelse(Effect_allele= str_split_i(qtl.bulk.hg38.map$SNP,'\\_',3), -beta, beta)) ##flip effect size

df.sc_merge<- data.frame()
for(ct in cts.list){
  df.ct<- fread(paste0(ct,'.apaQTL.signif.txt'),header=T) %>%
    filter(phenotype_id %in% df.longest$phenotype_id) ## include the longest isoforms
  
  df.ct<- merge(snp.map, df.ct) %>%
    group_by(phenotype_id) %>%
    arrange(pval_nominal) %>%
    slice_head(n=1)%>% ## include leadSNP for each gAPA
    mutate(celltype=ct, gene= str_split_i(phenotype_id, '\\|',2))
  df.sc_merge<- rbind(df.sc_merge, df.ct)
}

df.merge<- merge(df.bulk, df.sc_merge, by=c('CHR','position','gene')) ##merge SNP-gAPA pairs


## visualization
consistency= percent(length(subset(df.merge,slope*beta >0)$variant_id)/length(df.merge$variant_id),accuracy = 0.01)
x<- cor.test(df.merge$slope, df.merge$beta)
options(scipen = 10)
options(digits = 3)
y<-format(x$p.value, scientific = TRUE)

options(digits = 2)
r<- format(x$estimate, scientific = F)
ggplot(data=df.merge,aes(x=as.numeric(slope),y=as.numeric(beta), size=-log10(pval_nominal),col=`Cell type`))+
  geom_point()+
  theme_bw()+
  labs(x='Effect size in scRNA-Seq', y='Effect size in bulk RNA-Seq')+
  geom_smooth(method='lm',se=F, color='black',linetype='dotdash')+
  scale_color_manual(values=cts.cols)+
  scale_x_continuous(expand = c(0, 0))+
  geom_vline(xintercept=0,linetype='dashed')+
  geom_hline(yintercept=0,linetype='dashed') + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  theme(legend.justification = c(0.5, 0),
        legend.position = c(0.9, 0.02), ###set the lengend position
        legend.background = element_blank(),
        legend.key = element_blank())+
  xlim(c(-0.6,0.6))+
  ylim(c(-0.6,0.6))+
  annotate("text",label=paste('Consistency',consistency, sep=': '),x=-0.45,y=0.6,size=4)+
  annotate("text",label=paste0('R=',r,', P=',y),x=-0.45,y=0.53,size=4)
ggsave(file='06.all_cts.sc_bulk.slope_compare.lead.pdf',height=4,width=4)

