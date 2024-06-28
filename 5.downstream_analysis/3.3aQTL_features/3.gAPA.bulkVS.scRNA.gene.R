library(data.table)
library(ggplot2)
library(stringr)

cts.list=c("Exc","Inh","Ast","Oli","Opc","Mic")
g.bulk_all<- fread('2021.NG_all_tissues_bulk.apaGene.csv', header=T, sep=',') ## read gAPAs in all tissues from (Lei Li et al. 2019 Nature Genetics)
g.list.bulk_all<- unique(str_split_i(g.bulk_all$`APA genes`, '\\|',2))

g.list.bulk_sig<- fread("2021.NG.cortex.bulk.apaGene", header=F) ## read gAPAs in brain cortex from (Lei Li et al. 2019 Nature Genetics)
g<- fread('231106.cortext_sig.not_anno.UCSC.gene',header=F) ## liftover genes
g.list.bulk_sig<- unique(rbind(g$V2,subset(g.list.bulk_sig,!(V1 %in% g$V1)),use.names=FALSE))
names(g.list.bulk_sig)<-'V1'
g.list.bulk_all<- c(g.list.bulk_all[!(g.list.bulk_all%in% g$V1)],g$V2)
g.sc.sig<-fread('01.all_cts.perm_gAPA.filter.csv', header=T,sep=',')%>% dplyr::select(celltype, gene)%>% unique()

df.com<- data.frame()
g.apa_test<- c()
for(ct in cts.list){
  g.ct<- fread(paste0('00.',ct,'.apaGene.tested.gene'), header=F) ## all APAs used for 3'aQTL calling
  g.apa_test<- unique(c(g.apa_test, g.ct$V1))
  g.ct.sig<- subset(g.sc.sig, celltype==ct)$gene
  ##shared between bulk and scRNA-Seq
  g.share<- intersect (g.ct.sig, g.list.bulk_sig$V1)
  df.com<- rbind(df.com,data.frame(celltype=ct, group='shared', gene=g.share))
  ## scRNA-seq specific
  g.sc<- g.ct.sig[!(g.ct.sig %in% g.share)& g.ct.sig %in% g.list.bulk_all]
  df.com<- rbind(df.com,data.frame(celltype=ct, group='scRNA-specific', gene=g.sc))
  ## bulk-RNA-seq specific
  g.bulk<- g.list.bulk_sig$V1[!(g.list.bulk_sig$V1 %in% g.share)& g.list.bulk_sig$V1 %in% g.ct$V1]
  df.com<- rbind(df.com,data.frame(celltype=ct, group='bulk-RNA-specific', gene=g.bulk))
}

df.com$group<- factor(df.com$group, levels=c('scRNA-specific','shared','bulk-RNA-specific'))
df.com$celltype<- factor(df.com$celltype, levels=cts.list)
write.table(df.com, file='6.scRNAvs.bulk/01.sc_bulk_compare.res',sep='\t', quote=F, col.names=T, row.names=F)

ggplot(data=df.com, aes(x=celltype, fill=group))+
  geom_bar(stat='count', position = 'dodge',width = 0.8)+
  scale_fill_manual(values=c("#e4967d", "#bc6645",'#418f88'))+
  theme_classic()+
  geom_text(stat='count',aes(label=..count..),
            position = position_dodge(width = 0.8), 
            vjust=-0.3,size=3.5)
ggsave(file='06.scRNA_bulkRNA.gene_compare.pdf', height=3, width=4.5)  

###compare apaGene from scRNA of all celltypes and bulk 
g.test<- intersect(g.apa_test, g.list.bulk_all)
g.bulk<- intersect(g.list.bulk_sig$V1, g.test)
g.scRNA<- intersect(g.sc.sig$gene, g.test)

library(VennDiagram)
# Create the Venn diagram
venn.diagram(
  x = list(scRNA_seq = g.scRNA, bulk_RNA = g.bulk ),imagetype = 'svg',
  filename = "06.scRNA_bulk.venn_diagram.svg",area.vector=F,
  cex = 1, fontfamily = 'serif', fill = c("#CD853F", '#556B2F'),
  cat.cex = 1, cat.fontfamily = 'serif',margin=0.2
)

