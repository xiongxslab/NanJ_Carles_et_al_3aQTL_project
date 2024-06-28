library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggplotify)

ct.list = c('Ast','Exc','Inh','Mic','Oli','Opc')

# read in the permutation QTLs
df.perm = data.frame()

for(ct in ct.list){     
  df = read.table(paste0(ct,'.permute.out.qvalue.signif.txt'),head=T)
  df = df[,c('phenotype_id','variant_id','distance','pval_nominal','slope')]
  df$celltype1 = ct
  df.perm = rbind(df.perm,df)
}
df.perm$idx = paste0(df.perm$phenotype_id,':',df.perm$variant_id)

# read in all QTLs
df.all = data.frame()
for(ct in ct.list){
  print(ct)
  df = read.table(paste0(ct,'.apaQTL.sumstat.txt.gz'),head=F)
  df = df[,c(1,2,7,8)]
  colnames(df) = c('phenotype_id','variant_id','pval_nominal','slope')
  df$celltype2 = ct
  df$idx = paste0(df$phenotype_id,':',df$variant_id)
  df = df[df$idx %in% as.character(df.perm$idx),]
  df.all = rbind(df.all,df)
}

df.compare = merge(df.perm,df.all,by='idx',suffix=c('.discovery','.replication'))

df.compare$idx2 = paste0(df.compare$celltype2,":",df.compare$idx)
df.compare$is.signif.replication = 'N'
# annotate the significant QTLs in tissue 2
for(ct in ct.list){
  print(ct)
  df = read.table(paste0('../',ct,'.apaQTL.signif.txt'),head=T)
  df$idx2 = paste0(ct,':',df$phenotype_id,':',df$variant_id)
  df.compare$is.signif.replication[df.compare$idx2 %in% as.character(df$idx2)] = 'Y'
}

# Merge = readRDS(paste0("Cross.ct.compare.apaQTL.rds"))
df.compare$celltype1 = factor(df.compare$celltype1,levels = cts.list)
df.compare$celltype2 = factor(df.compare$celltype2,levels = cts.list)

df.compare = df.compare[order(df.compare$is.signif.replication),]
p<-ggplot(df.compare,aes(x= slope.discovery,y=slope.replication,colour= is.signif.replication,shape=is.signif.replication)) + geom_point(size=0.6) +  
  theme(aspect.ratio = 1) + scale_colour_manual(values=c("grey70", "brown3"),name='Sgnificance in Cell type2',breaks=c('N','Y'),labels=c('Insignificant','Significant')) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + facet_grid(celltype1~celltype2)+
  labs(x='Slope in Cell type1', y='Slop in Cell type2')+
  scale_shape_discrete(name='Sgnificance in Cell type2',breaks=c('N','Y'),labels=c('Insignificant','Significant'))+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

ggsave('07.apaQTL.across_cts.pdf', height=6,width=8)


###calculate the pearson correlation across cell types
cor.mat<- as.data.frame(matrix(0, 6,6))
names(cor.mat)<- cts.list
row.names(cor.mat)<- cts.list
meth='pearson'
for(ct1 in cts.list){
  for(ct2 in cts.list){
    merge.ct<- subset(Merge, celltype1==ct1 & celltype2==ct2|celltype1==ct2 & celltype2==ct1 , select=c('variant_id.discovery','variant_id.replication','slope.discovery','slope.replic
ation'))
    ### filter the duplicated pairs
    merge.ct.uniq<- unique(merge.ct)
    ### calculate the pearson correlation
    attach(merge.ct.uniq)
    cor.mat[rownames(cor.mat)==ct1, names(cor.mat)==ct2]= unlist(cor.test(slope.discovery,slope.replication, method = meth)$estimate)
    detach(merge.ct.uniq)
  }
}

p<- as.ggplot(pheatmap( cor.mat,cluster_rows = T,cluster_cols=T,fontsize=18,number_color='black',color=colorRampPalette(c("#EEDFCC","brown3"))( 10 ),
                        display_numbers = T, border=F))

ggsave(paste0('07.apaQTL.across_cts.cor',meth,'.pdf'),height=4,width=5)
