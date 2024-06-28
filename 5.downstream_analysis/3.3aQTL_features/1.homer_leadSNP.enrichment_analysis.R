library(data.table)
library(dplyr)
library(ggplot2)

# cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')
df.perm<-fread('01.all_cts.perm_gAPA.filter.csv', sep=',', header=T)[,-1] ## read gAPA-leadSNP files
df.bed <- df.perm %>%
  select(CHR,pos, phenotype_id) %>%
  mutate(CHR=paste0('chr',CHR),start=pos-1, end=pos, strand= str_split_i(phenotype_id, '\\|',4))%>%
  select(CHR,start, end, strand) %>% unique
write.table(df.bed, file='5.anno_homer/01.all_cts.gAPA_leadVar.bed', quote=F, col.names = F, row.names = F, sep='\t')

command='annotatePeaks.pl 5.anno_homer/01.all_cts.gAPA_leadVar.bed hg38 -annStats 5.anno_homer/02.all_cts.gAPA_leadVar.annStats > 5.anno_homer/02.all_cts.gAPA_leadVar.anndetails'
write.table(command, file='5.homer_annoPeak_leadVar.sh', quote=F, col.names = F, row.names = F)
system('bash 5.homer_annoPeak_leadVar.sh')

### Do plot!!!
Merge<- fread('5.anno_homer/02.all_cts.gAPA_leadVar.annStats',header = T)[1:13,]
Merge1<- subset(Merge, Annotation %in% c('3UTR','TTS','TSS','Exon','Intron','Intergenic','Promoter','5UTR') )
Merge1$`LogP enrichment (+values depleted)`<- as.numeric(Merge1$`LogP enrichment (+values depleted)`)
Merge1$`Log2 Ratio (obs/exp)`<- as.numeric(Merge1$`Log2 Ratio (obs/exp)`)
Merge1 <- Merge1 %>%
  mutate(fdr= p.adjust(10^-abs(`LogP enrichment (+values depleted)`)),
         Signif= ifelse(abs(Merge1$`LogP enrichment (+values depleted)`)>2& fdr<0.05,'Y','N'))

Merge1$Annotation<- factor(Merge1$Annotation, levels=rev(c('Intergenic','3UTR','TTS','Intron','Exon','5UTR','TSS','Promoter')))
ggplot(data=Merge1, aes(x=Annotation,  y=`Log2 Ratio (obs/exp)`))+
  geom_bar(stat='identity',fill='#418f88')+
  theme_classic()+
  # scale_fill_manual(values=c('grey','#8B2323'))+
  # scale_fill_gradient2(low='#6495ED',high='#CD2626', mid='white', midpoint = 0)+
  # facet_wrap(~celltype, scales='free_y', nrow=3)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 0, linetype='dashed')
ggsave(file=paste0('05.all_cts.leadVar.enrichment_plot.pdf'), height=3, width=2.5)

write.csv(Merge1, file='05.all_cts.leadVar_enrich.csv')
