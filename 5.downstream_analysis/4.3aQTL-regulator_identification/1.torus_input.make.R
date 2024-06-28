library(data.table)
library(dplyr)
library(stringr)

snpmap_make<- function(){
  x<- subset(df1, select=c('CHR','Pos','variant_id'))
  write.table(x, file=paste0('01.',ct,'.snpmap.txt'), quote=F, col.names=F, row.names=F,sep='\t')
  system(paste0('bgzip -f ','01.',ct,'.snpmap.txt'))
}

assc_make<- function(){
  x<- df1%>% 
    mutate(tstat=case_when(slope > 0 ~ qnorm(1-pval_nominal/2), slope < 0 ~ -qnorm(1-pval_nominal/2)))
  x<- subset(x, select=c('variant_id','phenotype_id','slope','tstat','pval_nominal'))
  names(x)<-c('SNP',     'gene',    'beta',    't-stat',  'p-value')
  write.table(x, file=paste0('01.',ct,'.assoc.txt'), quote=F, col.names=T, row.names=F,sep='\t')
  system(paste0('bgzip -f ','01.',ct,'.assoc.txt'))
}

snpanno_make<- function(){
  # x<- matrix(0,length(unique(df1$variant_id)),length(unique(df.ann$Antibody))+1) %>%
  #   as.data.frame()
  # names(x)<-c('SNP',unique(df.ann$Antibody))
  # x$SNP<- unique(df1$variant_id)
  # for(s in unique(df1$variant_id)){
  #   df.s<- subset(df1, variant_id==s, select=c('variant_id','CHR','Pos'))%>% unique()
  #   df.ann.s<- subset(df.ann, CHR==df.s$CHR & start <=df.s$Pos & End >= df.s$Pos)
  #   if(nrow(df.ann.s)>0){
  #     x[x$SNP==s,names(x) %in% df.ann.s$Antibody]<-1
  #   }
  # }
  dt1 <- data.table(df1)
  dt.ann <- data.table(df.ann)
  
  # Get unique SNP and Antibody column lengths
  num_snps <- length(unique(dt1$variant_id))
  num_antibodies <- length(unique(dt.ann$Antibody))
  
  # Create a data table with all zeros
  x <- data.table(SNP = unique(dt1$variant_id))
  x[, (unique(dt.ann$Antibody)) := 0L]
  
  # Set column names
  setnames(x, c('SNP', unique(dt.ann$Antibody)))
  
  # Perform conditional filtering and assignment
  df.ann.s <- dt.ann[dt1, on = .(CHR, start <= Pos, End >= Pos), nomatch = 0]
  df.ann.s <- unique(df.ann.s, by = c('CHR', 'variant_id', 'start', 'Antibody'))
  for(a in unique(df.ann.s$Antibody)){
    df.ann.s1<- subset(df.ann.s, Antibody==a)
    x[x$SNP %in% df.ann.s1$variant_id,which(names(x)==a)]<- 1
  }
  
  write.table(x, file=paste0('01.',ct,'.snpanno.',cell,'.txt'), quote=F, col.names=T, row.names=F,sep='\t')
  system(paste0('bgzip -f 01.',ct,'.snpanno.',cell,'.txt'))
}

genemap_make<- function(){
  x<- df.bed%>% filter( V4 %in% df1$phenotype_id) %>%
    select(paste0('V', seq(4)))
  write.table(x, file=paste0('01.',ct,'.genemap.txt'), quote=F, col.names=F, row.names=F,sep='\t')
  system(paste0('bgzip -f ','01.',ct,'.genemap.txt'))
}


args<- commandArgs(TRUE)
ct=args[1]
cell=args[2]
cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')

df<- fread('all_cts.apaQTL.signif.ps.csv', header=T, sep=',') %>%
  group_by(celltype,variant_id ,phenotype_id) %>%
  arrange(pval_nominal)%>%
  slice_head(n=1)

df1<- df %>% filter(celltype==ct)
snpmap_make()
assc_make()

df.ann1<-fread('TF.all.sorted.bed', header=T) ##download from ENCODE 
df.ann2<-fread('RBP.all.sorted.bed', header=T) ##download from ENCODE
df.ann<- rbind(data.frame(df.ann1, type='TF'), data.frame(df.ann2, type='RBP'))%>%
  mutate(Antibody=paste0(str_split_i(Antibody, '\\-',1),'.',type,'_D'))
snpanno_make()
df.bed<- fread('GencodeV32.20230414.hg38_3UTR_annotation.uniq.bed', header=F)
genemap_make()
