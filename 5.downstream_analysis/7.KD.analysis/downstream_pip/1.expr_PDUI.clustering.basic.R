library(dplyr)
library(readr)
library(data.table)
library(dplyr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(reshape2)
# library(ggplot2)
args = commandArgs(TRUE)
sample_list = args[1]
wkdir=args[2]
# sample_list = 'samples.list'
setwd(wkdir)
## step 1 generate expression matrix
df.id <- fread(sample_list,header=F)
############################################################################################
################ mRNA based-analysis starting ###############################################
############################################################################################
if(file.exists('05.expr_analys/01.all_samples.expr_mat.rds')){
  print('expression matrix is existed...')
  df.Merge = read_rds('05.expr_analys/01.all_samples.expr_mat.rds')
  df.Merge = df.Merge %>%
    select(names(df.Merge)[1:3],df.id$V1)
}else{
  df.Merge <- data.frame()
  for(id in df.id$V1){
    if(!file.exists(paste0('04.count/',id,'.stringtie.tab'))){
      print(paste(id,'does not exist...'))
    }else{
      df.x <- fread(paste0('04.count/',id,'.stringtie.tab'),header=T)%>%
        mutate(index=paste(Reference,Start,End, Strand,sep=':'))%>%
        select(index,`Gene ID`,`Gene Name`,'TPM')
      names(df.x)[4]<- id
      if(nrow(df.Merge)>0){
        df.Merge <- merge(df.Merge, df.x,by=c('index','Gene ID','Gene Name'),all.x=T, all.y=T)
      }else{
        df.Merge <- df.x
      }
    }
  }
  df.Merge = df.Merge[which(rowMeans(df.Merge[,-c(1:3)])>2),] ## filter the low expression genes
  if(!dir.exists('05.expr_analys')){
    dir.create('05.expr_analys')
  }
  # names(df.Merge)[4:6] = paste0(names(df.Merge)[4:6],'-3')
  # names(df.Merge) = gsub('PABPN2','PABPN1',names(df.Merge))
  # names(df.Merge)[grep('SF3BF-1-MIX',names(df.Merge))] = 'SF3B4-1-MIX'
  saveRDS(df.Merge, file='05.expr_analys/01.all_samples.expr_mat.rds')
  print('expression matrix has saved...')
}
# df.Merge = read_rds('05.expr_analys/02.all_samples.useful.merged.rds')

# # df.Merge.expr = read_rds('05.expr_analys/01.samples.last_version.expr_mat.rds')
# df.old = read_rds('../PLKO.ri.batch1/05.expr_analys/01.samples.last_version.expr_mat.rds')
# names(df.old)[4:5] = c('NTC-1-1','NTC-1-2')
# names(df.old)[-c(1:3)] = paste0(names(df.old)[-c(1:3)],'-old')
# df.Merge.expr = merge(df.Merge.expr,df.old[,-c(14:15)], by=names(df.Merge.expr)[1:3])


#### step 2 calculate mRNA-based correlation 
df.cor = cor(df.Merge[,-c(1:3)])
if(!(dir.create('0.imgs'))){
  dir.create('0.imgs')
}
if(!(dir.create(paste0('0.imgs/',sample_list)))){
  dir.create(paste0('0.imgs/',sample_list))
}

# pheatmap(df.cor, file=paste0('0.imgs/',sample_list,'/01.expr_cor.Merged.heatmap.pdf'),height=9, width=9.5,display_numbers = F,
#          number_format = "%.2f", main = 'mRNA expression based clustering',number_color = 'white',border_color='white',
#          color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100))
# 
# 
# #### step 3 show the gene * samples matrix
# df.expr.plot = df.Merge[,-c(1:3)]%>% 
#   filter(rowMeans(.) > 2)%>% t()%>%
#   as.data.frame()
# pheatmap(log2(df.expr.plot+1) , file=paste0('0.imgs/',sample_list,'/01.expr_cor.Merged.expr.heatmap.pdf'),height=5.5, width=9.5,display_numbers = F,
#          number_format = "%.2f", main = 'mRNA expression based clustering',number_color = 'white',border_color='white',
#          clustering_method ="average", show_colnames = F,cluster_rows = T, cluster_cols = T )
#### step 4 calculate log2FC-based correlation with the corresponding NTC (all reps) as the background (aim to control the batch effect)
df.sm = fread('sample_change.txt', header=T)
df.Merge.m = melt(df.Merge, c('index','Gene ID','Gene Name'), variable.name = 'shID', value.name ='TPM')
df.Merge.m = merge(df.Merge.m, df.sm, by.x='shID', by.y='names_raw')%>%
  mutate(shID=names_new,
    shGene=str_split_i(shID,'\\-',1), batchID = str_split_i(shID,'\\-',3))%>%
  select(-names_new)
df.ntc = subset(df.Merge.m, shGene=='NTC')%>%
  group_by(index,`Gene ID`,`Gene Name`, batchID)%>%
  summarise(TPM=mean(TPM))
df.Merge.m = merge(df.ntc, df.Merge.m %>% filter(shGene!='NTC'),
                   by=c('index','Gene ID','Gene Name','batchID'),
                   suffixes = c('.NTC','.sh'))%>%
  mutate(log2FC=log2(TPM.sh/TPM.NTC))
df.Merge.d = dcast(df.Merge.m,index+`Gene ID`+`Gene Name`~shID, value.var = 'log2FC', fun.aggregate = mean)
df.cor = cor(df.Merge.d[,-c(1:3)], method='pearson')

pheatmap(df.cor , file=paste0('0.imgs/',sample_list,'/01.expr_cor.Merged.expr_log2FC.heatmap.pdf'),height=5.5, width=9.5,display_numbers = F,
         number_format = "%.2f", main = 'mRNA expression based clustering',number_color = 'white',border_color='white',
         clustering_method ="average", show_colnames = F,cluster_rows = T, cluster_cols = T )
## step 5 calculate KD efficiencies (not used) 
df.gene = subset(df.Merge.m, `Gene Name` %in% unique(df.Merge.m$shGene))%>%
  mutate(FC=round(TPM.sh/TPM.NTC,3))
write.csv(df.gene, file=paste0('0.imgs/',sample_list,'/00.all_gene.KD.eff.csv'))

############################################################################################
################ APA based-analysi starting ################################################
############################################################################################

## step 1 merge the PDUI matrix
if(file.exists('06.PDUI_analys/01.all_samples.PDUI_mat.rds')){
  print('PDUI matrix is existed...')
  df.Merge = read_rds('06.PDUI_analys/01.all_samples.PDUI_mat.rds')
  df.Merge = df.Merge %>%
    select(names(df.Merge)[1:4],df.id$V1)
}else{
  df.Merge <- data.frame()
  for(i in c(seq(22),'X','Y')){
    df.x <- fread(paste0('Dapars2_chr',i,'/Dapars2_result_temp.chr',i,'.txt'),header=T)
    if(nrow(df.Merge)==0){
      df.Merge<- df.x
    }else{
      df.Merge <- rbind(df.Merge, df.x)
    }
  }
  names(df.Merge) = gsub("_PDUI","",gsub("04.bigwig/",'',names(df.Merge)))
  df.Merge = na.omit(df.Merge)
  if(!dir.exists('06.PDUI_analys')){
    dir.create('06.PDUI_analys')
  }
  
  saveRDS(df.Merge, file='06.PDUI_analys/01.all_samples.PDUI_mat.rds')
}


# #### step 2 calculate mRNA-based correlation 
# df.cor = cor(df.Merge[,-c(1:4)])
# if(!(dir.create(paste0('0.imgs/',sample_list)))){
#   dir.create(paste0('0.imgs/',sample_list))
# }
# 
# pheatmap(df.cor, file=paste0('0.imgs/',sample_list,'/01.PDUI_cor.Merged.heatmap.pdf'),height=9, width=9.5,display_numbers = F,
#          number_format = "%.2f", main = 'APA expression based clustering',number_color = 'white',border_color='white',
#          color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100))
# 
# #### step 3 show the gene * samples matrix
# pheatmap(df.Merge[,-c(1:4)]%>% na.omit() , file=paste0('0.imgs/',sample_list,'/01.PDUI_cor.Merged.PDUI.heatmap.pdf'),height=5.5, width=9.5,display_numbers = F,
#          number_format = "%.2f", main = 'PDUI based clustering',number_color = 'white',border_color='white',
#          show_colnames = F,cluster_rows = T, cluster_cols = T)
#### step 4 calculate log2FC-based correlation with the corresponding NTC (all reps) as the background
if(!file.exists('06.PDUI_analys/02.all_samples.PDUI_melt.rds')){
  df.melt.pdui = df.Merge %>%
    melt(c('Gene','Loci','fit_value','Predicted_Proximal_APA'), value.name='PDUI', variable.name='si')%>%
    mutate(siGene=str_split_i(si,'\\-',1) ,
           batchID = str_split_i(si,'\\-',4) ,
           PDUI = PDUI+ 0.001)
  
  df.Merge.m = melt(df.Merge, c('Gene','Loci','fit_value','Predicted_Proximal_APA'), variable.name = 'shID', value.name ='PDUI')%>%
    mutate(shGene=str_split_i(shID,'\\-',1), batchID = str_split_i(shID,'\\-',3),PDUI= PDUI+ 0.001)
  df.ntc = subset(df.Merge.m, shGene=='NTC')%>%
    group_by(Gene,Loci,fit_value,Predicted_Proximal_APA)%>%
    summarise(PDUI=mean(PDUI))
  df.Merge.m = merge(df.ntc, df.Merge.m %>% filter(shGene!='NTC'),
                     by=c('Gene','Loci','fit_value','Predicted_Proximal_APA'),
                     suffixes = c('.NTC','.sh'))%>%
    mutate(log2FC=log2(PDUI.sh/PDUI.NTC),delt_PDUI=PDUI.sh-PDUI.NTC, 
           Is_target =ifelse(abs(log2FC)>0.15&abs(delt_PDUI)>0.15,'Yes','No'))
  df.Merge.m$shID_merged= paste(df.Merge.m$shGene, str_split_i(df.Merge.m$shID,'\\-',2), str_split_i(df.Merge.m$shID,'\\-',4),sep ='-')
  saveRDS(df.Merge.m ,file = '06.PDUI_analys/02.all_samples.PDUI_melt.rds')
}else{
  df.Merge.m = read_rds('06.PDUI_analys/02.all_samples.PDUI_melt.rds')%>%
    filter(shID %in%  df.id$V1)
}
df.Merge.d = dcast(df.Merge.m, Gene+ Loci+ fit_value + Predicted_Proximal_APA ~ shID, value.var = 'log2FC', fun.aggregate=mean)
df.cor = cor(df.Merge.d[,-c(1:4)], method='spearman')
pheatmap(df.cor, file=paste0('0.imgs/',sample_list,'/01.PDUI_cor.Merged.PDUI_log2FC.heatmap.pdf'),height=5.5, width=9.5,display_numbers = F,
         number_format = "%.2f", main = 'PDUI based clustering',number_color = 'white',border_color='white',
         clustering_method ="average", show_colnames = F,cluster_rows = T, cluster_cols = T )         
         