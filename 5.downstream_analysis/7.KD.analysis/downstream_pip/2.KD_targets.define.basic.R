library(dplyr)
library(readr)
library(data.table)
library(dplyr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(UpSetR)
# library(ggplot2)
args = commandArgs(TRUE)

sample_list = args[1]
wkdir=args[2]
gene=args[3]
cutoff=0.25

setwd(wkdir)
df.pdui = read_rds('06.PDUI_analys/02.all_samples.PDUI_melt.rds')%>%
  filter(shGene==gene)%>%
  mutate(Is_target=ifelse(abs(log2FC)> cutoff & abs(delt_PDUI) > cutoff, 'Yes','No'))

### step1 define the intersection between reps as the shID_merged's targets
df.pdui.sum = df.pdui %>%
  filter(Is_target=='Yes')%>%
  group_by(Gene, Loci, fit_value,shID_merged)%>%
  summarise(N=n())%>%
  filter(N==2) %>% ## defined the intersection between replicates as the targets
  mutate(shGene = str_split_i(shID_merged, '\\-',1),
         index= paste(shID_merged, Gene))
df.pdui.sh = df.pdui %>%
  group_by(Gene, Loci, fit_value,shID_merged)%>%
  summarise(PDUI.NTC=mean(PDUI.NTC), PDUI.sh=mean(PDUI.sh),
            log2FC=log2(PDUI.sh/PDUI.NTC),
            delt_PDUI=PDUI.sh - PDUI.NTC)%>%
  mutate(Is_target = ifelse(paste(shID_merged, Gene) %in% unique(df.pdui.sum$index)& abs(log2FC)>cutoff&abs(log2FC)>cutoff,'Yes','No'))
saveRDS(df.pdui.sh, file=paste0('06.PDUI_analys/03.',gene,'.reps_merged.PDUI_melt.rds'))

### step2 define the intersection between shs as the shs_merged's targets
all_combos = combn(unique(df.pdui.sum$shID_merged),2)
sorted_combos <- apply(all_combos, 2, sort)

# 3. 转置后使用unique()去重
unique_sorted_combos <- unique(sorted_combos)
df.shs_merged = data.frame()
for(i in seq(ncol(unique_sorted_combos))){
  shs = unlist(unique_sorted_combos[,i])
  df.y = subset(df.pdui.sum, shID_merged %in% shs)%>%
    group_by(Gene, Loci, fit_value)%>%
    summarise(N=n())%>%
    filter(N==2)%>%
    mutate(shs_merged=paste(shs,collapse = '&'),
           index =paste( shs_merged,Gene))
  df.pdui.sh1 = df.pdui.sh %>%
    filter(shID_merged  %in% shs)%>%
    group_by(Gene, Loci, fit_value)%>%
    summarise(PDUI.NTC=mean(PDUI.NTC), PDUI.sh=mean(PDUI.sh),
              log2FC=log2(PDUI.sh/PDUI.NTC),
              delt_PDUI=PDUI.sh - PDUI.NTC)%>%
    mutate(shs_merged =paste(shs,collapse = '&'),
      Is_target = ifelse(paste(shs_merged, Gene) %in% unique(df.y$index)& abs(log2FC)>cutoff&abs(log2FC)>cutoff,'Yes','No'))
  
  df.shs_merged = rbind(df.shs_merged, df.pdui.sh1)
}

saveRDS(df.shs_merged, file=paste0('06.PDUI_analys/03.',gene,'.shs_merged.PDUI_melt.rds'))
