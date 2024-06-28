#!/bin/R
library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)

###compare APA between different groups

DAPA_identify<- function(df, ct){
  df <- df %>% as.data.frame()%>%
    mutate(ADdiag3types= case_when(ADdiag3types=="nonAD" ~ 0,
                                   ADdiag3types=="earlyAD"~1,
                                   ADdiag3types=="lateAD"~2)) ## recode as 0, 1, 2
  
  signif_test<- function(i,df){
    mat<- df
    mat1 <- mat %>% 
      select(library_id,ADdiag3types, names(mat)[i+7], pmi, msex, age_death)%>%
      na.omit()
    names(mat1)[3]<- 'PDUI_raw'
    model<- lm(mat1$PDUI_raw ~ mat1$pmi+mat1$msex+mat1$age_death) ## corrected by potential factors
    
    mat2<- data.frame(mat1,PDUI_adj=model$residuals)
    
    attach(mat2)
    df.cor.test<- cor.test(ADdiag3types, PDUI_adj, method='pearson')
    df.cor<- data.frame(phenotype_id=names(mat)[i+7],slope=df.cor.test$estimate, pval.cor=df.cor.test$p.value)
    detach(mat2)
    
    dignif_pair<- function(mat2, g1, g2){
      x<- mat2 %>% filter(ADdiag3types %in% c(g1, g2))
      x$PDUI_adj <- as.numeric(x$PDUI_adj)
      attach(x)
      x.pair <- wilcox.test(PDUI_adj ~ ADdiag3types)
      detach(x)
      
      df.pair.signif <- data.frame(log2FC=log2(mean(x[x$ADdiag3types==g2,]$PDUI_raw)/mean(x[x$ADdiag3types==g1,]$PDUI_raw)),
                                   pval.pair=x.pair$p.value)
      return(df.pair.signif)
    }
    df.pair.01 <- dignif_pair(mat2,0,1)
    names(df.pair.01)<- paste(names(df.pair.01),'nonAD_earlyAD',sep='.')
    df.pair.02 <- dignif_pair(mat2,0,2)
    names(df.pair.02)<- paste(names(df.pair.02),'nonAD_lateAD',sep='.')
    df.pair.12 <- dignif_pair(mat2,1,2)
    names(df.pair.12)<- paste(names(df.pair.12),'earlyAD_lateAD',sep='.')
    df.signif.all <- data.frame(df.cor, df.pair.01, df.pair.12, df.pair.02)
    return(df.signif.all)
  }
  
  df.signif.merge<- data.frame()
  for(i in seq(ncol(df)-7)){
    df.signif.i<- signif_test(i,df)
    df.signif.merge<- rbind(df.signif.merge, df.signif.i)
  }
  return(df.signif.merge)
}

###Main function#####
args=commandArgs(TRUE)
ct=args[1]

###for AD-DAPA identification
df<- fread(paste0(ct,'.raw_PDUI.csv') ,sep=",",header=T, check.names=F)
df.dapa<- DAPA_identify(df, ct)
pval.list<- names(df.dapa)[grep('pval',names(df.dapa) )]

df.fdr<- df.dapa %>% select(pval.list) 
df.fdr<-apply(df.fdr, 2, function(x){p.adjust(x,method='fdr')}) %>% as.data.frame()
names(df.fdr)<- gsub('pval','fdr',names(df.fdr))
df.dapa<- data.frame(df.dapa, df.fdr)
write.csv(df.dapa, file=paste0(ct,'.DAPA.cor_pair_test.csv'))

###for AD-DEG identification
df<- fread(paste0(ct,'.raw_expression.csv') ,sep=",",header=T, check.names=F)
df.deg<- DAPA_identify(df, ct)
pval.list<- names(df.deg)[grep('pval',names(df.deg) )]

df.fdr<- df.deg %>% select(pval.list) 
df.fdr<-apply(df.fdr, 2, function(x){p.adjust(x,method='fdr')}) %>% as.data.frame()
names(df.fdr)<- gsub('pval','fdr',names(df.fdr))
df.deg<- data.frame(df.deg, df.fdr)
write.csv(df.deg, file=paste0(ct,'.DEG.cor_pair_test.csv'))