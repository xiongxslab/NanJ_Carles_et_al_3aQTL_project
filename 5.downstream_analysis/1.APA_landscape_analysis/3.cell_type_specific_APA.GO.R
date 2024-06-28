library(clusterProfiler)
library(ggplot2)
library(data.table)
library(org.Hs.eg.db)
library(stringr)
#args<- commandArgs(TRUE)
#cutoff<- args[1]
cts.list = c("Neuron","Exc","Inh","Ast","Oli","Opc","Mic")
res<- fread('../all_cts.specific.res',header=T)
names(res)<- c("transcript","celltype","type",cts.list[-1])

###perform GO enrichment analysis
go_perform<- function(){
  res<- fread('all_cts.specific.res',header=T)
  gene.sum<- fread('../../01.all_cell-types.gene.sum',header=F, sep=" ")
  gene.fre<- as.matrix(table(gene.sum$V2))
  gene.useful<-names(gene.fre[gene.fre[,1]>=4, ])
  bggene<- unique(str_split_i(gene.useful, '\\|',2))

  ego.merge<- data.frame()
  for (tp in c("lengthening", "shortening" )){
    for (ct in cts.list){
      gen<- subset(res, celltype %in% ct & type %in% tp)
      gene<- unique(str_split_i(gen$transcript, '\\|',2))
      ego <- enrichGO(
        gene          = gene,
        keyType = "SYMBOL",
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1,
        readable      = TRUE,
        universe      = bggene
      )
      if(nrow(ego)>0){
        ego.merge<- rbind(ego.merge, data.frame(ct, tp, as.data.frame(ego)))
      }
    }
  }
  return(ego.merge)
}

###de-duplication for GO results
gene_compare<- function(a,b, cutoff){
  gen_split<- function(x){str_split(x, '\\/') %>% unlist}
  gen.dup<- intersect(gen_split(a), gen_split(b))
  if(length(gen.dup)/min(length(gen_split(a)), length(gen_split(b))) < cutoff){
    g.com<- TRUE
  }else{
    g.com<- FALSE
  }
}

go_select<- function(df, ct, tp){
  df.filter<- df %>%
    filter(celltype==ct, type==tp) %>%
    mutate(rank=frankv(pvalue, ties.method="min"))
  if(nrow(df.filter)==1){
    df.use<- df.filter
  }else{
    df.filter<- df.filter[order(df.filter$rank),] ## sort by rank
    df.use <- df.filter[1,]
    for(index in 2:nrow(df.filter)){
      df.use1<- df.filter[index,]
      num_rows <- nrow(df.use)
      go.compare<- list()
      for(i in 1:num_rows){
        go.compare1<- gene_compare(df.use1[[1, 11]],df.use[[i, 11]], 0.5)
        go.compare<- c(go.compare, go.compare1)
      }
      if(all(go.compare)){
        df.use<- rbind(df.use, df.use1)
      }else{
        df.use<- df.use
      }
    }
  }
  return(df.use)
}
##################Main function##############################
for (tp in c("lengthening", "shortening" )){
  for (ct in cts.list){
    go_perform(tp, ct)
  }
}

df<- go_perform() %>%
  filter(pvalue <0.01)
df.filter.merge<- data.frame()
for(ct in unique(df$celltype)){
  for(tp in unique(df$type)){
    print(paste(ct,tp))
    df.go<- go_select(df, ct, tp)
    print(nrow(df.go))
    if(nrow(df.go)>0){
      df.filter.merge<- rbind(df.filter.merge, df.go)
    }else{
      df.filter.merge<-df.filter.merge
    }
  }
}

write.csv(df.filter.merge, file='02.240402.all_cts.specific_APA.dedup.GO.csv')

