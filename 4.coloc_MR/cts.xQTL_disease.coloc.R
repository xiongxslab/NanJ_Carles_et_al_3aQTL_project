library(data.table)
library(dplyr)
library(coloc)
library(stringr)

## 1. make molQTL files
df.a_make<- function(ct){
  
  qtl.apa.sig<- fread(paste0('01.',ct,'.xQTL.signif.txt'),header=T)
  df.a<- fread(paste0(apa_dir,ct,'.xQTL.sumstat_allele.txt.gz'),header=T) %>%
    filter(phenotype_id %in%  qtl.apa.sig$phenotype_id)
  info.a<- read.table('xQTL_calling_sampleSize.sum',header=F, ##sample size for xQTL
                      col.names=c('celltype','N','no.gene'))
  df.a<- data.frame(Gene=str_split_i(df.a$phenotype_id,'\\|',2),
                    N=info.a[info.a$celltype==ct,]$N, df.a)
  df.p<- fread(paste0(disease,'.sumstats.gz'), header=T) %>%
    mutate(N=ncases+ncontrols)
  input<- merge(df.a, df.p, by.x='variant_id', by.y='SNP', suffixes = c('.mQTL','.GWAS')) %>%
    group_by(variant_id, phenotype_id)%>%
    filter(P==min(P),(A1==Affected_allele & A2==Reference_allele)|
             (A2==Affected_allele & A1==Reference_allele)) %>%
    mutate(A1=ifelse(A1== Affected_allele, A1, A2),
           A2=ifelse(A1== Affected_allele, A2, A1),
           beta=ifelse(A1==Affected_allele, beta, -beta))
  saveRDS(input,paste0('1.',ct,'.xQTL.',disease,'.input.rds'))
}

ct_coloc<- function(ct,disease, i,type.p){
  input<- readRDS(paste0('1.',ct,'.xQTL.',disease,'.input.rds')) %>%
    filter(CHR %in% i)
  if(nrow(subset(input, P<1e-4))>0){
    ##select the significant genes for coloc
    for (phen in unique(subset(input,P<1e-4)$phenotype_id)){
      input1<- input %>%
        filter(phenotype_id %in% phen) %>%
        select(-Gene) %>%
        na.omit()
      
      if(nrow(subset(input1, P<1e-4))>0){
        print(phen)
        if(type.p=='quanti'){
          result<- coloc.abf(dataset2=list(snp=input1$variant_id, beta=input1$slope, varbeta=(input1$slope_se)^2, 
                                           MAF=input1$maf.mQTL, N=input1$N.mQTL, type='quant'),
                             dataset1=list(snp=input1$variant_id, beta=input1$beta, varbeta=(input1$SE)^2, 
                                           MAF=input1$maf.GWAS, N=N.GWAS, type='quant'))
        }else{
          result<- coloc.abf(dataset1=list(snp=input1$variant_id, S=input1$ncases/input1$N.GWAS, beta=input1$beta, varbeta=(input1$SE)^2,
                                           MAF=input1$maf.GWAS, N=input1$N.GWAS, type='cc'),
                             dataset2=list(snp=input1$variant_id, beta=input1$slope, varbeta=(input1$slope_se)^2, 
                                           MAF=input1$maf.mQTL, N=input1$N.mQTL, type='quant'))
        }
        
        write.table(cbind(ct, phen, t(result$summary),result$results[which.max(result$results$SNP.PP.H4),1]), 
                    file=paste0('01.',ct,'.xQTL_',disease,'.coloc.res'),sep='\t',
                    quote=F, col.names=F, row.names=F,append=T)
        
      }
    }
  }

}
args<- commandArgs(T)
ct=args[1]
disease=args[2]
chr=args[3]
tp=args[4] #type of disease

df.a_make(args[1])
ct_coloc(args[1],args[2], args[3], args[4])

