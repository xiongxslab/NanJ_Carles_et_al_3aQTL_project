library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(foreach)
library(doParallel)

perform_mr <- function(phen,b_file, disease){
  mqtl_exp_dat <- input %>% filter(phenotype_id %in% phen ) %>%
    mutate(pheno=disease)
  if(length(mqtl_exp_dat$variant_id) >= 1){
    check <- 1
    if(length(mqtl_exp_dat$variant_id) > 1){
      check <- try(
        mqtl_exp_dat <- ld_clump_local(mqtl_exp_dat%>% dplyr::rename(rsid = variant_id, pval=pval_nominal), clump_kb = 100, clump_r2 = 0.01, clump_p = 1,
                                       bfile = b_file, plink_bin = get_plink_exe()) %>% dplyr::rename(variant_id = rsid),
        silent = FALSE
      )
    }
    if("try-error" %in% class(check)){
      next
    }
    if( nrow(mqtl_exp_dat)>0){
      mqtl_exp_dat<- mqtl_exp_dat %>%
        filter(variant_id %in% mqtl_exp_dat$variant_id)
      mqtl_exp_dat1 <- format_data(mqtl_exp_dat, type = "exposure", phenotype_col = "phenotype_id", snp_col = "variant_id", beta_col = "slope", se_col = "slope_se",
                                  effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "maf.mQTL", pval_col = "pval_nominal")
      gwas_out_dat <- format_data(mqtl_exp_dat, type = "outcome", phenotype_col = "AD", snp_col = "variant_id", beta_col = "beta", se_col = "SE",
                                  effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "maf.GWAS", pval_col = "P",
                                  samplesize_col = "N.GWAS", ncase = "ncases", ncontrol = "ncontrols")
      mr_dat <- harmonise_data(exposure_dat = mqtl_exp_dat1, outcome_dat = gwas_out_dat) %>% filter(mr_keep == TRUE)
      if(length(mr_dat$SNP) >= 1){
        mr_res <- mr(mr_dat, method_list = c("mr_wald_ratio", "mr_ivw"))
        mr_res <- mr_res %>% select(-id.exposure, -id.outcome) %>% select(exposure, everything())
        return(mr_res)
                   	
      }
    }
  }
  
}


args <- commandArgs(T)
ct <- args[1]
i<- args[2]
disease <- args[3]

b_file='EUR'
snp.1k<- fread('EUR.frq',header=T) %>%
  select('SNP')
qtl.apa.sig<- fread(paste0(ct,'.apaQTL.signif.txt'),header=T) %>%
  mutate(index=paste0(phenotype_id, variant_id, sep=':')) 
input<- readRDS(paste0('1.',ct,'.',disease,'.input.rds')) 
input$index<- paste0(input$phenotype_id, input$variant_id, sep=':')
if(length(intersect(input$CHR, i))>0){
  print(paste(i,'of', ct, 'analysis starting ...'))
  input<- subset(input, (CHR %in% i) & (index %in% qtl.apa.sig$index) &(variant_id %in% snp.1k$SNP))
  input$pval_nominal<- as.numeric(input$pval_nominal)
  input$P <- as.numeric(input$P)
  res<- data.frame()
  for(phen in unique(input$phenotype_id)){
    res1 <- perform_mr(phen, 'EUR',disease)
    print(paste(phen, 'EUR'))
    if(nrow(res1)>0){
      res<-data.frame(celltype=ct,res1)
      write.table(res1, file=paste0('02.',ct,'.', disease,'.MR'), quote = FALSE, sep = "\t", row.names = FALSE,
                               col.names=F,append=T)
    }
  }
  
}else{
  print(paste(i,'of', ct, 'has no significant loci, exiting ...'))
}



