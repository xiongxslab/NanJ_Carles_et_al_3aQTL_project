#!/bin/R
library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

#####################################################################################
##################function define####################################################
#####################################################################################
###step 1 aggregate the depth across all individuals
dep_prepare<- function(ct, gene){
  #### merge the genotype and bigwig files
  df<- fread("../all_samples.bedGraph", header=F)
  names(df)<- c('celltype','library_id','CHR','start','end','score')
  utr_bed<- fread(paste0('/data/slurm/nanjh/scAPA_project/05APAmatrix/01.sample_QC_matrix/','07',ct,'_imputation_MI_1.csv'),header=T, sep=',')
  utr_bed<- subset(utr_bed, str_split_i(Gene,'\\|',1) %in% gene1, select=c('Gene','Predicted_Proximal_APA','Loci'))##select the interested isoform
  chr=unique(str_split_i(utr_bed$Loci, '\\:',1))
  ps_start=min(as.numeric(str_split_i(str_split_i(utr_bed$Loci, '\\:',2), "-",1)))
  ps_end=max(as.numeric(str_split_i(str_split_i(utr_bed$Loci, '\\:',2), "-",2)))
  df<- subset(df,  celltype==ct)
  ##summary the mean value of 3'UTR across individuals
  Merge <- df %>%
    filter(start > ps_start-100 & start <ps_end+1000)%>%
    group_by(CHR,start) %>%
    # filter(score_adj > quantile(score_adj, 0.05) &score_adj< quantile(score_adj, 0.95)) %>%
    summarise(score_scale=mean(score))
  return(Merge)
}


## step 2 identify the start and end position of independent peak
find_peak<- function(inputfile, spacing_size,  window_size,cutoff){
  
  inputfile$smooth <- apply(as.matrix(inputfile$score_scale), 2, function(x) zoo::rollapply(x, window_size, mean, align = "center", fill = 0))
  
  inputfile$diff <- c(0,diff(inputfile$smooth))
  inputfile$diff.smooth<-apply(as.matrix(inputfile$diff), 2, function(x) zoo::rollapply(x, window_size, mean, align = "center", fill = 0))
  
  inputfile$diff.smooth[which(abs(inputfile$diff.smooth)<=cutoff)]<-0
  turn.index<-which(inputfile$diff.smooth==0)
  
  # Find the start and end positions of consecutive values
  start_positions <- c()
  end_positions <- c()
  start_val <- turn.index[1]
  end_val <- turn.index[1]
  
  for (i in 2:length(turn.index)) {
    if (turn.index[i] == end_val + 1) {
      end_val <- turn.index[i]
    } else {
      start_positions <- c(start_positions, start_val)
      end_positions <- c(end_positions, end_val)
      start_val <- turn.index[i]
      end_val <- turn.index[i]
    }
  }
  start_index<-end_positions[which(inputfile$diff.smooth[end_positions+1]>0)]
  if(start_positions[1]==1)start_positions[1]=2
  end_index<-start_positions[which(inputfile$diff.smooth[start_positions-1]<0)]
  if(length(start_index)!=length(end_index)){
    end_index0<-which(inputfile$smooth<=0.05)[which(which(inputfile$smooth<=0.05)>start_index[length(start_index)])[1]]
    end_index<-c(end_index,end_index0 )
  }
  num<-min(length(start_index),length(end_index))
  peak_bed<- data.frame(CHR=unique(inputfile$CHR)[1],start=inputfile$start[start_index[1:num]],end=inputfile$start[end_index[1:num]])
  
  ##save the peak position files
  peak_bed$peaks<- paste0('Peak',rownames(peak_bed))
  write.table(peak_bed, file=paste0(gene, '.',ct,'_multiple_peak.bed'),
              sep='\t', col.names=F, row.names=F, quote=F)
  #visualization of all peaks
  ggplot()+
    #scale_color_manual(values=c("#FFD700", "#4169E1",'#EE82EE'))+
    # geom_line(position = "identity", alpha=1)+
    # geom_smooth(se=F, size=0.8)+
    theme_classic()  +
    theme(axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          strip.background = element_blank())+
    xlab(unique(inputfile$CHR)[1])+ggtitle(paste0(ct, ':',gene))+
    xlim(c(min(c(ps_start,inputfile$start)), max(c(ps_end,inputfile$start))))+
    # ylim(c(-0.05,0.5))+
    geom_area(data=inputfile, aes(x=start, y= score_scale), alpha=0.55, position='identity')+
    #geom_line(data=merge,aes(x=start,y=smooth), linewidth=0.5)+
    #  geom_line(data=tmp1,aes(x=start,y=diff.smooth))+
    geom_vline(xintercept =inputfile$start[start_index])+
    geom_vline(xintercept = inputfile$start[end_index], linetype='dashed')
   ggsave(paste0(gene,'.',ct,'.mutiple_peaks.pdf'), height = 4.5, width = 6)
  #peak_bed<- data.frame(CHR=unique(df$CHR), start=inputfile$start[start_index], end=inputfile$start[end_index])
   return(peak_bed)

}
##step 3 calculate sequence abundance of different peaks across individuals
abundance_cal<- function(disease,gene1,peak, info,ct){
  df<- fread("../all_samples.bedGraph", header=F) 
  names(df)<- c('celltype','library_id','CHR','start','end','score')
  df<- subset(df, celltype==ct) 
  genoty<- fread(paste0('/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/',disease,'/',gene1,'/',gene1, '.ped'),sep=' ',header=F)
  genoty$genotype<- paste0(genoty$V7,genoty$V8)
  geno<- merge(genoty[,c(2,9)], info[,c(1,3)], by.x='V2',by.y='projid')
  df<- unique(merge( df,geno,by='library_id',allow.cartesian=TRUE))
  ##calculate the mean abundance of every peaks
  z<- data.frame()
  for(i in seq(1,nrow(peak))){
    chr=peak[i,1]
    Start=peak[i,2]
    End=peak[i,3]
    Peak=peak[i,4]
    x <- subset(df, CHR==chr & start>Start & end<End )
    ##using weighted mean score to represent abundance
    y<- data.frame()
    for(sm in unique(x$library_id)){
      x1<- subset(x, library_id== sm)
      mean_score<- sum((x1$end-x1$start)*x1$score)/sum(x1$end-x1$start)
      y1<- data.frame(library_id=sm, projid=unique(x1$V2),genotype=unique(x1$genotype),mean_score=mean_score)
      names(y1)[4]<-paste0('Peak',i)
      y<- rbind(y, y1)
    }
    if(i==1){
      z<- y
    }else{
      z<- merge(z,y,by=c('library_id','genotype','projid'))
    }
  }
  return(z)
}
## step 4 calculate the ratio between two peaks
ratio_cal<- function(df1, s_peaks,gene, ct){
  library(ggsignif)
  #s_peak=4
  #df1<- df1[,1:6]
  df1$sum<- apply(df1[,-c(1:3)], 1, sum)
  for(i in seq(s_peak)){
    ##A/total
    ratio<-data.frame(df1[,which(names(df1)==paste0('Peak',i))]/df1$sum)
    names(ratio)<- paste0('Peak',i,'/total')
    mat<- subset(df1,select=paste0('Peak', seq(s_peak)[!(seq(s_peak)%in% i)]))
    for(j in as.numeric(gsub('Peak','',names(mat)))){
      if(i < j){
        #A/B; A/(A+B)
        x<- data.frame(df1[,which(names(df1)==paste0('Peak',i))]/df1[,which(names(df1)==paste0('Peak',j))],
                       df1[,which(names(df1)==paste0('Peak',i))]/(df1[,which(names(df1)==paste0('Peak',j))]+df1[,which(names(df1)==paste0('Peak',i))]))
        names(x)<- c(paste0('Peak',i,'/','Peak',j),
                     paste0('Peak',i,'/','(Peak',j,'+','Peak',i,')'))
        ratio<- data.frame(ratio,x, check.names = F)
        
      }else{
        ratio<- ratio
      }
    }
    df1<- cbind(df1,ratio)
    
  }
  ##calculate significance of different peaks using linear regression method
  df1[df1$genotype=='CC',]$genotype<- 0
  df1[df1$genotype=='TC',]$genotype<- 1
  df1[df1$genotype=='TT',]$genotype<- 2
  signif<- data.frame()
  for(n in names(df1)[-c(1:2)]){
    x<- apply(subset(df1, select=c('genotype',n)),2,as.numeric)
    s<- cor.test(x[,1], x[,2])
    signif<- rbind(signif, data.frame(peaks=n, signif=s$p.value, slope=s$estimate))
  }
  signif$fdr<- p.adjust(signif$signif, method='fdr')
  write.table(signif, paste0(gene, '.',ct,'.all_peaks.comparison_linear.res'), 
              sep='\t', quote=F, row.names=F)
  return(df1)
}
### step 5 prepare matrix file for four independent peak apaQTL calling
fast_QTL<- function(df,gene,ct){
  df4<- df %>%
    filter(!(library_id %in% c('D19-2458','D19-4136')))
  #df4<- merge(df4, info)
  #df4<- subset(df4, select=c(names(df),'projid'))
  mat<- data.frame()
  for(n in names(df4)[!(names(df4)%in% c('library_id','genotype','projid'))]){
    x<- unique(subset(df4, select=c(names(df)[1:2],n,'projid')))
    y<- as.data.frame(cbind(gsub('chr','',chr),ps_start, ps_end,paste0("SNCA","_",n),t(x[,3])))
    names(y)<- c("#Chr","start","end","Gene",x$projid)
    mat<- rbind(mat, y)
  }
  cov<- fread('/data/slurm/xiongxs/Proj_APA/2.cisQTL/2.voom/4.PEER.cor/Oli.PEER_covariates.merged.txt', header=T, check.names = F)
  mat<- subset(mat,select=c("#Chr","start","end","Gene", intersect(names(cov), names(mat))) )
  write.table(mat, file=paste0(gene, '.multi_peaks.ratio.bed'),quote=F,
                               row.names=F, sep='\t')
  write.table(subset(cov, select=c('ID',intersect(names(cov), names(mat)) )),
                     file=paste0(ct,'.PEER_covariates.merged.txt'),
              row.names=F, sep='\t', quote=F)
  ### perform apaQTL calling
  if(file.exists(paste0(gene, '.multi_peaks.ratio.bed.gz'))){
    file.remove(paste0(gene, '.multi_peaks.ratio.bed.gz'))
    file.remove(paste0(gene, '.multi_peaks.ratio.bed.gz.tbi'))
  }
  system(paste0("bgzip ",paste0(gene, '.multi_peaks.ratio.bed')))
  system(paste0("tabix -p bed ",paste0(gene, '.multi_peaks.ratio.bed.gz')))
  if(!(file.exists('output.1M'))){
    dir.create('output.1M')
  }
  qtl_sh=paste("/data/dingk/software/fastqtl-master/bin/fastQTL.static --vcf /data/slurm/nanjh/scAPA_project/01genotype_data/liftover/split_chr/10WGSfromBraod_match_samples_merged_chrall_qc_hg38_dedup.chr4.vcf.gz",
                '--bed',paste0(gene, '.multi_peaks.ratio.bed.gz'),
                '--region',ch,paste0('--out output.1M/',ct,'.',ch,'.allCov.out.txt.gz'),
               '--cov', paste0(ct,'.PEER_covariates.merged.txt'),
               '--window 1e6',sep=' ')
               
  system(qtl_sh)
}

### step 6 perform clocalization analysis
col_perf<- function(ct, ch, disease ){
  trait<- disease
  library(data.table)
  library(coloc)
  library(locuscomparer)
  gwas_dir='/data/slurm/licy/QTL_integration/1.Database/GWAS/Sumstats/'
  apa<- fread(paste0('output.1M/',ct,'.',ch,'.allCov.out.txt.gz'), header=F)
  map<- fread('/data/slurm/nanjh/scAPA_project/06.3aQTL/SNP.affect_allele', header=F)
  apa<- merge(map, apa, by.x='V1',by.y='V2')[,-c(7:9)]
  names(apa)<- c( 'variant_id',  'CHR','position','Reference_allele','Affected_allele', 'phenotype_id', 'maf', 'pval_nominal','slope','slope_se')
  ## obtain the sample size of each cell types
  aqtl.info<- read.csv('/data/slurm/nanjh/scAPA_project/05APAmatrix/01.sample_QC_matrix/05samples_qc_sum_percelltype.csv')
  ## coloc for muptiple diseases
  gwas.info<- fread("/data/slurm/nanjh/scAPA_project/10.coloc/gwas.info", header=T)
  gwas.info.trait<- subset(gwas.info, tags_GWAS==trait)
  qtl.trait<-fread(paste0(gwas_dir, trait,'.sumstats.gz'),header=T) 
  qtl.trait$A1<- toupper(qtl.trait$A1)
  qtl.trait$A2<- toupper(qtl.trait$A2)
  qtl.apa_gwas<- merge(qtl.trait, apa, by.x='SNP',by.y='variant_id', suffixes = c('_gwas',"_aqtl"))
  ## A1 of gwas data is the affected allele
  qtl.apa_gwas<- subset(qtl.apa_gwas, A1==Reference_allele& A2==Affected_allele | A2==Reference_allele& A1==Affected_allele)
  qtl.apa_gwas[toupper(qtl.apa_gwas$A1)!=toupper(qtl.apa_gwas$Affected_allele),]$BETA=-qtl.apa_gwas[toupper(qtl.apa_gwas$A1)!=toupper(qtl.apa_gwas$Affected_allele),]$BETA
  
  ## include the APA driven by genetics
  #qtl.apa.sig<- fread(paste0('/data/slurm/nanjh/scAPA_project/06.3aQTL/3aQTL/',ct,'.apaQTL.signif.txt'))
  #    write.table(cbind('celltype','Diseases','transcript','nsnps',paste0('PP',seq(0,4))), file=paste0('01.',ct,'_',trait,'.coloc.res'),
  #               quote=F,col.names=F, row.names=F, append=F,sep='\t')
  for (phen in unique(qtl.apa_gwas$phenotype_id)){
    input<- na.omit(unique(subset(qtl.apa_gwas, phenotype_id %in% phen & P<1e-4 &maf<1 &MAF<1))) ## filter the insignificant loci in gwas or whose maf=1
    if(nrow(input>0)){
      print(phen)
      input<- na.omit(unique(subset(qtl.apa_gwas, phenotype_id %in% phen )))
      if(gwas.info.trait[[1,5]]=='quanti'){
        result<- coloc.abf(dataset1=list(snp=input$SNP, beta=input$BETA, varbeta=(input$SE)^2, 
                                         MAF=input$MAF, N=gwas.info.trait[[1,4]], type='quant'),
                           dataset2=list(snp=input$SNP, beta=input$slope, varbeta=(input$slope_se)^2, 
                                         MAF=input$maf, N=aqtl.info[aqtl.info$Cell.type==ct,3], type='quant'))
      }else if (gwas.info.trait[[1,5]]=='cc'){
        result<- coloc.abf(dataset1=list(snp=input$SNP, S=as.numeric(gwas.info.trait[[1,2]])/gwas.info.trait[[1,4]], beta=input$BETA, varbeta=(input$SE)^2,
                                         MAF=input$MAF, N=gwas.info.trait[[1,4]], type='cc'),
                           dataset2=list(snp=input$SNP, beta=input$slope, varbeta=(input$slope_se)^2, 
                                         MAF=input$maf, N=aqtl.info[aqtl.info$Cell.type==ct,3], type='quant'))
      }
    }
    write.table(cbind(ct, trait,phen, t(result$summary),result$results[which.max(result$results$SNP.PP.H4),1]), 
                file=paste0('01.',ct,'_',trait,'.',chr,'.coloc.res'),sep='\t',
                quote=F, col.names=F, row.names=F,append=T)
    if(t(result$summary)[1,6]>0.5){
      res2<- subset(input, select=c('SNP','pval_nominal'), phenotype_id %in% phen)
      gwas2<- subset(input, select=c('SNP','P') , phenotype_id %in% phen)
      names(gwas2)<- c('rsid','pval')
      names(res2)<- c('rsid','pval')
      p.list<- locuscompare( res2,gwas2, marker_col1 = 'rsid', pval_col1 = 'pval',
                             title1 = "apaQTL", marker_col2 = 'rsid', pval_col2 = 'pval',
                             title2 = "GWAS", snp = NULL, population = "EUR", combine = F,
                             legend = TRUE, legend_position = 'topleft', lz_ylab_linebreak = FALSE, genome =  "hg38")
      print(p.list[[2]])
      ggsave(paste0(gsub('/','_',phen),'.apaQTL.locuscompare_local.pdf'),height=3,width=6)
      print(p.list[[3]])
      ggsave(paste0(gsub('/','_',phen),'.gwas.locuscompare_local.pdf'),height=3,width=6)
    }
  }
}

#####################################################################################
##################Main function######################################################
#####################################################################################
## re-find peaks for SNCA
ct='Oli'
ct='Exc'
gene='SNCA'
gene1='ENST00000618500.4'
disease='Parkinson_Disease'
utr_bed<- fread(paste0('/data/slurm/nanjh/scAPA_project/05APAmatrix/01.sample_QC_matrix/','07',ct,'_imputation_MI_1.csv'),header=T, sep=',')
utr_bed<- subset(utr_bed, str_split_i(Gene,'\\|',1) %in% gene1, select=c('Gene','Predicted_Proximal_APA','Loci'))##select the interested isoform

chr=unique(str_split_i(utr_bed$Loci, '\\:',1))
ch=gsub('chr','',chr)
ps_start=min(as.numeric(str_split_i(str_split_i(utr_bed$Loci, '\\:',2), "-",1)))
ps_end=max(as.numeric(str_split_i(str_split_i(utr_bed$Loci, '\\:',2), "-",2)))
info<- fread('/data/slurm/nanjh/scAPA_project/05APAmatrix/AD430_fastq_projid_mapping.tsv',header=T)
for(ct in c('Inh','Mic','Opc','Ast')){
  wkdir=paste0('/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/',ct)
  if(!(file.exists(wkdir))){
    dir.create(wkdir)
    setwd(wkdir)
  }else{
    setwd(wkdir)
  }
  Merge<- dep_prepare(ct, gene) #prepare depth data
  peak_bed<-find_peak(Merge,50,10,0.01) #re-define pA peaks
  df1<- abundance_cal(disease,gene1, peak_bed, info,ct) #estimate aubundance of every pA peaks
  df<- ratio_cal(df1,4,gene,ct) # calculate pA peaks ratios
  fast_QTL(df,gene,ct) #re-perform QTL mapping 
  col_perf(ct,ch,disease)#perform colocalization analysis
}



##2. LD block plot
library(LDheatmap)
library(snpStats)
library(data.table)
library(pheatmap)
library(chopsticks)
wkdir='/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/LDplot/'
setwd(wkdir)
ld<- fread("03.tag.snp.ld" , header=F)
rs<- fread("02.tag.SNP.bim" ,header=F)
#rs1<- subset(rs, V4>44.7*10^6& V4< 45.1*10^6)
ld<- as.data.frame(ld)
names(ld)<- rs$V2
row.names(ld)<- rs$V2
#ld1<- as.matrix(subset(ld,rownames(ld) %in% rs1$V2, select=rs1$V2))
#### Use an RGB pallete for the color scheme ####
rgb.palette <- colorRampPalette(c("#fee0d2", "#fc9272", "#de2d26"), space = "rgb")

pdf("04.tagsnp_LDBlock.pdf",width=8,height = 6)
LDheatmap(ld, SNP.name = names(ld), newpage=FALSE,
          add.map=F, add.key=F,color=rgb.palette(38))
dev.off()

pheatmap(ld[-10,-10], cluster_rows=F, cluster_cols=F, display_numbers=T,
         color=rgb.palette(30), filename='tagSNP.LD.pheatmap.pdf',
         height=5, width=5.3)


##generate all the colcoc results
Merge<- data.frame()
for(ct in c('Exc','Oli','Inh')){
  if(ct =='Oli'){
    df<- fread(paste0('/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/', ct,'/01.',ct,'_Parkinson_Disease.4.coloc.res'),header=F)
  }else{
    df<- fread(paste0('/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/', ct,'/01.',ct,'_Parkinson_Disease.chr4.coloc.res'),header=F)
    
  }
  df1<- df[,c(1,3,9)]
  names(df1)<- c('celltype','peaks','PP4')
  Merge<- rbind(Merge, df1)
}
saveRDS(Merge, file='SNCA.3cts.PP4.Merge.rds')
Merge<- readRDS('SNCA.3cts.PP4.Merge.rds')

ggplot(data=Merge, aes(x=celltype, y=peaks, size=PP4, color=PP4))+
  geom_point()+
  scale_color_gradient(low = "#FFEBCD",  high = "red3",name='PP4')+
  theme_bw()
ggsave(file='SNCA.3cts.PP4.pdf',height=5, width=4.5)

### extract all the lead variants P from GWAS and 3'aQTL
setwd('/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/')
Merge<- data.frame()
for(ct in c('Exc','Oli','Inh')){
  x<- fread(paste0(ct, '/','output.1M/',ct,'.4.allCov.out.txt.gz' ), header=F,
            col.names = c('phenotype_id','variant_id','distance','ma_samples','ma_count','maf','pval_nominal','slope','slope_se' )) %>%
    mutate(celltype=ct) %>%
    filter(variant_id %in% rs$V2[-10])
   Merge<- rbind(Merge, x)
}
saveRDS(Merge, file='231127.all_cts.SNCA_mutiPa.rds')
Merge1<- subset(Merge, phenotype_id %in% paste0('SNCA_',c('Peak3/(Peak4+Peak3)',
                                                          'Peak1/(Peak4+Peak1)')))
ggplot(data=Merge1,aes(x=variant_id, fill=celltype, y=-log10(pval_nominal)))+
  geom_bar(stat='identity',position='dodge')+
  theme_bw()+
  facet_wrap(~phenotype_id, ncol=1)+
  scale_fill_manual(values=c("#33A02C", "#B2DF8A",  "#FDBF6F"))
ggsave(file='all_leadVar.apaQTL.pdf')

disease='Parkinson_Disease'
df<-fread(paste0('/data/slurm/licy/QTL_integration/1.Database/GWAS/Sumstats/', disease,'.sumstats.gz'),header=T) %>%
  filter(SNP %in% rs$V2[-10])
ggplot(data=df,aes(x=SNP, y=log10(P)))+
  geom_bar(stat='identity',position='dodge')+
  theme_bw()

ggsave(file='all_leadVar.PD_gwas.pdf')

### plot PDUI and expression of the most significant peaks 
wkdir='/data/slurm/nanjh/scAPA_project/02.bwfiles/disease_example/Parkinson_Disease/ENST00000618500.4/'
setwd(wkdir)
peak.sig<- c('SNCA_Peak1/(Peak4+Peak1)','SNCA_Peak3/(Peak4+Peak3)')
df.gen<- fread('ENST00000618500.4.ped', header=F)[,c(2,7,8)] %>%
  mutate(genotype=paste0(V7,V8), ID=as.character(V2)) %>%
  select(c('ID','genotype'))

Merge.sig<- data.frame()
for(ct in c('Exc','Inh','Oli')){
    x<- fread(paste0(ct, '/SNCA.multi_peaks.ratio.bed.gz'), header=T, check.names = F) %>%
      filter(Gene %in% peak.sig)
    for(p in peak.sig){
      x1<- data.frame(celltype=ct, Peak=p,ID=names(x)[-c(1:4)],PDUI=t(x[which(Gene %in% p),-c(1:4)]))
      x1<- merge(df.gen,x1)
      Merge.sig<- rbind(Merge.sig,x1)
      ####plot the PDUI
      model<- lm( x1$PDUI~x1$genotype)
      signif<- format(signif(glance(summary(model))[[1,5]],3),scientific = TRUE)
      ggplot(data=x1, aes(x=genotype, y=as.numeric(PDUI),fill=genotype))+
        geom_violin(trim = FALSE,outlier.shape = NA, alpha=0.65)+
        #geom_point()+
        geom_boxplot(width=0.2, fill = "white", color = "black",outlier.shape = NA)+
        ggtitle(paste0('P value =',signif))+
        theme_bw()+
        ylab(paste0(ct,p, sep=':'))
      ggsave(paste0(ct,'.',gsub('/','_',p),'.boxplot.pdf'), width=4,height=4)
    }
}

saveRDS(Merge.sig, file='231201.all_cts.SNCA_peaks.rds')

Merge.sig<- data.frame()
for(ct in c('Exc','Inh','Oli')){
  x<- fread(paste0('/data/slurm/nanjh/scAPA_project/04.expression_matrix/re-QC/02.',ct,'_expression_QC.voom.bed.gz'), header=T, check.names = F) %>%
    filter(Symbol %in% 'SNCA')
  x1<- data.frame(celltype=ct,ID=names(x)[-c(1:4)],Expression=t(x[which(Symbol %in% 'SNCA'),-c(1:4)]))
  x1<- merge(df.gen,x1)
  Merge.sig<- rbind(Merge.sig,x1)
  model<- lm(x1$Expression~x1$genotype)
  signif<- format(signif(glance(summary(model))[[1,5]],3),scientific = TRUE)
  ggplot(data=x1, aes(x=genotype, y=as.numeric(Expression),fill=genotype))+
    geom_violin(trim = FALSE,outlier.shape = NA, alpha=0.65)+
    #geom_point()+
    geom_boxplot(width=0.2, fill = "white", color = "black",outlier.shape = NA)+
    ggtitle(paste0('P value =',signif))+
    theme_bw()+
    ylab(ct)
  ggsave(paste0(ct,'.','expression.boxplot.pdf'), width=4,height=4)
}
saveRDS(Merge.sig, file='231201.all_cts.SNCA_expression.rds')
