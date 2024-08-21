library(ggplot2)
library(stringr)
library(data.table)
library(mice)
cts.list=c('Exc','Inh','Ast','Oli','Opc','Mic')

##Merge results of all chromosomes 
for (ct in cts.list){
  mat<- fread( paste0(matdir,"Dapars2_",ct,"_result_temp.chr1",".txt"),header=T) # see demoData for the format
  for ( i in c(seq(2,22),"X","Y")){
    mat1<-fread(paste0(matdir,"Dapars2_",ct,"_result_temp.chr",i,".txt"),header=T)
    mat<- rbind(mat, mat1)
  }
  write.table(mat, file=paste0(matdir,"Dapars2_",ct,"_result_temp.chrall.csv"),
              sep=",",quote=F, row.names=F)
}


for (ct in cts.list){
  apa<- fread(paste0(matdir,"Dapars2_",ct,"_result_temp.chrall.csv"),
              sep=",",header=T)
  nam<-gsub(paste0('.',ct,".possorted_genome_bam_PDUI"),"", gsub('wig\\/',"",names(apa)))
  names(apa)<- nam
  # write the APA matrix 
  write.csv(apa, file=paste0(matdir,"01Dapars2_",ct,"_result_temp.chrall.csv"))
  ###4.summarize the NA number per samples or per transcripts
  ###step1 summarize the NA number per gene
  
  sum_gen<-as.data.frame(matrix(0,0,8))
  names(sum_gen)<-c(names(summary(as.numeric(unlist(apa[1,-(1:4)])))),
                    "missing_rate")
  pdf("02.APA_pergene_miss_density.pdf", height=5, width=6,useDingbats=F)
  for (i in seq(1, nrow(apa))){
    sum1<- as.data.frame(t(as.matrix(summary(as.numeric(unlist(apa[i,-(1:4)]))))))
    if(ncol(sum1)==7){
      sum2<-cbind(sum1,sum1[,7]/(ncol(apa)-4))
      names(sum2)[8]<-"missing_rate"
      sum_gen<- rbind(sum_gen,sum2 )
    }else if(ncol(sum1)==6){
      sum2<-cbind(sum1,0,0)
      names(sum2)[c(7,8)]<-c("NA's","missing_rate")
      sum_gen<- rbind(sum_gen,sum2 )
    }
  }
  sum_gen<- as.data.frame(cbind(apa[,1:4],sum_gen))
  write.csv(sum_gen, 
            file=paste0("02",ct,"_APA_pergene_sum.csv"))
  p= ggplot(data=sum_gen)+
    geom_density(aes(x=missing_rate))+
    theme_bw()+
    ggtitle(ct)
  print(p)
  
}
dev.off()
##visualization of all pseudobulk cell-types
sum<- fread(paste0("02Exc","_APA_pergene_sum.csv"),sep=",",header=T)
sum1<- cbind("Exc",sum[,-1])
names(sum1)[1]<-"Cell-type"
for (ct in cts.list[-1]){
  sum2<- fread(paste0("02",ct,"_APA_pergene_sum.csv"),sep=",",header=T)
  sum3<- cbind(ct,sum2[,-1])
  names(sum3)[1]<-"Cell-type"
  sum1<- rbind(sum1,sum3)
}
pdf("02all_pseudo_celltype_misssing.pdf",height=5, width=6)
p<-ggplot(data=sum1)+
  geom_density(aes(x=missing_rate, fill=`Cell-type`),alpha=0.4)+
  geom_vline(xintercept=0.5,linetype="dashed")+
  theme_bw()
print(p)
dev.off()
write.table(sum1, file="02APA_qc_sum_allcelltypes.res",sep="\t",
            quote=F,col.names=F,row.names=F )
#5. filter the transcript whose missing rate higher than 0.5
sum_qc<- subset(sum1, missing_rate<=0.5)
sum_qc1<- as.data.frame(table(sum_qc$`Cell-type`))
names(sum_qc1)<-"Cell-types"
write.csv(sum_qc1,"03apa_qc_sum_percelltype.csv")
#6. extract the APA with high calling rate
for (ct in cts.list){
  apa<- fread(paste0(matdir,"01Dapars2_",ct,"_result_temp.chrall.csv"),sep=",",header=T)
  gene_qc<- subset(sum_qc, `Cell-type`==ct)
  apa1<- subset(apa, Gene %in% gene_qc$Gene)
  #7. summarize the missing rate per samples
  for (i in seq(6,ncol(apa1))){
    sum<- as.matrix(summary(as.numeric(unlist(
      subset(apa1, select=names(apa1)[i])[,1]
    ))))
    if( nrow(sum)==7){
      sum1<- cbind(names(apa1)[i],ct,sum[[7]],sum[[7]]/nrow(gene_qc))
    }else if(nrow(sum)==6){
      sum1<- cbind(names(apa1)[i],ct,0,0)
    }
    write.table(sum1, file="04APA_qc_samples_sum.res", append=T,sep="\t",
                quote=F,col.names=F,row.names=F)
  }
}
#8. visualization of missing rate of per samples across different cell-types
sum<- fread("04APA_qc_samples_sum.res",header=F)
names(sum)<-c("library_id","Cell-types","missing-number","missing-rate")
#9. filter the samples whose missing rate higher than 0.5
sum_qc<- subset(sum, `missing-rate`<=0.5)
sum_qc1<- as.data.frame(table(sum_qc$`Cell-type`))
names(sum_qc1)<-"Cell-types"
write.csv(sum_qc1,"05samples_qc_sum_percelltype.csv")
#10.extract the samples whose missing rate higher than 0.5
gene_qc<-fread("02APA_qc_sum_allcelltypes.res",header=F)
sample_qc<- fread("04APA_qc_samples_sum.res",header=F)
for (ct in cts.list){
  gene_qc1<- subset(gene_qc,V1==ct&V13<=0.5,select=V2)
  sample_qc1<- subset(sample_qc,V2==ct&V4<=0.5,select=V1)
  apa<- fread(paste0(matdir,"01Dapars2_",ct,"_result_temp.chrall.csv"),sep=",",header=T)
  apa1<- subset(apa, Gene %in% unlist(gene_qc1), select=c(names(apa)[2:5],unlist(sample_qc1$V1)))
  write.csv(apa1, paste0("06",ct,"_apa_sampleqc_matrix.csv"))
}

###using multiple imputation (MI) to imputate APA matrix

for (ct in cts.list){
  apa<- read.table(paste0("06",ct,"_apa_sampleqc_matrix.csv"),sep=",",header=T)
  apa1<- mice(apa[,-c(1:5)], seed=1234)
  completeData.2 <- cbind(apa[,2:5],mice::complete(apa1,1))
  write.csv(completeData.2,
            file=paste0("07",ct,"_imputation_MI_",i,".csv"))
  
}
