library(ggplot2)

samp.list = c('Ast','Exc','Inh','Mic','Oli','Opc')
pf.list = c(1:20)
start.list = c(1,2)

sum.df = data.frame()
for(samp in samp.list){
 for(start in start.list){
  for(pf in pf.list){
   if(start>pf){next}
   df.merge = data.frame()
   for(chr in 1:22){
    out = read.table(paste0('output.1M/',samp,'.chr',chr,'.pf',start,'-',pf,'.out.txt.gz'),head=F)
    if(ncol(out)==1){next}
    df.merge = rbind(df.merge,out)
    print(chr)
   }
   sum.df = rbind(sum.df,data.frame(samp=samp,pf_start=start,pf=pf,p_cut=0.0001,nqtl=nrow(df.merge),n_gTRNA=nrow(df.merge[!duplicated(df.merge$V1),]),min_P=min(df.merge$V7)))
   df.merge = df.merge[df.merge$V7<0.00001,]
   sum.df = rbind(sum.df,data.frame(samp=samp,pf_start=start,pf=pf,p_cut=0.00001,nqtl=nrow(df.merge),n_gTRNA=nrow(df.merge[!duplicated(df.merge$V1),]),min_P=min(df.merge$V7)))
  }
 }
}

saveRDS(sum.df,'APA_QTL.voom.20230625.fastqtl_test.rds')

pdf('APA_QTL.voom.20230625.fastqtl_test.pdf',width=10,height=5,useDingbats=F)

for(cutoff in c(1e-04,1e-05)){
sum.df.cutoff =  sum.df[sum.df$p_cut == cutoff,]
p = ggplot(sum.df.cutoff,aes(x=pf,y=nqtl,color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(samp~.,scales="free_y") + ggtitle(paste0('# QTL; p.cut=',cutoff))
print(p)

p = ggplot(sum.df.cutoff,aes(x=pf,y=n_gTRNA,color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(samp~.,scales="free_y") + ggtitle(paste0('# g-APA; p.cut=',cutoff))
print(p)

p = ggplot(sum.df.cutoff,aes(x=pf,y=-log10(min_P),color=as.character(pf_start))) + geom_point() + geom_line() + theme_bw() +
  facet_wrap(samp~.,scales="free_y") + ggtitle('Min p-value')
print(p)

}
dev.off()






