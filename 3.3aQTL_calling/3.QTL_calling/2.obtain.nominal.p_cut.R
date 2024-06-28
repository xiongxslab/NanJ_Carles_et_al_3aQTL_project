library(qvalue)

ct.list = c('Ast','Exc','Inh','Mic','Oli','Opc')

for(ct in ct.list){
        df = read.table(paste0(ct,'.permute.out.fdr_qvalue.txt'),head=T)
        df.sign = df[df$pval_perm_qvalue <= 0.05,]
        bpval.cut = max(df.sign$pval_beta)
        df.sign$npval.cut = qbeta(bpval.cut,df.sign$beta_shape1,df.sign$beta_shape2)
        write.table(df.sign,paste0(ct,".permute.out.qvalue.signif.txt"),quote=F,sep="\t",row.names = F)
}




