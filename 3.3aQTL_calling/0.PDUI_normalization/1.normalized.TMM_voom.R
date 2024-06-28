library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(qvalue)

ct.list = c('Ast','Exc','Inh','Mic','Oli','Opc')

library(stringr)
for(ct in ct.list){
	print(ct)
	mat = readRDS(paste0('../0.matrix/08',ct,'_imputation_MI_3_naomit.parse.rds'))
	info = mat[,1:4]
	mat = mat[,-c(1:4)]
	write.table(mat,'tmp.txt',quote=F,sep="\t") # a trick to make the matrix numeric
	mat = read.table('tmp.txt',head=T) 
	system('rm tmp.txt')
	y <- DGEList(counts=mat)
	## perform M-value normalization
	y <- calcNormFactors(y)
	y = voom(y)
	y.mtx = y$E

	all(rownames(y.mtx) == rownames(info))
	mat.t = cbind(info,y.mtx)
	colnames(mat.t) = gsub("X",'',colnames(mat.t))
	write.table(mat.t,paste0('./08',ct,'_imputation_MI_3_naomit.parse.voom.txt'),quote=F,sep="\t",row.names=F)
}




