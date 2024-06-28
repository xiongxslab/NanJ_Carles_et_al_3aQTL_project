## calculate 3'UTR length for each gene in each cell
samp.list = read.table('samp.txt',head=F)
samp.list = as.character(samp.list$V1)

if(0){   # run once
# read in APA annotation
anno_loc = read.table('All.pasites.genome_loc.tsv',head=F)
anno = read.table('All.pasites.annotation.tsv',head=F)
anno_loc$idx = paste0(anno_loc$V1,":",anno_loc$V2)
anno$idx = paste0(anno$V6,":",anno$V2)
anno = merge(anno_loc,anno,by="idx",suffix=c(".gppy",".scape"))
# filter for 3'UTR
anno = anno[anno$V10 %in% c('3UTRs','3UTRs_1kb','3UTRs_2kb','LastExon1Kb') & anno$V4.gppy %in% c('downstream','exon'),]

tx_info = read.table('/home/xiongxs/data/Proj_mix/scAPA/3.SCAPE/Ref/hg38.10x.genes.info.txt',head=T,sep="\t")
# consider coding-genes only
tx_info = tx_info[tx_info$gene_type == "protein_coding",]
tx_info.cut = tx_info[,c('tx_name','utr5_len','cds_len','utr3_len')]

anno.info = merge(anno,tx_info.cut,by.x='V1.gppy',by.y='tx_name')
anno.info.utr = anno.info[anno.info$V4.gppy == "exon",]
anno.info.dstr = anno.info[anno.info$V4.gppy == "downstream",]
anno.info.utr$pA_UTR3_len = anno.info.utr$V3.gppy - anno.info.utr$utr5_len - anno.info.utr$cds_len
anno.info.dstr$pA_UTR3_len = anno.info.dstr$utr3_len + anno.info.dstr$V3.gppy
anno.info = rbind(anno.info.utr,anno.info.dstr)
anno.info.cut = anno.info[,c('V1.gppy','idx','V3.gppy','V4.gppy','V12','V8','utr5_len','cds_len','utr3_len','pA_UTR3_len')]
colnames(anno.info.cut)[1:6] = c('TxName','Tx.pa_idx','Tx.pos','pa_Type','Genome.pa_idx','GeneName')
rownames(anno.info.cut) = anno.info.cut$Genome.pa_idx
write.table(anno.info.cut,'All.pasites.annotation.pa_UTR_length.tsv',quote=F,sep="\t",row.names=F)
}

library(argparser, quietly=TRUE)
p <- arg_parser("sample for pa 3'UTR length;")
p <- add_argument(p, "samp", help="samp for analysis")

arg <- parse_args(p)
samp=arg$samp

anno.info.cut = read.table('All.pasites.annotation.pa_UTR_length.tsv',head=T)
rownames(anno.info.cut) = anno.info.cut$Genome.pa_idx
library(stringr)

pa = read.table(paste0(samp,'.APA/pasite.csv.gz'),head=T,sep=",")
pa.name = data.frame(str_split_fixed(pa[,1],':',4))
rownames(pa) = paste0(pa.name[,1],':',pa.name[,2],':',pa.name[,4])
pa = pa[,-1]
pa.name$idx = paste0(pa.name[,1],':',pa.name[,2],':',pa.name[,4])
pa.name = pa.name[pa.name$idx %in% anno.info.cut$Genome.pa_idx,]
anno.info.cut.samp = anno.info.cut[as.character(pa.name$idx),]
pa = pa[as.character(anno.info.cut.samp$Genome.pa_idx),]
lengths = anno.info.cut.samp$pA_UTR3_len
# calculate weighted length for each gene in each cell
pa2 = apply(pa, 2, function(x) x*lengths)
pa.bygene.total = aggregate(pa2, by=list(as.character(anno.info.cut.samp$TxName)), FUN=sum)
pa.bygene.npa = aggregate(pa, by=list(as.character(anno.info.cut.samp$TxName)), FUN=sum)
all(rownames(pa.bygene.total) == rownames(pa.bygene.npa)) # check order
all(colnames(pa.bygene.total) == colnames(pa.bygene.npa))
pa.bygene.total.mat = as.matrix(pa.bygene.total[,-1])
pa.bygene.npa.mat = as.matrix(pa.bygene.npa[,-1])
pa.bygene.weighted.len = pa.bygene.total.mat/pa.bygene.npa.mat
rownames(pa.bygene.weighted.len) = pa.bygene.total[,1]
saveRDS(pa.bygene.weighted.len,paste0('pa_UTR3_avg_length/',samp,'.pa_UTR3_avg_length.bygene.rds'))

