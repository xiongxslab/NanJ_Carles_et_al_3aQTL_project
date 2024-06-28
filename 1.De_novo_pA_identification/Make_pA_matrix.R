## get names
samps = read.table('samp.txt',head=F)
samps = as.character(samps$V1)
files = paste0(samps,".APA/pasite.csv.gz")
names(files) <- gsub(".APA","",basename(dirname(files)))

collapse_pa = './collapse_pa.tmp.tsv.gz'
# get matrix
source('~/software/SCAPE/SCAPE.R/R/loadData.R')

pa_mtx <- loadData(
  fileList = files,
  collapsePa = collapse_pa,
  matrix = TRUE,
  cores = 12
)

saveRDS(pa_mtx,'pa_mtx.rds')

