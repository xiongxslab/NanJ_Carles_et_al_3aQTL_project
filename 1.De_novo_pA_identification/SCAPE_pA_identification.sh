## The script was used to carry out de novo pA site identification using SCAPE 

# 1. gtf from 10X; /home/Genomes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf
conda activate scape
module load bedtools/2.26.0  # required
# 2. prepare gtf
sed -e '1,5d' genes.gtf | bedtools sort -i /dev/stdin | bgzip >  genes.sorted.gtf.gz
python ~/software/SCAPE/main.py prepare --gtf genes.sorted.gtf.gz --prefix hg38.10x

# 3. run scape for each pseudobulk 
while read samp
do
python ~/software/SCAPE/main.py apamix --bed hg38.10x_utr.bed --bam ${samp}.possorted_genome_bam.bam --out ${samp}.APA --cores 8 --cb ${samp}.barcodes.tsv.gz
done < samp.txt   

# 4. get grouped pA sites
python ~/software/SCAPE/scripts/group_pa.py --files $file --labels $labels --outfile collapse_pa.tsv.gz
# where:
# $file lists all the output ${samp}.APA/pasite.csv.gz files, separated by ","
# $labels lists the labels corresponding to each sample id

# 5. make pA x cell matrix
Rscript Make_pA_matrix.R

# 6. calculate 3' UTR length for each gene in each cell
Rscript Calculate.UTR3_length.R


