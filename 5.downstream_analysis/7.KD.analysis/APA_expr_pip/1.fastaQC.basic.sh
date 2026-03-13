clean_data_dir=./02.clean/
log_dir=./log
#bowtie2_index=./1.reference/susScr11.fa
data_dir='01.rawdata/'
d=$1
# step1: Filtering
#fastqc -t 10 -o ./fastqc $data_dir/$d'_R1.fastq.gz' $data_dir/$d'_R2.fastq.gz'
#cd $data_dir
#mv $d'_S1_L001_R1_001.fastq.gz' $d'_R1.fastq.gz'
#mv $d'_S1_L001_R2_001.fastq.gz' $d'_R2.fastq.gz'
#cd ../
mkdir $clean_data_dir
mkdir $log_dir
mkdir 01.fastqc_trimmed 
#source /data1/jhnan/skin_color/chIP_test/cutadapt_venv/bin/activate
trim_galore -q 30 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o ./02.clean \
            $data_dir/$d'_R1.fastq.gz' \
            $data_dir/$d'_R2.fastq.gz'
fastqc -t 10 -o ./01.fastqc_trimmed ./02.clean/$d'_R1_val_1.fq.gz' ./02.clean/$d'_R2_val_2.fq.gz' &
printf "Step 1: Filtering finished at `eval date +%Y%m%d"_"%H:%M:%S`\n" >> $log_dir/$d'.calling.log'
