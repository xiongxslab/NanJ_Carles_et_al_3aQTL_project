file1=$1
file2=$2
sm_new=$3
wkdir=$4

cd $wkdir
zcat ${file1}._R1.fastq.gz ${file2}_R1.fastq.gz >> 01.rawdata/${sm_new}_R1.fastq.gz
zcat ${file1}._R2.fastq.gz ${file2}_R2.fastq.gz >> 01.rawdata/${sm_new}_R2.fastq.gz

echo $sm_new >> sm.merged.list
