#!/bin/bash
genome_path="/data/slurm/dingk/07.tRNA-QTL/00reference/hg38"
sample_name=$1
bamdir="./03.bam"
mkdir -p $bamdir 
cd $bamdir

  hisat2 \
  -p 30 \
        -x ${genome_path} \
        -1  ../02.clean/${sample_name}'_R1_val_1.fq.gz' -2 ../02.clean/${sample_name}'_R2_val_2.fq.gz' \
        -S ${sample_name}.sam \
        --no-mixed \
        --no-discordant \
        --no-unal \
        --time \
        --omit-sec-seq \
  2>&1 | tee -a ../log/${sample_name}.log
#  echoLog hisat2 Done
  #rm $read1 $read2
  #echoLog sam to bam
  samtools view \
  -b \
  -S \
  -@ 30 \
  -o ${sample_name}.bam \
  ${sample_name}.sam \
  2>&1 | tee -a ../log/${sample_name}.log

 # echoLog sort bam by name
  samtools \
    sort \
    -n \
    -@ 30 \
    -o ${sample_name}.nsort.bam \
    ${sample_name}.bam \
    2>&1 | tee -a ../log/${sample_name}.log
  #  echoLog mark mate tag
  samtools \
    fixmate \
    -@ 30\
    -m ${sample_name}.nsort.bam \
    ${sample_name}.fixmate.bam \
    2>&1 | tee -a ../log/${sample_name}.log

 # echoLog sort bam by position
  samtools \
    sort \
    -@ 30 \
    -o ${sample_name}.sort.bam \
    ${sample_name}.fixmate.bam \
    2>&1 | tee -a ../log/${sample_name}.log
 #summary the mapping rate
 total_reads=$(samtools view -c ${sample_name}.sort.bam)
 mapped_reads=$(samtools view -c -F 4 ${sample_name}.sort.bam)
 mapping_rate=$(echo "scale=2; ($mapped_reads / $total_reads) * 100" | bc)
 echo $sample_name "Total Reads: $total_reads">> ${sample_name}.mappingRatio
 echo $sample_name "Mapped Reads: $mapped_reads">> ${sample_name}.mappingRatio
 echo $sample_name "Mapping Rate: $mapping_rate%" > ${sample_name}.mappingRatio
#samtools flagstat ${sample_name}.sort.bam > ${sample_name}.mappingRate
 # echoLog mark and remove duplication
  samtools \
    markdup \
    -r \
    -s \
    -@ 30 \
    ${sample_name}.sort.bam \
    ${sample_name}.rmdup.bam \
    2>&1 | tee -a ../log/${sample_name}.log
  #  echoLog index bam
  samtools \
    index \
    -@ 30 \
    ${sample_name}.rmdup.bam \
    2>&1 | tee -a ../log/${sample_name}.log

  rm \
    ${sample_name}.sam \
    ${sample_name}.bam \
    ${sample_name}.nsort.bam \
    ${sample_name}.fixmate.bam \
    ${sample_name}.sort.bam
