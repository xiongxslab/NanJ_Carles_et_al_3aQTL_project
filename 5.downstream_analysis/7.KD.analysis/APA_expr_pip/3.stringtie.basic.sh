#!/bin/bash
sample_name=$1
gff_file='/data/slurm/nanjh/scAPA_project/00.ref_bed/gencode.v32.annotation.gtf'
mkdir -p 04.count
cd 03.bam
#samtools sort -@ 50 -n ${sample_name}.rmdup.bam -o ${sample_name}.nsort.bam
#samtools view -h ${sample_name}.nsort.bam >${sample_name}.nsort.sam

stringtie -G $gff_file -e -p 20 -o ../04.count/${sample_name}.stringtie.gtf -A ../04.count/${sample_name}.stringtie.tab ${sample_name}.rmdup.bam
