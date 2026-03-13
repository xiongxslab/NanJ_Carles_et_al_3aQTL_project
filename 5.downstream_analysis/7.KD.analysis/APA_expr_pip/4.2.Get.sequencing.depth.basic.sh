#for id in $(cat sample.list)
#do
samp=$1
if [ -d 04.depth ];then
        mkdir -p 04.depth
fi
samtools view -F 4 03.bam/$samp.rmdup.bam | wc -l |awk '{print "04.bigwig/"v".wig",$1}' v=$samp > 04.depth/$samp.mapping_wig_location_with_depth.txt
#perl Get.sequencing.depth.pl $id
#done

