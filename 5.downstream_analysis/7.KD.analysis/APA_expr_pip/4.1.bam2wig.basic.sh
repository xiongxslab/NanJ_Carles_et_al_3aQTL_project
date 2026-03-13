#for id in $(cat sample.list)
#do
id=$1
if [ -d 04.bigwig ];then
	mkdir -p 04.bigwig
fi
	bedtools genomecov -ibam 03.bam/${id}.rmdup.bam -bga -split -trackline > 04.bigwig/${id}.wig
#done
