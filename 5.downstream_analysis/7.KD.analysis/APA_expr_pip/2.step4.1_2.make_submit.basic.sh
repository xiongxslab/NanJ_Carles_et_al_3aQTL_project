sample_file=$1
wkdir=$2
pip_dir='/data/slurm/nanjh/scAPA_project/34.NG.revision2/1.regulator.validation/0.pipline/'
cd $wkdir
for sm in $(cat $sample_file)
do
	echo "bash "${pip_dir}"4.1.bam2wig.basic.sh "${sm}" && bash "${pip_dir}"4.2.Get.sequencing.depth.basic.sh "${sm} >> "1.2.big2wig_depth.submit.sh"
done
