sample_file=$1
wkdir=$2
pip_dir='/data/slurm/nanjh/scAPA_project/34.NG.revision2/1.regulator.validation/0.pipline/'
cd $wkdir
for sm in $(cat $sample_file)
do
	echo "bash "${pip_dir}1.fastaQC.basic.sh ${sm}"&& bash "${pip_dir}"2.mapping.basic.sh "${sm}" && bash ${pip_dir}3.stringtie.basic.sh" ${sm} >> 1.1.fastqc_mapping_TPM.submit.sh
done
