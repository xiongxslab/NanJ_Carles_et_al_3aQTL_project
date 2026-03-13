#sample_file=$1
wkdir=$1
cd $wkdir
cat 04.depth/*|sed 's/ /\t/g' > mapping_wig_location_with_depth.txt
python3 /data/slurm/nanjh/scAPA-project/tools/DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file /data/slurm/nanjh/scAPA_project/00.ref_bed/chrList.txt
