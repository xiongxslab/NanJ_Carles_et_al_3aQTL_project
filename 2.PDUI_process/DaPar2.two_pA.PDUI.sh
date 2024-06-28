## The script was used to carry out PDUI quantification for each gene, based on two-pA site model by DaPars2 at the (sample + cell type) pseudobulk level 

# 1. get the sequencing depth of mappable reads for each pseudobulk
for ct in Ast Exc Inh Mic Oli Opc Vas
do
perl Get.sequencing.depth.pl $ct
done 

# 2. make the bash scripts for generating Dapars2_configure_file for each cell type
perl Make_DaPars.script.pl
# this step will generate Dapars2_configure_file.$ct files that will be used for the next step

# 3. run DaPars2, this will generate a sample x APA matrix with PDUI values for each cell type
for ct in Ast Exc Inh Mic Oli Opc Vas
do
module swap python python3/3.6.4
python3 ~/software/DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py Dapars2_configure_file.$ct Ref/chrList.txt
done 





