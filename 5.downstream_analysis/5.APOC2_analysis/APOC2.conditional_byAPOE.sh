plink \
--bfile 00.reference_pannel/EUR \
--chr 19 \
--maf 0.01 \
--cojo-file Alzheimer_disease.summstats.ma \
--cojo-cond 01.APOE_conditonal.snp \
--out 02.APOE_conditonal_chr19

