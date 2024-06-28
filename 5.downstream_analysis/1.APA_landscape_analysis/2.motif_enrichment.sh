#!/bin/bash
##miRNA and RBP motif datasets were downloaded from MEME.
for ct in Ast Exc Opc Oli Inh Mic Neuron
do
  for tp in lengthening shortening
  do
    sea --p 03${ct}_${tp}.RNA.fasta --m ${miRNA_database_dir}/Homo_sapiens_hsa.meme --o ${ct}_${tp}_miRNA
    sea --p 03${ct}_${tp}.RNA.fasta --m ${RBP_database_dir}/Ray2013_rbp_Homo_sapiens.meme --o ${ct}_${tp}_RBP

  done
done

