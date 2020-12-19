#!/bin/bash

DIR=/home/users/a/avalosma/scratch/5_CRDgene/significants

analysis_file=signif_analysis.txt
echo ' ' > $analysis_file

for file in $DIR/*"0.05"*.txt ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file
done;

analysis_file=signif_analysis_mean.txt
echo ' ' > $analysis_file

for file in $DIR/*"0.05_hist"*"mean"*.txt ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file
done;


# compare with G's values
DIRG=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS

analysis_file2=signif_analysis_Guillaume.txt

echo ' ' > $analysis_file2

for file in $DIRG/*/mapping_aCRD_gene/gene_CRD_mean_permutations* ; do
        f=$(echo $file | rev | cut -d "/" -f1-3 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file2
done


for file in $DIRG/*/mapping_aCRD_gene_0.01/gene_CRD_mean_permutations* ; do
        f=$(echo $file | rev | cut -d "/" -f1-3 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file2
done

