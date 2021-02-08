#!/bin/bash
DIR=/home/users/a/avalosma/scratch/5_CRDgene/significants_1000
analysis_file=4.1_sizes.txt
echo ' ' > $analysis_file

for file in $DIR/*"mean"*"significant"*.txt ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(cat $file | wc -l)
        echo $length $f >> $analysis_file
done;

DIR=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000
for file in $DIR/*"mean"*.txt.gz ; do
        f=$(echo $file | rev | cut -d "/" -f1 | rev)
        length=$(zcat $file | wc -l)
        echo $length $f >> $analysis_file
done;
