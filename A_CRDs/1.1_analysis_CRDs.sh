#!/bin/bash

OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD
ANALYSIS_FILE=$OUT_FOLDER/analysis_CRD_counts.txt

echo > $ANALYSIS_FILE
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
	echo -e $cell_type' \t '$(zcat $OUT_FOLDER/${cell_type}_ALL.CRDs.txt.gz | wc -l) >> $ANALYSIS_FILE
done 
