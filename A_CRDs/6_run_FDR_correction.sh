#!/bin/bash
INDIR=/home/users/a/avalosma/scratch/5_CRDgene/merged
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene/significants
mkdir -p $OUTDIR

for file in $INDIR/*.txt.gz; do
	filename=$(echo $file | rev | cut -d "/" -f1 | rev)
	Rscript 6bis_runFDR_cis.R $INDIR/$filename 0.05 $OUTDIR/$filename
done
