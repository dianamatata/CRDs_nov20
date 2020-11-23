#!/bin/bash
DIR_script=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/B_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
INDIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_gene_CRDs/merged
OUTDIR=/home/users/a/avalosma/scratch/6_CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_gene_CRDs/significants
mkdir -p $OUTDIR

for file in $INDIR/*.txt.gz; do
	filename=$(echo $file | rev | cut -d "/" -f1 | rev)
	Rscript $DIR_script/6bis_runFDR_cis.R $INDIR/$filename 0.05 $OUTDIR/$filename
done
