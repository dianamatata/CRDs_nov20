#!/bin/bash
DIR=/home/users/a/avalosma/scratch/5_CRDgene
DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
DATACRD=/home/users/a/avalosma/scratch/2_CRD
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
CLOMICs=/srv/beegfs/scratch/groups/funpopgen/Tools/clomics1/bin/clomics
OUT=/home/users/a/avalosma/scratch/12_TRIPLETS/PC1
DIR_2col=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000/2_columns

module='mean'
for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		for c in $(seq 1 22); do
                        LI=$DATADIR_0/quantif_M_${data_type}_${cell_type}.bed.gz
			LM=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
			LT=$DATADIR/${data_type}_${cell_type}.chr$c\.module.txt.gz
			LO=$OUT/quantify/${data_type}_${cell_type}.chr$c
			L2=$DIR_2col/${data_type}_${cell_type}_mean_mapping_CRD_gene_ALL.txt
			cmd_maris="$CLOMICs quantify --bed $LI $LM --grp $L2 --pca 1 --normal --out $LO.pc1.txt"
			eval $cmd_maris
		done
        done
done

#                       cmd1="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.pc1.txt.gz --pca 1 --normal"


<< 'MULTILINE-COMMENT'

awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' $LO.pc1.txt | bgzip -c > $LO.pc1.txt

tabix -p bed $LO.pc1.txt


MULTILINE-COMMENT

