#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/2_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT


# histone CRDs
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	LI=$DATADIR_0/${cell_type}_merged_residuals.bed.gz
        LM=$DATADIR/${cell_type}_ALL.CRDs.txt.gz
	for c in $(seq 1 22); do
		echo "$cell_type $c"
		LT=$DATADIR/$cell_type\.chr$c\.module.txt.gz
		LO=$OUT_FOLDER/$cell_type\.chr$c\

		cmd1="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.pc1.txt.gz --pca 1 --normal"
		cmd2="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.mean.txt.gz --mean --normal"
		cmd3="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.loom.txt.gz --loo 0 --normal"

                eval $cmd1; eval $cmd2; eval $cmd3
                JOB=${cell_type}_chr${c}_quantify
		#wsbatch -J $JOB\.job --partition=mono-EL7 --time=00:10:00 -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd && $cmd2 && $cmd3" 
	done
done


# methyl CRDs, only LI changes
for cell_type in 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
        LI=$DATADIR_0/${cell_type}.bed.gz
        LM=$DATADIR/${cell_type}_ALL.CRDs.txt.gz
        for c in $(seq 1 22); do
                echo "$cell_type $c"
                LT=$DATADIR/$cell_type\.chr$c\.module.txt.gz
                LO=$OUT_FOLDER/$cell_type\.chr$c\

                cmd1="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.pc1.txt.gz --pca 1 --normal"
                cmd2="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.mean.txt.gz --mean --normal"
                cmd3="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.loom.txt.gz --loo 0 --normal"

                eval $cmd1; eval $cmd2; eval $cmd3
                JOB=${cell_type}_chr${c}_quantify
                #wsbatch -J $JOB\.job --partition=mono-EL7 --time=00:10:00 -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd && $cmd2 && $cmd3"
        done
done
