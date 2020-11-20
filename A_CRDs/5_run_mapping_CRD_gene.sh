#!/bin/bash

K=100
DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene
mkdir -p $OUTDIR $OUTDIR/mapping_aCRD_gene $OUTDIR/inverse_mapping_aCRD_gene $OUTDIR/mapping_aCRD_gene_nominal $OUTDIR/inverse_mapping_aCRD_gene_nominal

# also add sCRD-genes in the pipeline
mkdir -p $OUTDIR $OUTDIR/mapping_sCRD_gene $OUTDIR/inverse_mapping_sCRD_gene $OUTDIR/mapping_sCRD_gene_nominal $OUTDIR/inverse_mapping_sCRD_gene_nominal

# look for equivalent names of cell types from chip seq and methyl marks to RNA data
declare -A cell_types
cell_types[EGAD00001002670]="EGAD00001002675"
cell_types[EGAD00001002672]="EGAD00001002674"
cell_types[EGAD00001002673]="EGAD00001002671"

cell_types[methyl_neut]="EGAD00001002675"
cell_types[methyl_mono]="EGAD00001002674"
cell_types[methyl_tcell]="EGAD00001002671"

PC=10  #PC depends on cell type: change that!!!!

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
	RNA=$DATADIR/${cell_types[${cell_type}]}_RNA.PC$PC\.bed.gz
	MOD=$DATACRD/${cell_type}.ALLchr.mean.txt.gz

	for k in $(seq 1 $K); do
		OUT1=$OUTDIR/mapping_aCRD_gene/${cell_type}_mapping_aCRD_gene_mean_chunk$k
		if [ ! -f "$OUT1" ]
		then
			cmd1="QTLtools cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out ${OUT1}.txt"
			eval $cmd1
                	OUT2=$OUTDIR/inverse_mapping_aCRD_gene/${cell_type}_inverse_mapping_aCRD_gene_mean_chunk$k
                	cmd2="QTLtools cis --vcf $RNA --bed $MOD --permute 200 --chunk $k $K --out ${OUT2}.txt"
			eval $cmd2
		fi
	done
done


#<<'COMMENTS'

# the one used at the end

DATADIR0=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
	TH=$DATADIR0/${cell_type}_CLOMICS_v3.0/mapping_aCRD_gene/gene_CRD_mean_permutations.thresholds.txt
	RNA=$DATADIR/${cell_types[${cell_type}]}_RNA.PC$PC\.bed.gz
        MOD=$DATACRD/${cell_type}.ALLchr.mean.txt.gz
	for k in $(seq 1 $K); do
		OUTFILE1=$OUTDIR/mapping_aCRD_gene_nominal/${cell_type}_mapping_aCRD_gene_mean_nominal_chunk$k
        	cmd_nominal1="QTLtools cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUTFILE1.txt"
		eval $cmd_nominal1
                OUTFILE2=$OUTDIR/inverse_mapping_aCRD_gene_nominal/${cell_type}_inverse_mapping_aCRD_gene_mean_nominal_chunk$k
                cmd_nominal2="QTLtools cis --vcf $RNA --bed $BED --nominal $TH --chunk $k $K --out $OUTFILE2.txt"
                #eval $cmd_nominal2
	done
done

##COMMENTS

