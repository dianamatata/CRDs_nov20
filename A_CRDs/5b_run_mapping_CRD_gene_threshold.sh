#!/bin/bash

K=100
DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene

mkdir -p $OUTDIR $OUTDIR/mean $OUTDIR/loom 
for module in 'mean' 'loom' ; do
        mkdir -p $OUTDIR/${module}/mapping_CRD_gene \
	$OUTDIR/${module}/inverse_mapping_CRD_gene \
	$OUTDIR/${module}/mapping_CRD_gene_nominal \
	$OUTDIR/${module}/inverse_mapping_CRD_gene_nominal
done


declare -A rna_file 
rna_file[neut]="EGAD00001002675"
rna_file[mono]="EGAD00001002674"
rna_file[tcell]="EGAD00001002671"
 
PC=10  
#PC depends on cell type: change that!!!!

##<<'COMMENTS'
# TH files copied from Guillaume, see in folder $OUTDIR/permutations_thresholds the commands to find it
#  DATADIR0=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS
# where is TH defined? Need to run the non nominal (permut) version, then the FDR, then get the threshold

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
 
                        RNA=$DATADIR/${rna_file[${cell_type}]}_RNA.PC$PC\.bed.gz
                        MOD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
	                #TH=$OUTDIR/permutations_thresholds/${data_type}_${cell_type}_mean_perm_threshold.txt
 			TH=$OUTDIR/significants/FDR_0.05_${data_type}_${cell_type}_${module}_mapping_CRD_gene_permuts.thresholds.txt
##### ERROR: Cannot find thresholds for 99 phenotypes!			
                        for k in $(seq 1 $K); do
                                OUT1=$OUTDIR/${module}/mapping_CRD_gene_nominal/${data_type}_${cell_type}_${module}_CRD_gene_chunk$k
                                cmd1="QTLtools cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out ${OUT1}.txt"
                                eval $cmd1
                                OUT2=$OUTDIR/${module}/inverse_mapping_CRD_gene_nominal/${data_type}_${cell_type}_${module}_inverse_CRD_gene_chunk$k
                                cmd2="QTLtools cis --vcf $RNA --bed $MOD --nominal $TH --chunk $k $K --out ${OUT2}.txt"
                                eval $cmd2
                        done
                done
        done
done

##COMMENTS

### concatenate at the end
mkdir -p $OUTDIR/merged

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
			for condition in 'mapping_CRD_gene_nominal' 'inverse_mapping_CRD_gene_nominal'; do
				name=${data_type}_${cell_type}_${module}
                        	cat $OUTDIR/${module}/${condition}/${name}_*.txt | gzip -c > $OUTDIR/merged/${name}_${condition}_permuts.txt.gz
			done
                done
        done
done

# COMMENTS

