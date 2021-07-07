#!/bin/bash
permuts=1000
K=100

DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/5_CRDgene

mkdir -p $OUTDIR $OUTDIR/mapping_CRD_gene_permuts_$permuts $OUTDIR/merged_$permuts

declare -A rna_file 
rna_file[neut]="EGAD00001002675"
rna_file[mono]="EGAD00001002674"
rna_file[tcell]="EGAD00001002671"
 
PC=10  
#PC depends on cell type: change that!!!!

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
 
	                RNA=$DATADIR/${rna_file[${cell_type}]}_RNA.PC$PC\.bed.gz
        	        MOD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
 
                	for k in $(seq 1 $K); do
                        	OUT1=$OUTDIR/mapping_CRD_gene_permuts_$permuts/${data_type}_${cell_type}_${module}_CRD_gene_chunk$k
                         	cmd1="QTLtools cis --vcf $MOD --bed $RNA --permute $permuts --chunk $k $K --out ${OUT1}.txt"
                              	eval $cmd1
			done
                done
        done
done




### concatenate at the end
mkdir -p $OUTDIR/merged

                
for module in 'mean' 'loom' ; do                        
	for condition in 'mapping_CRD_gene' ; do
		for data_type in  'methyl' 'hist' ; do
			for cell_type in 'neut' 'mono' 'tcell' ; do
				name=${data_type}_${cell_type}_${module}
				cat $OUTDIR/mapping_CRD_gene_permuts_$permuts/${name}_*.txt | gzip -c > $OUTDIR/merged_$permuts/${name}_${condition}_ALL.txt.gz
			done
		done
	done
done
# COMMENTS

