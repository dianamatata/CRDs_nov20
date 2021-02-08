#!/bin/bash
permuts=1000
K=100

VCF=/home/users/a/avalosma/scratch/10_CRD_QTLs/All_chr.BPWP10_13_12_15.vcf.gz
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/12_TRIPLETS

mkdir -p $OUTDIR $OUTDIR/mapping_$permuts $OUTDIR/merged_$permuts

declare -A rna_file
rna_file[neut]="EGAD00001002675"
rna_file[mono]="EGAD00001002674"
rna_file[tcell]="EGAD00001002671"

PC=10
#PC depends on cell type: change that!!!!

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do
                        BED=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
                        for k in $(seq 1 $K); do
                                OUT1=$OUTDIR/mapping_$permuts/${data_type}_${cell_type}_${module}CRD_gene_var_chunk$k
                                cmd1="QTLtools cis --vcf $VCF --bed $BED --permute $permuts --chunk $k $K --out ${OUT1}.txt"
				wsbatch -J 1000permuts.job --partition=shared-cpu --time=04:00:00 -o 1000permuts.out -e 1000permuts.err --wrap="$cmd1"
                        done
                done
        done
done


#<< 'MULTILINE-COMMENT'

### concatenate at the end

for module in 'mean' 'loom' ; do
        for condition in 'mapping_CRD_gene' ; do
                for data_type in  'methyl' 'hist' ; do
                        for cell_type in 'neut' 'mono' 'tcell' ; do
                                name=${data_type}_${cell_type}_${module}
                                cat $OUTDIR/mapping_$permuts/${name}_*.txt | gzip -c > $OUTDIR/merged_$permuts/${name}_${condition}_ALL.txt.gz
                        done
                done
        done
done

# MULTILINE-COMMENT
