#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
CTM=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/trans/signif/centromeres.range.txt
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans
mkdir -p $OUTDIR $OUTDIR/OUT

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do                        
			MOD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
			for c0 in $(seq 1 22); do
				for c1 in $(seq $(($c0 + 1)) 22); do
					CTM0=$(cat $CTM | awk -v c0=$c0 '{ if ($1 == "chr"c0) print $2, $3 }')
					CTM1=$(cat $CTM | awk -v c1=$c1 '{ if ($1 == "chr"c1) print $2, $3 }')
					OUT=$OUTDIR/${data_type}_${cell_type}_${module}_chr_${c0}_${c1}_trans
					cmd="$CLOMICs trans --bed $MOD --regions $c0 $c1 --out ${OUT}.txt.gz --centromere $CTM0 $CTM1 --full"
					echo $cmd
					eval $cmd
					#one cmd is around 5sec. 2*3*2*22*22*4/(60*24)=16 days 
				done
			done
		done
	done
done

# verification 
cd $OUTDIR
count=$(ls | cut -d '_' -f1-3 | uniq -c)
count2=0
for c0 in $(seq 1 22); do
	for c1 in $(seq $(($c0 + 1)) 22); do
		count2=$((count2+1))	
	done
done

echo "$count" 
echo "expected: $count2" # somme de 1 a n= n*(n+1)/2 with n=21, result=231

#HIC_VAL=/data/unige/funpopgen/odelanea/SGX/V2/data_hic/chr$c0\_chr$c1\/MAPQG0/chr$c0\_$c1\_5kb.RAWobserved.gz
#HIC_KR0=/data/unige/funpopgen/odelanea/SGX/V2/data_hic/chr$c0\_chr$c1\/MAPQG0/chr$c0\_5kb.KRnorm.gz
#HIC_KR1=/data/unige/funpopgen/odelanea/SGX/V2/data_hic/chr$c0\_chr$c1\/MAPQG0/chr$c1\_5kb.KRnorm.gz
#LOUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/TRANS/EGAD00001002670_ALL.$c0\.$c1
