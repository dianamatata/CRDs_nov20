#!/bin/bash

# for each gene, find it and extract the type and the name

inFolder1='/Users/dianaavalos/Programming/CRDs_nov20/A_CRDs/trans/pathways'
inFolder='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:TRHs/trans_hubs_genes'
queryFile='/home/users/a/avalosma/scratch/RNA_rpkm/gene_info_71_74_75.txt'

# query file obtained with "  sort EGAD00001002675_gene_info.txt EGAD00001002674_gene_info.txt EGAD00001002671_gene_info.txt | uniq > gene_info_71_74_75.txt " to avoid dupliactes in the merge

outFolder= inFolder

outFile=$outFolder/gene_name.txt
outFile2=$outFolder/gene_name_type.txt
''>$outFile
''>$outFile2

for file in "$inFolder"/*; do
        printf $file
	file_nbr=$(echo $file | awk -F '[_.]'  '{print $4"_"$5}')
	outFile=$outFolder/gene_name_${file_nbr}.txt
	outFile2=$outFolder/gene_name_type_${file_nbr}.txt
	'' > $outFile
	'' > $outFile2
	l=0

	while IFS= read -r line
	do

		if [ $l -gt 0 ] #remove header line
		then
  			gene=$(echo "$(echo $line | cut -f1)")
			a=$(cat $queryFile | grep $gene )
			if [ ${#a} -gt 1 ]
			then
				echo $a | awk -F '[;=]' 'BEGIN {OFS="\t"}; {print $4" "$8}' >> $outFile2
				echo $a | awk -F '[;=]' '{print $8}' >> $outFile
			fi
		fi

mkdir -p $outFolder	((l+=1))
	done <  $file

done


        while IFS= read -r line
        do
		li=$(echo $line)
        done <  $file


gene='ENSG00000204113.4'
# for loop on all the genes, grep $gene

a=$(cat sample_EGAD00001002675_gene_info.txt | grep 'ENSG00000155959.6')
# len(a)= ${#a}
if [ ${#a} -gt 1 ]
then
	echo $a | awk -F '[;=]' 'BEGIN {OFS="\t"}; {print $4" "$8}'
	echo $a | awk -F '[;=]' '{print $8}'
fi
