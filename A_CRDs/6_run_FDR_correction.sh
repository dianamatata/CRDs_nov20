#cat gene_CRD_mean_chunk*.txt | gzip -c > gene_CRD_mean_permutations_full.txt.gz
Rscript 6bis_runFDR_cis.R gene_CRD_mean_permutations_full.txt.gz 0.05 gene_CRD_mean_permutations
