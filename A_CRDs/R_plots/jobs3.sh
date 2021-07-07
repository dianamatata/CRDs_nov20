for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/10_correlation_peaks_vs_PCHiC_2.R $data_type $cell_type"
		wsbatch -J peak.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o peak.out -e peak.err --wrap="$cmd"
	done
done



for file in /home/users/a/avalosma/scratch/11_GSEA/trans_hubs_genes/genes_* ; do
	echo $file
	cmd="genewalk --project p1_$file --genes $file --id_type ensembl_id"
	wsbatch -J gene.job --partition=shared-bigmem  --time=5:00:00 --mem-per-cpu=100000 -o gene.out -e gene.err --wrap="$cmd"
done

for file in /home/users/a/avalosma/scratch/11_GSEA/trans_hubs_genes/genes_* ; do
	echo $file
	cmd="genewalk --project p1_$file --genes $file --id_type ensembl_id"
	wsbatch -J gene1.job --partition=shared-cpu  --time=12:00:00 --mem-per-cpu=10000 -o gene1.out -e gene1.err --wrap="$cmd"
done