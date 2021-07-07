#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5_run_mapping_CRD_gene_permutation_pass.sh"
wsbatch -J 5p200_sh.job --partition=mono-shared-EL7 --time=04:00:00 -o 5p200_sh.out -e 5p200_sh.err --wrap="$cmd"



cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5_run_mapping_CRD_gene_nominal_pass.sh"
wsbatch -J 5nominal_sh.job --partition=mono-shared-EL7 --time=04:00:00 -o 5nominal_sh.out -e 5nominal_sh.err --wrap="$cmd"


cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5_run_mapping_CRD_gene_permutation_pass.sh"
wsbatch -J 5p1000_sh.job --partition=shared-bigmem-EL7 --time=04:00:00 -o 5p1000_sh.out -e 5p1000_sh.err --wrap="$cmd"



cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/5_run_mapping_CRD_gene_nominal_pass.sh"
wsbatch -J 5nom.job --partition=mono-shared-EL7 --time=11:00:00 -o 5nom.out -e 5nom.err --wrap="$cmd"


cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/10_correlation_peaks_vs_PCHiC.R "
wsbatch -J 10_1.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o 10_1.out -e 10_1.err --wrap="$cmd"


cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R"
wsbatch -J CRDQTL4.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o CRDQTL4.out -e CRDQTL4.err --wrap="$cmd"
