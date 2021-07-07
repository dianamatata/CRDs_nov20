
shared-bigmem 4 Days	250GB days-hours:minutes:seconds 12:00:00


Goals:

R script modulable avec les args and launch all the different params

wsbatch -J m1.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m1.out -e m1.err --wrap="$cmd"



cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean tcell neut"
wsbatch -J m2.job --partition=shared-bigmem --time=12:00:00 --mem-per-cpu=500000 -o m2.out -e m2.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean tcell mono"
wsbatch -J m3.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m3.out -e m3.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean mono neut"
wsbatch -J m4.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m4.out -e m4.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean mono tcell"
wsbatch -J m5.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m5.out -e m5.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean neut tcell"
wsbatch -J m6.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m6.out -e m6.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R methyl mean neut mono"
wsbatch -J m7.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o m7.out -e m7.err --wrap="$cmd"




cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean tcell neut"
wsbatch -J hh2.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh2.out -e hh2.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean tcell mono"
wsbatch -J hh3.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh3.out -e hh3.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean mono neut"
wsbatch -J hh4.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh4.out -e hh4.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean mono tcell"
wsbatch -J hh5.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh5.out -e hh5.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean neut tcell"
wsbatch -J hh6.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh6.out -e hh6.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name_2.R hist mean neut mono"
wsbatch -J hh7.job --partition=shared-bigmem  --time=12:00:00 --mem-per-cpu=500000 -o hh7.out -e hh7.err --wrap="$cmd"



