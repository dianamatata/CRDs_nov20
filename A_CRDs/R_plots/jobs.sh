
public-bigmem 4 Days	250GB days-hours:minutes:seconds 4-00:00:00


Goals:

R script modulable avec les args and launch all the different params

wsbatch -J m1.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m1.out -e m1.err --wrap="$cmd"



cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean tcell neut"
wsbatch -J m2.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m2.out -e m2.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean tcell mono"
wsbatch -J m3.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m3.out -e m3.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean mono neut"
wsbatch -J m4.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m4.out -e m4.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean mono tcell"
wsbatch -J m5.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m5.out -e m5.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean neut tcell"
wsbatch -J m6.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m6.out -e m6.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R methyl mean neut mono"
wsbatch -J m7.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o m7.out -e m7.err --wrap="$cmd"




cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean tcell neut"
wsbatch -J mh2.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh2.out -e mh2.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean tcell mono"
wsbatch -J mh3.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh3.out -e mh3.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean mono neut"
wsbatch -J mh4.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh4.out -e mh4.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean mono tcell"
wsbatch -J mh5.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh5.out -e mh5.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean neut tcell"
wsbatch -J mh6.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh6.out -e mh6.err --wrap="$cmd"

cmd="Rscript /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/R_plots/CRD-QTLs_sharing_name.R hist mean neut mono"
wsbatch -J mh7.job --partition=public-bigmem  --time=4-00:00:00 --mem-per-cpu=250000 -o mh7.out -e mh7.err --wrap="$cmd"



