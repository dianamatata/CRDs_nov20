#!/bin/bash

cmd="source /home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/CRD-QTLs/1_run_permutation_pass.sh"
wsbatch -J qtl.job --partition=bigmem-EL7 --time=36:50:00 -o 1.out -e 1.err --wrap="$cmd"
