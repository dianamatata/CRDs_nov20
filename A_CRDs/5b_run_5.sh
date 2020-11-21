#!/bin/bash

cmd="source 5_run_mapping_CRD_gene.sh"
wsbatch -J 5_script.job --partition=mono-EL7 --time=10:10:00 -o 5_script.out -e 5_script.err --wrap="$cmd"

