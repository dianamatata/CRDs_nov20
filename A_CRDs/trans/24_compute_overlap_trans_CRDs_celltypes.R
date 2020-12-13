# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(tidyverse)
library(data.frame)


#############################################################################################
#
# FUNCTION
#
#############################################################################################

compute_shared_transCRD <- function(shared_crds,trans_crd_cell1,trans_crd_cell2,filename){
  # function written for cell1_vs_cell2_sharedCRDs.txt, applicable for all pairs afterwards
  
  
  # select CRD IDs in significant trans associations
  trans_crd_cell1_CRD_pairs=trans_crd_cell1 %>% select(4,8)
  trans_crd_cell2_CRD_pairs=trans_crd_cell2 %>% select(4,8)
  
  # select PAIRS in trans_crd_cell1_CRD_pairs that are also in shared_crds
  cell1shared=shared_crds$V1
  cell2shared=shared_crds$V2
  
  trans_crd_cell1_CRD_pairs_shared=trans_crd_cell1_CRD_pairs[(trans_crd_cell1_CRD_pairs$V4 %in% cell1shared) & (trans_crd_cell1_CRD_pairs$V8 %in% cell1shared), ]
  trans_crd_cell2_CRD_pairs_shared=trans_crd_cell2_CRD_pairs[(trans_crd_cell2_CRD_pairs$V4 %in% cell2shared) & (trans_crd_cell2_CRD_pairs$V8 %in% cell2shared), ]
  
  file.create(filename)
  # replace in trans_crd_cell1_CRD_pairs_shared with the equivalent in cell2 -B
  for (i in 1:length(trans_crd_cell1_CRD_pairs_shared[,1])) {
    crd_pair_cell1=c(trans_crd_cell1_CRD_pairs_shared[i,]$V4,trans_crd_cell1_CRD_pairs_shared[i,]$V8)
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V4)]
    cell2_equivalent_crd2 = shared_crds$V2[which(shared_crds$V1 == trans_crd_cell1_CRD_pairs_shared[i,]$V8)]
    crd_pair_cell1_in_cell2=c(cell2_equivalent_crd1 ,cell2_equivalent_crd2) # equivalent of cell1_pair in cell2_pair
    # pair_test=c("10_internal_7905",  "19_internal_6283") # idx  5160 : 10_internal_7905  19_internal_6283, from trans_crd_cell2_CRD_pairs_shared
    # find cell2_pair in trans_crd_cell2_CRD_pairs_shared
    indexes1=which(trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[2])
    indexes2=which(trans_crd_cell2_CRD_pairs_shared$V4 == crd_pair_cell1_in_cell2[1] &  trans_crd_cell2_CRD_pairs_shared$V8 == crd_pair_cell1_in_cell2[2])
    for (index in c(indexes1,indexes2)){
      if (length(trans_crd_cell2_CRD_pairs_shared[index,1]) != 0){
        line=paste0(trans_crd_cell2_CRD_pairs_shared[index,]$V4,' ',trans_crd_cell2_CRD_pairs_shared[index,]$V8)
        write(line,file=filename,append=TRUE)
      }
    }
  }
  result=as.data.frame(data.table::fread(filename), head=FALSE, stringsAsFactors=FALSE)
  result 
}


#################################### Folders and Files


path_transCRDs_signif = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant"
path_shared_crds="/Users/dianaavalos/Programming/CRDs_nov20/debug/analysis_files"

### loop on pairs of cells
FDRthreshold=0.01
trans_crd_neut = as.data.frame(data.table::fread(paste0(path_transCRDs_signif,'/hist_neut_mean_trans.significant_',FDRthreshold,'.txt'), head=FALSE, stringsAsFactors=FALSE))
trans_crd_mono = as.data.frame(data.table::fread(paste0(path_transCRDs_signif,'/hist_mono_mean_trans.significant_',FDRthreshold,'.txt'), head=FALSE, stringsAsFactors=FALSE))
trans_crd_tcell = as.data.frame(data.table::fread(paste0(path_transCRDs_signif,'/hist_tcell_mean_trans.significant_',FDRthreshold,'.txt'), head=FALSE, stringsAsFactors=FALSE))


#################################### Main

trans_crd_cell1=trans_crd_mono
trans_crd_cell2=trans_crd_neut
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_cell1,'_vs_',trans_crd_cell2,'.txt'), head=FALSE, stringsAsFactors=FALSE))
# shared crd computed from 23.R


# mono neut 
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_mono,'_vs_',trans_crd_neut,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_mono_neut=compute_shared_transCRD(shared_crds,trans_crd_mono,trans_crd_neut,paste0(path_shared_crds,'/mono_vs_neut_transCRDs_shared.txt'))

# neut mono
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_neut,'_vs_',trans_crd_mono,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_neut_mono=compute_shared_transCRD(shared_crds,trans_crd_neut,trans_crd_mono,paste0(path_shared_crds,'/neut_vs_mono_transCRDs_shared.txt'))

# neut tcell
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_neut,'_vs_',trans_crd_tcell,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_neut_tcell=compute_shared_transCRD(shared_crds,trans_crd_neut,trans_crd_tcell,paste0(path_shared_crds,'/neut_vs_tcell_transCRDs_shared.txt'))

# tcell neut
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_tcell,'_vs_',trans_crd_neut,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_tcell_neut=compute_shared_transCRD(shared_crds,trans_crd_tcell,trans_crd_neut,paste0(path_shared_crds,'/tcell_vs_neut_transCRDs_shared.txt'))

# mono tcell
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_mono,'_vs_',trans_crd_tcell,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_mono_tcell=compute_shared_transCRD(shared_crds,trans_crd_mono,trans_crd_tcell,paste0(path_shared_crds,'/mono_vs_tcell_transCRDs_shared.txt'))

# tcell mono
shared_crds = as.data.frame(data.table::fread(paste0(path_shared_crds,'/hist_shared_aCRDs_',trans_crd_tcell,'_vs_',trans_crd_mono,'.txt'), head=FALSE, stringsAsFactors=FALSE))
sharedtransCRDs_tcell_mono=compute_shared_transCRD(shared_crds,trans_crd_tcell,trans_crd_mono,paste0(path_shared_crds,'/tcell_vs_mono_transCRDs_shared.txt'))




# hist_aCRD_tcel_vs_neut_sharedCRDs.txt
# debug
trans_crd_cell1=trans_crd_neut
trans_crd_cell2=trans_crd_mono
filename=paste0(path_shared_crds,'/mono_vs_neut_transCRDs_shared.txt')
# File '/Users/dianaavalos/Programming/CRDs_nov20/debug/analysis_files/mono_CRDs_shared_with_neut.txt' has size 0. Returning a NULL data.table.


