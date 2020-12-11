# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(UpSetR)
library(data.table)
library(corrplot)

#############################################################################################
#
# FUNCTION
#
#############################################################################################

# if 50% of the peak overlap, it is a shared CRD
# save the CRDs that are shared

compare_CRD <- function(query,reference,name,threshold=0.5){
  filename=paste0(name,"_sharedCRDs.txt")
  file.create(filename)

  CRD_IDs = unique(query$V2)
  n_replicated = 0
  for(j in 1:length(CRD_IDs)){
    current_CRD = CRD_IDs[j]
    current_peaks = query$V1[which(query$V2 == current_CRD)] # take peaks from query that are in query CRD
    idx_overlap = which(reference$V1 %in% current_peaks)
    if(length(idx_overlap)>0){
      overlapping_peaks_in_same_CRD = sort(table(reference$V2[idx_overlap]),decreasing=T)[1]
    } else {
      overlapping_peaks_in_same_CRD = 0
    }
    if(overlapping_peaks_in_same_CRD/length(current_peaks)>threshold){
      n_replicated = n_replicated + 1
      associated_CRDs=unique(reference$V2[reference$V1 %in% current_peaks])
      for (associated in associated_CRDs){
        line=paste(current_CRD ,associated,sep=" ")
        #cat(current_CRD, associated,sep="\t")
        write(line,file=paste0(name,"_sharedCRDs.txt"),append=TRUE)
      }
    }
  }
  list(total= n_replicated,fraction=n_replicated/length(CRD_IDs))
}


#############################################################################################
#
# LOAD DATA GUILLAUME
#
#############################################################################################

directory='/Users/dianaavalos/Programming/CRDs_nov20/debug/analysis_files/'
# ideal: update with FDR=0.05 and 0.01
# histone aCRDs
peakset_neut = as.data.frame(fread(paste0(directory,'EGAD00001002670.ALLchr.peaksID.txt'),header=F))
peakset_mono = as.data.frame(fread(paste0(directory,'EGAD00001002672.ALLchr.peaksID.txt'),header=F))
peakset_tcel = as.data.frame(fread(paste0(directory,'EGAD00001002673.ALLchr.peaksID.txt'),header=F))

mono_vs_neut = compare_CRD(peakset_mono,peakset_neut,paste0(directory,'hist_aCRD_mono_vs_neut'))
mono_vs_tcel = compare_CRD(peakset_mono,peakset_tcel,paste0(directory,'hist_aCRD_mono_vs_tcel'))

neut_vs_mono = compare_CRD(peakset_neut,peakset_mono,paste0(directory,'hist_aCRD_neut_vs_mono'))
neut_vs_tcel = compare_CRD(peakset_neut,peakset_tcel,paste0(directory,'hist_aCRD_neut_vs_tcel'))

tcel_vs_mono = compare_CRD(peakset_tcel,peakset_mono,paste0(directory,'hist_aCRD_tcel_vs_mono'))
tcel_vs_neut = compare_CRD(peakset_tcel,peakset_neut,paste0(directory,'hist_aCRD_tcel_vs_neut'))

#############################################################################################
#
# LOAD DATA DIANA
#
#############################################################################################
directory='/Users/dianaavalos/Programming/A_CRD_plots/quantify_ALL/'
# need to get peaksID
# histone aCRDs
peakset_neut = as.data.frame(fread(paste0(directory,'EGAD00001002670.ALLchr.peaksID.txt'),header=F))
peakset_mono = as.data.frame(fread(paste0(directory,'EGAD00001002672.ALLchr.peaksID.txt'),header=F))
peakset_tcel = as.data.frame(fread(paste0(directory,'EGAD00001002673.ALLchr.peaksID.txt'),header=F))


# debug
query=peakset_mono
reference=peakset_neut
name=paste0(directory,'mono_vs_neut')
threshold=0.5

### input data
# add sCRD and methyl CRDs

#############################################################################################
#
# PLOT
#
#############################################################################################


pdf("CRD_pairwise_comparisons_between_cell_types_Method_3.pdf")
M = matrix(c(1,neut_vs_mono$fraction,neut_vs_tcel$fraction,
             mono_vs_neut$fraction,1,mono_vs_tcel$fraction,
             tcel_vs_neut$fraction,tcel_vs_mono$fraction,1),
           ncol=3,byrow=T)
colnames(M) = c("Neutrophils","Monocytes","T cells")
rownames(M) = c("Neutrophils","Monocytes","T cells")
corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()
