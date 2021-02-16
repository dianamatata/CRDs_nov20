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
# FUNCTION: COMPARE CRDs
#
#############################################################################################

# if 50% of the peak overlap, it is a shared CRD
# save the CRDs that are shared

# query=peakset_neut
# reference=peakset_mono
# name=paste0(out_directory,paste0(name,'_neut_vs_mono'))
compare_CRD <- function(query,reference,name,threshold=0.5){
  filename=paste0(name,"_sharedCRDs1.txt")
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
        write(line,file=filename,append=TRUE)
      }
    }
  }
  list(total= n_replicated,fraction=n_replicated/length(CRD_IDs))
}


#############################################################################################
#
# FUNCTION: COMPUTE CORRELATIONS
#
############################################################################################# 

compute_correlation_matrix_CRD_sharing <- function(peakset_neut,peakset_mono,peakset_tcell, name, out_directory){
  
  neut_vs_mono = compare_CRD(peakset_neut,peakset_mono,paste0(out_directory,paste0(name,'_neut_vs_mono')))
  neut_vs_tcell = compare_CRD(peakset_neut,peakset_tcell,paste0(out_directory,paste0(name,'_neut_vs_tcell')))
  mono_vs_neut = compare_CRD(peakset_mono,peakset_neut,paste0(out_directory,paste0(name,'_mono_vs_neut')))
  mono_vs_tcell = compare_CRD(peakset_mono,peakset_tcell,paste0(out_directory,paste0(name,'_mono_vs_tcell')))
  tcell_vs_mono = compare_CRD(peakset_tcell,peakset_mono,paste0(out_directory,paste0(name,'_tcell_vs_mono')))
  tcell_vs_neut = compare_CRD(peakset_tcell,peakset_neut,paste0(out_directory,paste0(name,'_tcell_vs_neut')))
  
  overlap_array_input=c(neut_vs_mono$fraction, neut_vs_tcell$fraction, mono_vs_neut$fraction, mono_vs_tcell$fraction , tcell_vs_mono$fraction, tcell_vs_neut$fraction)
}


#############################################################################################
#
# FUNCTION: PLOT
#
#############################################################################################


plot_correlation_matrix_CRD_sharing <- function(overlap_array_input, name, plot_directory){
  
  neut_vs_mono=overlap_array_input[1]
  neut_vs_tcell=overlap_array_input[2]
  mono_vs_neut=overlap_array_input[3]
  mono_vs_tcell=overlap_array_input[4]
  tcell_vs_mono=overlap_array_input[5]
  tcell_vs_neut=overlap_array_input[6]
  
  pdf(paste0(plot_directory,name,"_pairwise_comparisons.pdf"))
  M = matrix(c(1,neut_vs_mono,neut_vs_tcell,
               mono_vs_neut,1,mono_vs_tcell,
               tcell_vs_neut,tcell_vs_mono,1),
             ncol=3,byrow=T)
  colnames(M) = c("Neutrophils","Monocytes","T cells")
  rownames(M) = c("Neutrophils","Monocytes","T cells")
  corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
  corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
}


#############################################################################################
#
# MAIN
#
#############################################################################################

setwd('/Users/dianaavalos')
directory='Programming/A_CRD_plots/trans_files/7_CRD_peaks/peaks/'
out_directory='/Users/dianaavalos/Programming/A_CRD_plots/CRD_sharing/'
plot_directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_peaks/plots/'

condition='mean'
data_type='hist'
for(data_type in c('hist','methyl')){
  name=paste0(data_type,"_", condition)
  
  #### MAIN
  peakset_neut = as.data.frame(fread(paste0(directory,data_type,"_neut_",condition,".ALLchr.peaks.txt"),header=F))
  peakset_mono = as.data.frame(fread(paste0(directory,data_type,"_mono_",condition,".ALLchr.peaks.txt"),header=F))
  peakset_tcell = as.data.frame(fread(paste0(directory,data_type,"_tcell_",condition,".ALLchr.peaks.txt"),header=F))
  
  overlap_array_input=compute_correlation_matrix_CRD_sharing(peakset_neut,peakset_mono,peakset_tcell, name, out_directory)
  # plot_correlation_matrix_CRD_sharing(overlap_array_input, name, plot_directory)
  
}

#############################################################################################
#
# check up these are actually shared
#
#############################################################################################


dir_CRD_QTL_signif='/Users/dianaavalos/Programming/A_CRD_plots/10_CRD_QTL/conditional_merged/'

file="/Users/dianaavalos/Programming/A_CRD_plots/CRD_sharing/hist_mean_neut_vs_mono_sharedCRDs.txt"
shared_file=as.data.frame(fread(file),head=FALSE)
data_type='hist'
cell1='neut'
cell2='mono'
module='mean'
crd_qtl_cell1_signif=as.data.frame(data.table::fread(paste0(dir_CRD_QTL_signif,data_type,'_',cell1,'_',module,'_ALL.txt.gz'), head=FALSE, stringsAsFactors=FALSE))
crd_qtl_cell2_signif=as.data.frame(data.table::fread(paste0(dir_CRD_QTL_signif,data_type,'_',cell2,'_',module,'_ALL.txt.gz'), head=FALSE, stringsAsFactors=FALSE))

colnames(crd_qtl_cell1_signif) =colnames(crd_qtl_cell2_signif) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                   "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","rank",
                                   "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig","bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")

head(crd_qtl_cell1_signif)


for (i in seq(100))
{

  crd1_cell1=shared_file[i,1]
  crd1_cell2=shared_file[i,2]
  # cat(crd1_cell1, ' ',crd1_cell2)
  
  sub1=crd_qtl_cell1_signif %>% filter(phenotype_ID %in%  crd1_cell1)   
  # cat(sub1[1,]$phenotype_ID, ' ', sub1[1,]$phenotype_ID_start, '  ', sub1[1,]$phenotype_ID_end)
  sub2=crd_qtl_cell2_signif %>% filter(phenotype_ID %in%  crd1_cell2)   
  # cat(sub2[1,]$phenotype_ID, ' ', sub2[1,]$phenotype_ID_start, '  ', sub2[1,]$phenotype_ID_end)
  
  dist=sub2[1,]$phenotype_ID_start-sub1[1,]$phenotype_ID_start
  cat(dist, '  ')
}

