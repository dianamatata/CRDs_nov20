# Clean environment ---------------------------------

rm(list=ls())
gc()


# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(corrplot)


# Cluster Directories ---------------------------------

directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/merged_1000/'
plot_directory='/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/plots_pval/'

# Mac Directories ---------------------------------

directory='/Users/dianaavalos/Programming/A_CRD_plots/10_CRD_QTL/mixed/'
plot_directory='/Users/dianaavalos/Programming/A_CRD_plots/10_CRD_QTL/mixed/'

# Main Loop ---------------------------------

file_list=list.files(path=directory,pattern=c('permuts.txt.gz'))

for (file in file_list) {
  name=substr(file, 1, nchar(file)-15)
  df = as.data.frame(fread(paste0(directory,file),header=F))
  # Hist P val ---------------------------------
  
  qobj=qvalue(df$V20)
  pi0=qobj$pi0
  pi1=1-pi0
  pdf(paste0(plot_directory,"histogram_p_",name,".pdf"))
  hist(df$V20, main=paste0(' pi0=',round(pi0,digits = 3),'  pi1=',round(pi1,digits = 3)), xlab='p values', breaks=20,cex.main=1.5,  cex.lab=1.2,cex.axis=1.2)
  dev.off()
}

# Corr plot ---------------------------------

overlap_array_input=c(neut_vs_mono$fraction, neut_vs_tcell$fraction, mono_vs_neut$fraction, mono_vs_tcell$fraction , tcell_vs_mono$fraction, tcell_vs_neut$fraction)


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

plot_correlation_matrix_CRD_sharing(overlap_array_input, name, plot_directory)

