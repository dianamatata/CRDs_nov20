# Clean environment ---------------------------------

rm(list=ls())
gc()

# Packages ---------------------------------

library(ggplot2)
library(corrplot)
library(reshape2)
library(data.table)
library(ggmosaic)

# Function ---------------------------------

plot_correlation_matrix_CRD_sharing <- function(Matrix, name, plot_directory){
  
  pdf(paste0(plot_directory,name,"_sharing.pdf"))
  colnames(Matrix) = c("Neutrophils","Monocytes","T cells")
  rownames(Matrix) = c("Neutrophils","Monocytes","T cells")
  corrplot(Matrix,is.corr=F,cl.lim = c(0, 100),p.mat = Matrix,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
}

# Plots ---------------------------------
plot_directory='/Users/dianaavalos/Desktop/'
# pi1 estimate

histmean = matrix(c(1,85.3 ,86.3,
                  80.5,1,86,
                  92.8,95.3,1),ncol=3,byrow=T)

methylmean = matrix(c(1,83.6 ,85.7,
                    84.3,1,86.8 ,
                    89.6,88.6,1),ncol=3,byrow=T)

histloom = '0'

methylloom ='0'

name='histmean'
plot_correlation_matrix_CRD_sharing(histmean, name, plot_directory)


name='methylmean'
plot_correlation_matrix_CRD_sharing(methylmean, name, plot_directory)

