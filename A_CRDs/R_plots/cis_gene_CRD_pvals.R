# plot cis gene-CRD p values
# Clean environment
rm(list=ls())
gc()


directory_s='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/significants/'
directory='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/nominal_merged/'
out_directory='Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots'

files=list.files(path = directory, pattern = ".txt")
files=intersect(list.files(path = directory_s, pattern = "significant.txt"),
                list.files(path = directory_s, pattern = "0.05"))

file=files[1]

for (file in files){
  crd_file=as.data.frame(data.table::fread(paste0(directory,file), head=FALSE, stringsAsFactors=FALSE))
  
  colnames(crd_file) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                  "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","dof1",
                                  "dof2","bml1","bml2","nom_pval","r_squared","slope","slope_se","adj_pvalue")
  
  colnames(crd_file) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                         "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end",
                         "nom_pval","r_squared","slope","besthit")
  head(crd_file)

  pdf(paste0(out_directory,"histogram_",file,".pdf"))
  
  hist(crd_file$nom_pval,xlab="adj_pvalue",main=paste0("title"), breaks=20)
  lines(density(crdsize.log), lty=2, lwd=2)
  abline(v=  median(crdsize.log) ,lty=2,lwd=3, col='red')
  # text(x=median(crdsize.log)+1,y=max(crdsize.log)*40,paste0("median:\n",round(median(crdsize)/1000,1)," kb"), col='red')
  
  dev.off()
}


#### median pval
hist_transCRD_sharing=c(0,45, 0.35, 0.46,0.58,0.27, 0.25)
methyl_transCRD_sharing=c(0.93,0.93,0.97,0.9,0.8,0.81)
median(hist_transCRD_sharing)
median(methyl_transCRD_sharing)