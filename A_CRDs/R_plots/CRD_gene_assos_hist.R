# pval crd gene assoc
rm(list=ls())
gc()


directory_s='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/significants/'
directory='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/nominal_merged/'
out_directory='Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/plots'

files=intersect(list.files(path = directory, pattern = "mapping_CRD_gene_ALL.txt.gz"),
                list.files(path = directory, pattern = "mean"))

file=files[1]
for (file in files){
  crd_file=as.data.frame(data.table::fread(paste0(directory,file), head=FALSE, stringsAsFactors=FALSE))
  
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

colnames(crd_file) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                       "nb_variants","distance","var_ID","var_ID_chr","var_ID_start","var_ID_end","dof1",
                       "dof2","bml1","bml2","nom_pval","r_squared","slope","slope_se","adj_pvalue")
