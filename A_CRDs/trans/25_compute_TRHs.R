
#############################################################################################
#
# NOT ORDERED THINGS
#
#############################################################################################


# 
df = data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
df <-df[order(-df$SIZE),]
df['GENE']=0

# plot gene number and plot the nbr of CRDs and genes associated per TRH

path2='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes'
for (i in 1:N_TRH) {
  file=file.path(path2,sprintf("genesonly_70_cluster%s.txt", i) ) ######## cell type
  if (file.exists(file)){
    df$GENE[i]=length(read.table(f, header = FALSE, sep = "", dec = ".")$V1)
    #sprintf("70_cluster%s.txt", i)
    #print(length(read.table(file, header = FALSE, sep = "", dec = ".")$V1))
  }
}

