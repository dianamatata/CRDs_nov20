# in Hi-C analysis
library(qvalue)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(tidyverse)

# path:/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/mapping_aCRD_gene%
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CORRELATION/expression

cell_types = list('70','72','73')
cell_type='70'

if (cell_type=='70'){
  # for neutrophils
  allCRDs = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz',header=F)
  mapdata = read.table('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/mapping_gene_CRD_mean_ALL_70.txt',stringsAsFactors=F)
  d = read.table("/Users/dianaavalos/Programming/Hi-C_correlated_peaks/gene_CRD_mean_permutations_full_70.txt.gz", hea=F, stringsAsFactors=F)
  cordata = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002675_RNA.ALL.txt.gz')
  data='neutrophils'
} else if( cell_type=='72'){
  # for monocytes
  allCRDs = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002672_ALL.modules.MOD1.NRE2.txt.gz',header=F)
  mapdata = read.table('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/mapping_gene_CRD_mean_ALL_72.txt',stringsAsFactors=F)
  d = read.table("/Users/dianaavalos/Programming/Hi-C_correlated_peaks/gene_CRD_mean_permutations_full_72.txt.gz", hea=F, stringsAsFactors=F)
  cordata = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002674_RNA.ALL.txt.gz')
  data='monocytes'
} else if( cell_type=='73'){
  # for t-cells
  allCRDs = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002673_ALL.modules.MOD1.NRE2.txt.gz',header=F)
  mapdata = read.table('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/mapping_gene_CRD_mean_ALL_73.txt',stringsAsFactors=F)
  d = read.table("/Users/dianaavalos/Programming/Hi-C_correlated_peaks/gene_CRD_mean_permutations_full_73.txt.gz", hea=F, stringsAsFactors=F)
  cordata = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/EGAD00001002671_RNA.ALL.txt.gz')
  data='tcells'
} 


protein_coding_genes = scan("/Users/dianaavalos/Programming/reference_files/gencode.v15.annotation.protein_coding.gene_id.txt",what="")
long_nc_RNA_genes = scan("/Users/dianaavalos/Programming/reference_files/gencode.v15.annotation.long_noncoding_RNAs.gene_id.txt",what="")
PCHiC = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/PCHiC_peak_matrix_cutoff5.tsv')
genelist = c(protein_coding_genes,long_nc_RNA_genes)

colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
 "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")
#Filter protein coding and long nc RNA genes
mapdata = mapdata[mapdata$phenotype_ID %in% genelist,]
mapdata = mapdata[order(mapdata[,2],mapdata[,3]),]
nb_CRD_not_associated = nrow(allCRDs) - length(unique(mapdata$CRD_ID))

colnames(d) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
"nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","degree_freedom",
"Dummy","1st_param_beta","2nd_param_beta","nominal_pval","slope","empirical_pval","beta_pval")
d = d[d$phenotype_ID %in% genelist,]

cordata$distance = cordata$V4-cordata$V3
cordata = cordata[,c(3,4,5,6,8,9)]
colnames(cordata)[1:5] = c("pos1","pos2","gene1","gene2","pval")
cordata$pval.adj = p.adjust(cordata$pval,method='fdr')
cordata$sameCRD = -1
setkeyv(cordata,c("gene1","gene2"))

cordata = cordata[cordata$gene1 %in% genelist,]
cordata = cordata[cordata$gene2 %in% genelist,]

cordata$rowid = c(1:nrow(cordata))
nb_genes_not_associated = length(unique(cordata$gene1)) - length(unique(mapdata$phenotype_ID))

colnames(PCHiC)[1] = "baitChr"

baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))

genebed <- GRanges(seqnames=mapdata$phenotype_ID_chr,ranges=IRanges(start=mapdata$phenotype_ID_start, end=mapdata$phenotype_ID_end))
CRDbed <- GRanges(seqnames=mapdata$CRD_ID_chr,ranges=IRanges(start=mapdata$CRD_ID_start, end=mapdata$CRD_ID_end))

#fwd
x = findOverlaps(baitbed,genebed)
y = findOverlaps(oebed,CRDbed)
tmp = rbind(as.data.frame(x),as.data.frame(y))
validated.fwd = tmp[which(duplicated(tmp)),]

#bwd
x = findOverlaps(baitbed,CRDbed)
y = findOverlaps(oebed,genebed)
tmp = rbind(as.data.frame(x),as.data.frame(y)) 
validated.bwd = tmp[which(duplicated(tmp)),]

validated = unique(rbind(validated.fwd,validated.bwd))
# 
# hic_validated = rep(1,nrow(PCHiC))
# for(i in 1:nrow(validated)){
#     currenthic = validated$queryHits[i]
#     currentpval = mapdata$nominal_pval[validated$subjectHits[i]]
#     if(currentpval<hic_validated[currenthic]){ # <1
#         hic_validated[currenthic] = currentpval
#     }
# }

mapdata_validated = rep(0,nrow(mapdata))
for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentmap = validated$subjectHits[i]
    if (cell_type=='70'){
      currentHiCScore = PCHiC$Neu[currenthic] 
    }else if (cell_type=='72'){
      currentHiCScore = PCHiC$Mon[currenthic]
    }else if (cell_type=='73'){
      currentHiCScore = PCHiC$nCD4[currenthic] ##### depend on cell type so maybe redo plots with HiC?
    }
    if(currentHiCScore>mapdata_validated[currentmap]){
        mapdata_validated[currentmap] = currentHiCScore
    }
}

mapdata$HIC=mapdata_validated

crd_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)

compute_ratio_hic_mapdata <-function(CRDmindist,CRDmaxdist,cutoff=5){
  up = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist & mapdata_validated>cutoff,na.rm=T)
  down = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist,na.rm=T)
  up/down*100
}

mat_hic = matrix(0,nrow=(length(crd_dist_bins)-1),ncol=1)

for(i in 1:(length(crd_dist_bins)-1)){
  mat_hic[i,1] = compute_ratio_hic_mapdata(crd_dist_bins[i],crd_dist_bins[i+1])
}



# plots
pdf("Distance_analysis.pdf") # Relative positive of genes within CRDs histogram
hist(mapdata$distance[which(mapdata$distance!=0)],n=101,main='Distance between genes and CRDs',xlab='Distance [bp]',cex.lab=1.3,cex.axis=1.3)
relative_position = (mapdata$phenotype_ID_start[which(mapdata$distance==0)]-mapdata$CRD_ID_start[which(mapdata$distance==0)])/(mapdata$CRD_ID_end[which(mapdata$distance==0)]-mapdata$CRD_ID_start[which(mapdata$distance==0)])
hist(relative_position,xlim=c(0,1),n=50, main='Relative positive of genes within CRDs',xlab='Relative position',cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Histogram_effect_size.pdf")
hist(mapdata$slope,n=101,main='Effect size of the associations between genes and CRDs',xlab='Slope',cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Histogram_pvalue.pdf")
hist(d$beta_pval,main=paste0("Histogram of adjusted p-values%\n",nrow(mapdata)," significant associations at 5% FDR"),xlab="P-values",cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Boxplot_pvalue_by_distance.pdf")
within_CRD = -log10(mapdata$nominal_pval[which(mapdata$distance==0)])
less_1kb = -log10(mapdata$nominal_pval[which(mapdata$distance>0 & mapdata$distance<1e03)])
less_10kb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e03 & mapdata$distance<1e04)])
less_100kb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e04 & mapdata$distance<1e05)])
less_1Mb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e05 & mapdata$distance<1e06)])
boxplot(within_CRD,less_1kb,less_10kb,less_100kb,less_1Mb,names=c("within_CRD","<1kb","<10kb","<100kb","<1Mb"),ylab="-log10(p-value)",cex.lab=1.3,cex.axis=1.3)
dev.off()
#stop()

pdf("Connectivity_CRD_gene.pdf")
CRDhist = hist(table(mapdata$CRD_ID),breaks=c(0,1,2,3,4,100),plot=F)
genehist = hist(table(mapdata$phenotype_ID),breaks=c(0,1,2,3,100),plot=F)
barplot(c(nb_CRD_not_associated,CRDhist$counts),names=c("0","1","2","3","4","5+"),main="Number of genes associated with each CRD",cex.names=1.5,cex.axis=1.5)
barplot(c(nb_genes_not_associated,genehist$counts),names=c("0","1","2","3","4+"),main="Number of CRDs associated with each gene",cex.names=1.5,cex.axis=1.5)

ggplot(data.frame(counts = c(nb_CRD_not_associated,CRDhist$counts),Number = c("0","1","2","3","4","5+")), aes(x = Number, y = counts))+ ggtitle("Genes associated with each CRD") +
  geom_bar(stat = "identity",fill="#E69F00") +
  geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated genes",y = "CRD counts") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

 ggplot(data.frame(counts = c(nb_genes_not_associated,genehist$counts),Number = c("0","1","2","3","4+")), aes(x = Number, y = counts))+ ggtitle("CRDs associated with each gene") +
    geom_bar(stat = "identity",fill="#56B4E9") +
    geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated CRDs",y = "Gene counts") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

dev.off()

pdf("HiC_support_gene_CRD_associations.pdf")
toplot = data.frame(counts = mat_hic[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene-CRD associations") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

#######################
#Co-expression analysis
#######################

#Annotate cordata with CRD mapping

CRD_IDs = unique(mapdata$CRD_ID)
IDXs = c()
CRDdistance = c()
CRDdistance_gene1=c()
CRDdistance_gene2=c()
hicsupport = c()
cordata$dist1=-1
cordata$dist2=-1

for(i in 1:length(CRD_IDs)){
  # cat(i,"\n")
  idx = which(mapdata$CRD_ID == CRD_IDs[i])
  if(length(idx)>1){
    for(k in 1:(length(idx)-1)){
    	for(l in (k+1):length(idx)){
      	geneA = mapdata$phenotype_ID[idx[k]]
      	geneB = mapdata$phenotype_ID[idx[l]]
        myhit = c(unlist(cordata[gene1 == geneB & gene2 == geneA,9]),unlist(cordata[gene1 == geneA & gene2 == geneB,9]))
        if(length(myhit)>0){
                #cat(CRD_IDs[i],geneA,geneB,sep='\n')
          IDXs = c(IDXs,myhit)
          tmp = mean(abs(mapdata$distance[idx[k]]),abs(mapdata$distance[idx[l]])) # distance of the gene to the CRD CRD_IDs[i]
          tmp2 = mean(mapdata_validated[idx[k]],mapdata_validated[idx[l]])
          hicsupport = c(hicsupport,tmp2)
          
          CRDdistance = c(CRDdistance,tmp) # mean distance with the 2 genes
          CRDdistance_gene1 = c(CRDdistance_gene1,abs(mapdata$distance[idx[k]])) # dist of gene 1
          CRDdistance_gene2 = c(CRDdistance_gene2,abs(mapdata$distance[idx[l]])) # dist of gene 2
        }
      }
    }
  }
}

tmpdf = data.table(idx=IDXs,dist=CRDdistance,hicsupport=hicsupport, dist1=CRDdistance_gene1, dist2=CRDdistance_gene2)
tmpdf <- tmpdf %>%
    group_by(idx) %>%
    summarize(mindist=min(dist),hic=mean(hicsupport), mindist1=min(dist1), mindist2=min(dist2))

cordata$sameCRD[tmpdf$idx] = tmpdf$mindist # by default cordata$sameCRD=-1. replace it with the min distance to the CRDs
cordata$hic = 0
cordata$hic[tmpdf$idx] = tmpdf$hic
cordata$dist1[tmpdf$idx] = tmpdf$mindist1
cordata$dist2[tmpdf$idx] = tmpdf$mindist2

compute_ratio <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,pval.cutoff=0.01){
    if(coexpressed){
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    } else {
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    }
    up/down*100
}

# new 
compute_ratio_n <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,pval.cutoff=0.01){
  if(coexpressed){
      # cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist
      up_both = sum( cordata$sameCRD==0 & (cordata$dist1 + cordata$dist2 ==0) & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
      up_one = sum( cordata$sameCRD==0 & (cordata$dist1 + cordata$dist2 !=0) & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
      down = sum(cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    } else {
    up_both = sum( cordata$sameCRD==0 & (cordata$dist1 + cordata$dist2 ==0) & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    up_one = sum( cordata$sameCRD==0 & (cordata$dist1 + cordata$dist2 !=0) & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    down = sum(cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
  }
  c(up_both/down*100, up_one/down*100)
}

compute_OR <-function(mindist,maxdist,pval.cutoff=0.01){
    coexpressed_same_CRD = sum(cordata$sameCRD>=0 & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    coexpressed_notsame_CRD = sum(!(cordata$sameCRD>=0) & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    notcoexpressed_same_CRD = sum(cordata$sameCRD>=0 & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    notcoexpressed_notsame_CRD = sum(!(cordata$sameCRD>=0) & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    m = matrix(c(coexpressed_same_CRD,coexpressed_notsame_CRD,notcoexpressed_same_CRD,notcoexpressed_notsame_CRD),ncol=2)
    pval = fisher.test(m)$p.value
    OR = fisher.test(m)$estimate
    list(OR=OR,pval=pval)
}

OR_less_10kb = compute_OR(0,1e04)
OR_less_20kb = compute_OR(1e04,2e04)
OR_less_50kb = compute_OR(2e04,5e04)
OR_less_100kb = compute_OR(5e04,1e05)
OR_less_200kb = compute_OR(1e05,2e05)
OR_less_500kb = compute_OR(2e05,5e05)
OR_less_1Mb = compute_OR(5e05,1e06)
OR_all = c(OR_less_10kb[[1]],OR_less_20kb[[1]],OR_less_50kb[[1]],OR_less_100kb[[1]],OR_less_200kb[[1]],OR_less_500kb[[1]],OR_less_1Mb[[1]])
pval_all = c(OR_less_10kb[[2]],OR_less_20kb[[2]],OR_less_50kb[[2]],OR_less_100kb[[2]],OR_less_200kb[[2]],OR_less_500kb[[2]],OR_less_1Mb[[2]])


gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06) # dist btw genes
crd_dist_bins = c(0,1,1e04,1e05,1e06) # dist to CRD
coexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))
notcoexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))
for(i in 1:(length(gene_dist_bins)-1)){
    for(k in 1:(length(crd_dist_bins)-1)){
        coexpressed_mat[i,k] = compute_ratio(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
        notcoexpressed_mat[i,k] = compute_ratio(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
    }
}

# new version
gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06) # dist btw genes
crd_dist_bins = c(0,1,1e04,1e05,1e06) # dist to CRD, here to change
coexpressed_mat_n = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)))
notcoexpressed_mat_n = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)))
for(i in 1:(length(gene_dist_bins)-1)){
  for(k in 1:(length(crd_dist_bins)-1)){
    if (k>1){
      km=k+1
      coexpressed_mat_n[i,km] = compute_ratio(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])   
      notcoexpressed_mat_n[i,km] = compute_ratio(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
    } else {
      coexpressed_mat_n[i,c(k,k+1)]=compute_ratio_n(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
      notcoexpressed_mat_n[i,c(k,k+1)] = compute_ratio_n(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
    }
  }
}
coexpressed_mat_n2=coexpressed_mat_n
notcoexpressed_mat_n2=notcoexpressed_mat_n
coexpressed_mat2=coexpressed_mat
notcoexpressed_mat2=notcoexpressed_mat

coexpressed_mat = data.frame(coexpressed_mat)
colnames(coexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
coexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
coexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
coexpressed_table <- melt(coexpressed_mat, id.vars = "row")
colnames(coexpressed_table) = c("Gene","CRD","value")

notcoexpressed_mat = data.frame(notcoexpressed_mat)
colnames(notcoexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
notcoexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
notcoexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
notcoexpressed_table <- melt(notcoexpressed_mat, id.vars = "row")
colnames(notcoexpressed_table) = c("Gene","CRD","value")

coexpressed_mat_n = data.frame(coexpressed_mat_n)
colnames(coexpressed_mat_n) = c("inside2","inside1","<10kb","<100kb","<1Mb")
coexpressed_mat_n$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
coexpressed_mat_n$row <- factor(coexpressed_mat_n$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
coexpressed_table_n <- melt(coexpressed_mat_n, id.vars = "row")
colnames(coexpressed_table_n) = c("Gene","CRD","value")

notcoexpressed_mat_n = data.frame(notcoexpressed_mat_n)
colnames(notcoexpressed_mat_n) = c("inside2","inside1","<10kb","<100kb","<1Mb")
notcoexpressed_mat_n$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
notcoexpressed_mat_n$row <- factor(notcoexpressed_mat_n$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
notcoexpressed_table_n <- melt(notcoexpressed_mat_n, id.vars = "row")
colnames(notcoexpressed_table_n) = c("Gene","CRD","value")


# ### Fig 3.7 new
Paper=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")
color_palette= Paper
# size 572 666
pdf("Co-expression_analysis_CRD_gene.pdf",paper="a4r")

ggplot(coexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + geom_bar(stat = "identity") + 
  scale_fill_manual(values = color_palette)  + ylim(0, 65) + 
  geom_text(  aes(   y=rep (  apply(coexpressed_mat_n[,1:5],1,sum  ),5)  ,  label = c(paste0("(",round(OR_all),")"),rep("",28))  ),vjust = -.5,  size=5) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

ggplot(notcoexpressed_table_n, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes") + ylim(0, 10) +
  geom_bar(stat = "identity") +scale_fill_manual(values = color_palette) +
  labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

dev.off()

# ### Fig 3.7: 

pdf("Co-expression_analysis_CRD_gene.pdf",paper="a4r")


ggplot(coexpressed_table, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + geom_bar(stat = "identity") + scale_fill_manual(values = color_palette)  + ylim(0, 65) + 
  geom_text(aes(y=rep(apply(coexpressed_mat[,1:4],1,sum),4),label = c(paste0("(",round(OR_all),")"),rep("",21))),vjust = -.5,  size=5) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

ggplot(notcoexpressed_table, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes") + ylim(0, 10) +
  geom_bar(stat = "identity") +scale_fill_manual(values = color_palette) +
  labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))

dev.off()




#HiC data

compute_ratio_hic <-function(coexpressed=T,CRDmindist,CRDmaxdist,mindist,maxdist,pval.cutoff=0.01,threshold=5){
    if(coexpressed){
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$hic>threshold & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    } else {
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$hic>threshold & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    }
    up/down*100
}

#HiC support by gene to gene distance

gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06,3e09)
coexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
notcoexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)

for(i in 1:(length(gene_dist_bins)-1)){
        coexpressed_mat_hic[i,1] = compute_ratio_hic(T,0,1e09,gene_dist_bins[i],gene_dist_bins[i+1])
        coexpressed_mat_hic[i,2] = compute_ratio_hic(T,-1,0,gene_dist_bins[i],gene_dist_bins[i+1])
}

pdf("HiC_support_gene_pairs_with_same_CRD.pdf")
toplot = data.frame(counts = coexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene pairs associated with CRD") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


#HiC support by gene to CRD distance

CRD_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
coexpressed_mat_hic_crddist = matrix(0,nrow=(length(CRD_dist_bins)-1),ncol=1)

for(i in 1:(length(CRD_dist_bins)-1)){
        coexpressed_mat_hic_crddist[i,1] = compute_ratio_hic(T,CRD_dist_bins[i],CRD_dist_bins[i+1],0,1e09)
}

pdf("HiC_support_gene_pairs_with_same_CRD_geneCRD_dist.pdf")
toplot = data.frame(counts = coexpressed_mat_hic_crddist[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene pairs associated with CRD") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

save.image('Gene_CRD_associations.RDAta')


