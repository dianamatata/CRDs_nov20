library(UpSetR) # Creates visualizations of intersecting sets using a novel matrix design
library(data.table) # like pandas in python
library(corrplot) # graphical display of a correlation matrix, confidence interval

# neut 70, mono 72, tcell 73

#### Method 1: Look at the overlap in the surface covered by the CRDs, with the intersect method in bedtools
# this method counts the bp inside CRDs that overlap in CRDs in other tissue.
# We compute the bedtools intersect btw 2 cell types, we merge and order in case there are fragments overlapping, and 
# in 11_compute_bp_intersect.py we compute how many bp in total are overlapping
# (bp intersect mono and neut/bp neut= %shared in neut by mono)

# plot
# plot(rnorm(50), rnorm(50))

pdf("CRD_pairwise_comparisons_between_cell_types.pdf Method 1")
M0 = matrix(c(1,0.5559791012801089, 0.29371452801106906,
              0.517517715921954, 1, 0.2821886091864699,
              0.4647242130510857, 0.4796700438785532, 1),ncol=3,byrow=T)
colnames(M0) = c("Neutrophils","Monocytes","T cells")
rownames(M0) = c("Neutrophils","Monocytes","T cells")
corrplot(M0, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M0,is.corr=F,cl.lim = c(0, 1),p.mat = M0,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()

#### Method 2: Look at the shared peaks belonging to CRDs among cell types
# very small results, up to 3% of shared peaks
# these are the values for the bar plot

peakset_neut = as.data.frame(fread('../peaks/EGAD00001002670.ALLchr.peaksID',header=F))
peakset_mono = as.data.frame(fread('../peaks/EGAD00001002672.ALLchr.peaksID',header=F))
peakset_tcel = as.data.frame(fread('../peaks/EGAD00001002673.ALLchr.peaksID',header=F))

neut = length(peakset_neut$V1)
mono = length(peakset_mono$V1)
tcel = length(peakset_tcel$V1)
all = length(intersect(intersect(peakset_neut$V1,peakset_mono$V1),peakset_tcel$V1))
neut_and_mono = length(intersect(peakset_neut$V1,peakset_mono$V1))
neut_and_tcel = length(intersect(peakset_neut$V1,peakset_tcel$V1))
tcel_and_mono = length(intersect(peakset_tcel$V1,peakset_mono$V1))

expressionInput2 <- c(neutrophil = (neut-all-(neut_and_mono-all)-(neut_and_tcel-all)), monocyte = (mono-all-(tcel_and_mono-all)-(neut_and_mono-all)), tcell = (tcel-all-(neut_and_tcel-all)-(tcel_and_mono-all)), `neutrophil&monocyte` = neut_and_mono-all,
                     `neutrophil&tcell` = neut_and_tcel-all,`monocyte&tcell` = tcel_and_mono-all, `neutrophil&monocyte&tcell` = all)

# plot
pdf("CRD_pairwise_comparisons_between_cell_types.pdf")
M2 = matrix(c(1,neut_and_mono/neut,neut_and_tcel/neut,
              neut_and_mono/mono,1,tcel_and_mono/mono,
              neut_and_tcel/tcel,tcel_and_mono/tcel,1),ncol=3,byrow=T)
colnames(M2) = c("Neutrophils","Monocytes","T cells")
rownames(M2) = c("Neutrophils","Monocytes","T cells")
corrplot(M2, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M2,is.corr=F,cl.lim = c(0, 1),p.mat = M2,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()

pdf("Overlap_CRD_beds2.pdf",paper='a4r')
upset(fromExpression(expressionInput2), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput2), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")

dev.off()

#### Method 3: Look at the shared peaks belonging to CRDs among cell types
# take the peaks of 1 CRD, what is the % of shared in the other cell type? if more than 50% shared, the CRD is shared amon the 2 cell types

peakset_neut = as.data.frame(fread('../peaks/EGAD00001002670.ALLchr.peaksID',header=F))
peakset_mono = as.data.frame(fread('../peaks/EGAD00001002672.ALLchr.peaksID',header=F))
peakset_tcel = as.data.frame(fread('../peaks/EGAD00001002673.ALLchr.peaksID',header=F))

query=peakset_neut
reference=peakset_mono
threshold=0.5
j=1


compare_CRD <- function(query,reference,threshold=0.5){
  CRD_IDs = unique(query$V2)
  n_replicated = 0
  for(j in 1:length(CRD_IDs)){
    current_CRD = CRD_IDs[j]
    current_peaks = query$V1[which(query$V2 == current_CRD)]
    idx_overlap = which(reference$V1 %in% current_peaks)
    if(length(idx_overlap)>0){
      overlapping_peaks_in_same_CRD = sort(table(reference$V2[idx_overlap]),decreasing=T)[1]
    } else {
      overlapping_peaks_in_same_CRD = 0
    }
    if(overlapping_peaks_in_same_CRD/length(current_peaks)>threshold){
      n_replicated = n_replicated + 1
    }
  }
  list(total= n_replicated,fraction=n_replicated/length(CRD_IDs))
}

mono_vs_neut = compare_CRD(peakset_mono,peakset_neut)
mono_vs_tcel = compare_CRD(peakset_mono,peakset_tcel)

neut_vs_mono = compare_CRD(peakset_neut,peakset_mono)
neut_vs_tcel = compare_CRD(peakset_neut,peakset_tcel)

tcel_vs_mono = compare_CRD(peakset_tcel,peakset_mono)
tcel_vs_neut = compare_CRD(peakset_tcel,peakset_neut)


pdf("CRD_pairwise_comparisons_between_cell_types.pdf")
M = matrix(c(1,neut_vs_mono$fraction,M$fraction,mono_vs_neut$fraction,1,mono_vs_tcel$fraction,
             tcel_vs_neut$fraction,tcel_vs_mono$fraction,1),ncol=3,byrow=T)
colnames(M) = c("Neutrophils","Monocytes","T cells")
rownames(M) = c("Neutrophils","Monocytes","T cells")
corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()


