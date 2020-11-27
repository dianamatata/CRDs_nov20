library(qvalue)
library(ggplot2)
library(gplots)
library(data.table)
library(GenomicRanges)
library(RColorBrewer)
library("circlize")
library(igraph)

# keep DATA when qvalue < 0.01
DATA = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002670_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))
colnames(DATA) = c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval")
Q = qvalue(DATA$pval)
DATA$qval = Q$qvalues
DATAs = DATA[Q$qvalue < 0.01, ]

#WRITE SIGNIFICANT HITS
write.table(DATAs, "test.significant1FDR.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')
colnames(PCHiC)[1] = "baitChr"
interchromosomal = which(PCHiC$baitChr!=PCHiC$oeChr)
PCHiC = PCHiC[interchromosomal,]

baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))

CRD1bed <- GRanges(seqnames=DATA$chr1,ranges=IRanges(start=DATA$start1, end=DATA$end1))
CRD2bed <- GRanges(seqnames=DATA$chr2,ranges=IRanges(start=DATA$start2, end=DATA$end2))

#fwd
x = findOverlaps(baitbed,CRD1bed)
y = findOverlaps(oebed,CRD2bed)
tmp = rbind(as.data.frame(x),as.data.frame(y))
validated.fwd = tmp[which(duplicated(tmp)),]

#bwd
x = findOverlaps(baitbed,CRD2bed)
y = findOverlaps(oebed,CRD1bed)
tmp = rbind(as.data.frame(x),as.data.frame(y))
validated.bwd = tmp[which(duplicated(tmp)),]

validated = unique(rbind(validated.fwd,validated.bwd))

hic_validated = rep(1,nrow(PCHiC))

for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentqval = DATA$qval[validated$subjectHits[i]]
    if(currentqval<hic_validated[currenthic]){
        hic_validated[currenthic] = currentqval
    }
}


myhist_bg = hist(PCHiC$Neu,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif = hist(PCHiC$Neu[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)

pdf("HiC_validation.pdf",5,5)
toplot = data.frame(counts = myhist_signif$counts/myhist_bg$counts*100,Number = c("0-10","11-20","21-50","51-100",">100"))
toplot$Number = factor(toplot$Number,levels = c("0-10","11-20","21-50","51-100",">100"))
g <- ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("HiC contacts with CRD associations") +
  geom_bar(stat = "identity",fill="#E69F00") +
  labs(x = "PC HiC Score",y = "Fraction (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
print(g)
dev.off()


# annotateHiC <- function(x, celltype= "Neu"){
#     if(celltype == "Neu"){
#         #forward
#         DATA.tmp = subset(DATA, chr1 == x$baitChr & ((start1<=x$baitEnd & start1>=x$baitStart) | (end1<=x$baitEnd & end1>=x$baitStart)))
#         interactiondata.fwd = subset(DATA.tmp, chr2 == x$oeChr & ((start2<=x$oeEnd & start2>=x$oeStart) | (end2<=x$oeEnd & end2>=x$oeStart)))
#         #backward
#         DATA.tmp = subset(DATA, chr1 == x$oeChr & ((start1<=x$oeEnd & start1>=x$oeStart) | (end1<=x$oeEnd & end1>=x$oeStart)))
#         interactiondata.bwd = subset(DATA.tmp, chr2 == x$baitChr & ((start2<=x$baitEnd & start2>=x$baitStart) | (end2<=x$baitEnd & end2>=x$baitStart)))
#     }
#     interactiondata = rbind(interactiondata.fwd,interactiondata.bwd)
#     if(nrow(interactiondata)>=1){
#         n = nrow(interactiondata)
#         mean = mean(interactiondata$Neu)
#     } else {
#         n = 0
#         mean = 0
#     }
#     c(n,mean)
# }
#
#
# #TO DO: change the direction of the analysis, annotate HiC data with correlation data (all not only the significant bit)
# getHiCcontact <- function(x, celltype= "Neu"){
#     CRD1Chr = as.character(x[2])
#     CRD1Start = as.numeric(x[3])
#     CRD1End = as.numeric(x[4])
#     CRD2Chr = as.character(x[7])
#     CRD2Start = as.numeric(x[8])
#     CRD2End = as.numeric(x[9])
#     if(celltype == "Neu"){
#         #forward
#         promoterdata = subset(PCHiC, baitChr == CRD1Chr & ((baitStart<=CRD1End & baitStart>=CRD1Start) | (baitEnd<=CRD1End & baitEnd>=CRD1Start)) & Neu >5)
#         interactiondata.fwd = subset(promoterdata, oeChr == CRD2Chr & ((oeStart<=CRD2End & oeStart>=CRD2Start) | (oeEnd<=CRD2End & oeEnd>=CRD2Start)))
#         #backward
#         promoterdata = subset(PCHiC, baitChr == CRD2Chr & ((baitStart<=CRD2End & baitStart>=CRD2Start) | (baitEnd<=CRD2End & baitEnd>=CRD2Start)) & Neu >5)
#         interactiondata.bwd = subset(promoterdata, oeChr == CRD1Chr & ((oeStart<=CRD1End & oeStart>=CRD1Start) | (oeEnd<=CRD1End & oeEnd>=CRD1Start)))
#     }
#     interactiondata = rbind(interactiondata.fwd,interactiondata.bwd)
#     if(nrow(interactiondata)>=1){
#         n = nrow(interactiondata)
#         mean = mean(interactiondata$Neu)
#     } else {
#         n = 0
#         mean = 0
#     }
#     c(n,mean)
# }
#
#
# stop()
#
# hicontacts = t(apply(mapdata,1,getHiCcontact))

DATAs$chr1 = paste0("chr",DATAs$chr1)
DATAs$chr2 = paste0("chr",DATAs$chr2)


#PLOT1: HISTOGRAM + PI1
pdf("Histogram.pdf",10,5)
par(mfrow=c(1, 2))
hist(DATA$pval, xlab="Nominal P-values", main="", breaks=100)
abline(h=Q$pi0 * nrow(DATA) / 100, col="red")
legend("topright", legend=c(paste("pi1=", signif(100-Q$pi0*100, 3), "%", sep=""), paste("#hits=", nrow(DATAs), " (1% FDR)", sep="")), bty="n")
hist(DATAs$corr, breaks=100, main="", xlab="Correlation coefficient")
hist(Q$qvalue,xlab="Q-values", main="", breaks=100)
dev.off()

#PLOT2: EXAMPLE HEATMAP
#TD = as.data.frame(data.table::fread(paste("zcat ~/VitalIT/SGX/V2/data_correlation/trans/plot/LCL_ALL.",C1,".",C2,".plot.txt.gz", sep="")), head=FALSE)
#M1 = as.data.frame(data.table::fread(paste("zcat ~/VitalIT/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr",C1,".subset.txt.gz", sep="")), head=FALSE)
#M2 = as.data.frame(data.table::fread(paste("zcat ~/VitalIT/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr",C2,".subset.txt.gz", sep="")), head=FALSE)

#C1=10; C2=16; C1c=c(4000, 4500); C2c=c(3250, 3750);
#TDs = TD[TD$V1 > C1c[1] & TD$V1 < C1c[2] & TD$V3 > C2c[1] & TD$V3 < C2c[2], ]
#M1s = M1[M1$V1 > C1c[1] & M1$V5 < C1c[2], ]
#M2s = M2[M2$V1 > C2c[1] & M2$V5 < C2c[2], ]

#pdf("figureB.4.2.pdf", 8,8)
#plot(TDs$V1-C1c[1], TDs$V3-C2c[1], xlim=c(-125, 500), ylim=c(-125, 500), xaxt="n", yaxt="n", xlab=paste("Chromosome", C1), ylab=paste("Chromosome", C2), col=rgb(0,0,1, TDs$V5^2), pch=15, main="", cex=0.5)
#abline(h=0); abline(v=0);
#points((M1s$V1 + M1s$V5)/2-C1c[1], (M1s$V1-M1s$V5)*0.5, col=rgb(0,0,1, M1s$V10^2), pch=15, cex=0.5)
#points((M2s$V1-M2s$V5)*0.5, (M2s$V1 + M2s$V5)/2-C2c[1], col=rgb(0,0,1, M2s$V10^2), pch=15, cex=0.5)
#axis(1, at=seq(0, 500, 100), labels=seq(C1c[1], C1c[2], 100))
#axis(2, at=seq(0, 500, 100), labels=seq(C2c[1], C2c[2], 100))
#dev.off()

#PLOT3: CIRCLE PLOT FOR BEST 1000 LINKS
COL = brewer.pal(9,"Set1")
DATAss = DATAs[order(DATAs$pval),]
DATAss = DATAss[1:1000, ]

pdf("Circle_plot.pdf", 8, 8)
par(mar=c(1, 1, 1, 1))
bed1 = data.frame(chr=DATAss$chr1, from=DATAss$start1, to=DATAss$end1)
bed2 = data.frame(chr=DATAss$chr2, from=DATAss$start2, to=DATAss$end2)
dir = (DATAss$corr > 0)
circos.clear()
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 22:1))
circos.genomicLink(bed1, bed2, col=rgb(0,0,1,0.05))
dev.off()

#PLOT4: CONNECTIVITY
pdf("Connectivity.pdf", 6, 6)
hist(table(c(DATAs$id1, DATAs$id2)), breaks=40, xlab="Number of connected modules (module degree)", main="")
conhist = hist(table(c(DATAs$id1, DATAs$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
frac = conhist$counts/sum(conhist$counts)*100
barplot(frac,names = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
toplot = data.frame(counts = frac,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
toplot$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("Connectivity of CRD trans associations") +
  geom_bar(stat = "identity",fill="#E69F00") +
  labs(x = "Connectivity",y = "Fraction (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))

dev.off()

#PLOT4: CONNECTIVITY 3 
pdf("Connectivity.pdf", 6, 6)
hist(table(c(DATAs.NEU$id1, DATAs.NEU$id2)), breaks=40, xlab="Number of connected modules (module degree)", main="")
conhistNEU = hist(table(c(DATAs.NEU$id1, DATAs.NEU$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
conhistMON = hist(table(c(DATAs.MON$id1, DATAs.MON$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
conhistTCL = hist(table(c(DATAs.TCL$id1, DATAs.TCL$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)

fracNEU = conhistNEU$counts/sum(conhistNEU$counts)*100
fracMON = conhistMON$counts/sum(conhistMON$counts)*100
fracTCL = conhistTCL$counts/sum(conhistTCL$counts)*100

barplot(frac,names = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
toplot = data.frame(counts = frac,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
toplot$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
ggplot(toplot, aes(x = Number, y = counts))+ ggtitle("Connectivity of CRD trans associations") +
  geom_bar(stat = "identity",fill="#E69F00") +
  labs(x = "Connectivity",y = "Fraction (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))

dev.off()


#PLOT5: NETWORK EXAMPLE OF THE 1000 BEST LINKS
library('igraph')
library(RColorBrewer)
COL1 = brewer.pal(8,"Set1")
COL2 = brewer.pal(8,"Set2")
COL3 = brewer.pal(8,"Set3")
COL=c(COL1, COL2, COL3)
REF = data.frame(i=1:22, s=paste("", 1:22, sep=""))
# DATA with qvalue < 0.01, ordered by pvalue, take first 1000 nodes
DATAss = DATAs[order(DATAs$pval),]
DATAss = DATAss[1:1000, ]

LINKS = data.frame(from=DATAss$id1, to=DATAss$id2, weigth=ifelse(DATAss$corr>0, 1, 2), stringsAsFactors=FALSE)
NODES = data.frame(id=names(table(c(LINKS$from, LINKS$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS$from, LINKS$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
DATAnet = graph_from_data_frame(d=LINKS, vertices=NODES, directed=FALSE)
deg <- degree(DATAnet, mode="all")

# added by me
write.table(LINKS, "LINKS.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

V(DATAnet)$size <- log2(deg)+1
V(DATAnet)$frame.color <- "white"
V(DATAnet)$color <- COL[match(V(DATAnet)$chr, REF$s)]
V(DATAnet)$border <- "n"
V(DATAnet)$label <- ""
E(DATAnet)$arrow.mode <- 0
E(DATAnet)$color <- ifelse(E(DATAnet)$weigth == 1, rgb(0,0,1, 0.4), rgb(1,0,0, 0.4))

pdf("igraph_network.pdf", 8, 8)
par(mar=c(1,1,1,1))
plot(DATAnet)
legend("topleft", fill=COL[1:22], legend=REF$s, bg="white", title="Nodes:", ncol=11, cex=0.8)
legend("bottomleft", fill=c(rgb(0,0,1, 0.4), rgb(1,0,0, 0.4)), legend=c("Positively correlated", "Negatively correlated"), bg="white", title="Links:")
dev.off()

###### PLOT6: CALL TRHs using greedy algorithm

LINKS.TRH = data.frame(from=DATAs$id1, to=DATAs$id2, weigth=ifelse(DATAs$corr>0, 1, 2), stringsAsFactors=FALSE)
NODES.TRH = data.frame(id=names(table(c(LINKS.TRH$from, LINKS.TRH$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS.TRH$from, LINKS.TRH$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
DATAnet.TRH = graph_from_data_frame(d=LINKS.TRH, vertices=NODES.TRH, directed=FALSE)
communities  = fastgreedy.community(DATAnet.TRH)
N_TRH = max(communities$membership)
N_TRH_TOSHOW = 50
pdf("TRH_group_sizes.pdf")
df = data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
df <-df[order(-df$SIZE),]
df <- data.frame(head(df$SIZE,50),seq(1,50))
names(df)[1]='SIZE'
names(df)[2]='ID'

p2<- ggplot(df, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") 
p2 + theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() + scale_y_continuous(trans = 'log10')

#head(df,10)

#### g version
pdf("TRH_group_sizes.pdf")
p<- ggplot(df,aes(x = ID, y = SIZE))
p+labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") +
  geom_bar(stat="identity", fill="steelblue", width=0.7)+
  theme_minimal() + xlim(0, N_TRH_TOSHOW) + scale_y_continuous(trans = 'log10') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))
print(p)
dev.off()



#PLOT7: CALL chromosome pair frequencies

DATAsig = DATA[DATA$pval<1e-06,]
chromFreq = matrix(0,22,22)
for(i in 1:21){
    for(j in (i+1):22){
        frac_signif = sum(DATAsig$chr1 == i & DATAsig$chr2 == j)/sum(DATA$chr1 == i & DATA$chr2 == j)
        chromFreq[i,j] = frac_signif
        chromFreq[j,i] = frac_signif
    }
}

# chromFreq = sweep(chromFreq,MARGIN=1,FUN="/",STATS=rowSums(chromFreq))
coul <- colorRampPalette(brewer.pal(9, "Blues"))(15)

pdf("Frequency_by_chrom_pairs.pdf")
heatmap.2(chromFreq,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col=coul)
dev.off()

save(PCHiC,DATA,DATAs,DATAss,hic_validated,Q,file="TRANS_analysis.rda")
