# Clean environment
rm(list=ls())
gc()

#############################################################################################
#
# PACKAGES
#
#############################################################################################

library(qvalue)
library(ggplot2)
library(gplots)
library(data.frame)
library(GenomicRanges)
library(RColorBrewer)
library("circlize")
library(igraph)


#############################################################################################
#
# DIRECTORIES AND FILES
#
#############################################################################################

directory='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:significant/'
out_dir='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:other_plots/'

files <- list.files(path=directory, pattern="0.0*.txt", full.names=TRUE, recursive=FALSE)


#############################################################################################
#
# MAIN
#
#############################################################################################

Paper=c("#00AFBB", "#E7B800", "#FC4E07" )
mycolors= Paper

for(data_type in c('hist','methyl')){
  for(module in c('mean','loom')){
    for(FDR in c('0.01','0.05')){
      
      files <- intersect(list.files(path=directory, pattern=data_type, full.names=TRUE, recursive=FALSE),
                         intersect(list.files(path=directory, pattern=module, full.names=TRUE, recursive=FALSE),
                                   list.files(path=directory, pattern=FDR, full.names=TRUE, recursive=FALSE)))

      transCRD_NEU = as.data.frame(data.table::fread(files[grepl("neu", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      transCRD_MON = as.data.frame(data.table::fread(files[grepl("mon", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      transCRD_TCL = as.data.frame(data.table::fread(files[grepl("tc", files)==TRUE], head=TRUE, stringsAsFactors=FALSE))
      colnames(transCRD_NEU) = colnames(transCRD_MON) = colnames(transCRD_TCL) = c("idx1","chr1","midplace","id1","idx2","chr2","midplace2","id2","corr","pval","qvalue")
      name=paste0(data_type,'_',module,'_',FDR)
      cat (name)
      try(plot_Connectivity_ALL_CellTypes(transCRD_NEU,transCRD_MON,transCRD_TCL, name))
    }
  }
}




data_type="hist"
module='mean'
FDR='0.05'

## need for the 3 cell types of 1 condition, the trans CRD association id1 and id2

#############################################################################################
#
# Fig 5C FUNCTION
#
#############################################################################################

plot_Connectivity_ALL_CellTypes <- function(transCRD_NEU,transCRD_MON,transCRD_TCL, name){
  
  pdf(paste0(out_dir,"5c_Connectivity_ALL_CellTypes_",name,".pdf"), 6, 6)
  hist.NEU = hist(table(c(transCRD_NEU$id1, transCRD_NEU$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  hist.MON = hist(table(c(transCRD_MON$id1, transCRD_MON$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  hist.TCL = hist(table(c(transCRD_TCL$id1, transCRD_TCL$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
  
  frac.NEU = hist.NEU$counts/sum(hist.NEU$counts)*100
  frac.MON = hist.MON$counts/sum(hist.MON$counts)*100
  frac.TCL = hist.TCL$counts/sum(hist.TCL$counts)*100
  
  toplot = data.frame(NEU = frac.NEU,MON = frac.MON,TCL = frac.TCL,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  toplot.melted = reshape2::melt(toplot,id.vars = "Number", measure.vars = c("NEU","MON","TCL"))
  colnames(toplot.melted)[2] = "CellType"
  toplot.melted$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
  
  ### now plot Fig 5C
  g <- ggplot(toplot.melted, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
    geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07" )) +
    labs(x = "Connectivity",y = "Fraction of CRDs (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
  
  print(g)
  dev.off()
}


### then not finished, would like to have 5e with PCHIC score and fraction of associations too
### problem. we need start and stop for Trans CRD associations and so far we have chr and middle point?
### need to integrate start and stop from another file containing these CRDs. the one that is inputed to compute the trans associations

#############################################################################################
#
# GET HI-C DATA
#
#############################################################################################

#### which ones for tcells? CD4 not CD8!!!
PCHiC = fread('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/PCHiC_peak_matrix_cutoff5.tsv')
PCHiC_Neu=PCHiC$Neu
PCHiC_Mon=PCHiC$Mon
PCHiC_Tcl=PCHiC$tCD4

#######  get PCHiC data

colnames(PCHiC)[1] = "baitChr"
interchromosomal = which(PCHiC$baitChr!=PCHiC$oeChr)
PCHiC = PCHiC[interchromosomal,]
baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
CRD1bed <- GRanges(seqnames=trans_data$chr1,ranges=IRanges(start=trans_data$start1, end=trans_data$end1))
CRD2bed <- GRanges(seqnames=trans_data$chr2,ranges=IRanges(start=trans_data$start2, end=trans_data$end2))

#fwd
x = findOverlaps(baitbed,CRD1bed)
y = findOverlaps(oebed,CRD2bed)
tmp = rbind(as.trans_data.frame(x),as.trans_data.frame(y))
validated.fwd = tmp[which(duplicated(tmp)),]

#bwd
x = findOverlaps(baitbed,CRD2bed)
y = findOverlaps(oebed,CRD1bed)
tmp = rbind(as.trans_data.frame(x),as.trans_data.frame(y))
validated.bwd = tmp[which(duplicated(tmp)),]

validated = unique(rbind(validated.fwd,validated.bwd))

hic_validated = rep(1,nrow(PCHiC))

for(i in 1:nrow(validated)){
  currenthic = validated$queryHits[i]
  currentqval = trans_data$qval[validated$subjectHits[i]]
  if(currentqval<hic_validated[currenthic]){
    hic_validated[currenthic] = currentqval
  }
}



#############################################################################################
#
# OTHERS
#
#############################################################################################


myhist_bg.NEU = hist(PCHiC$Neu,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.NEU = hist(PCHiC$Neu[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)

myhist_bg.MON = hist(PCHiC$Mon,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.MON = hist(PCHiC$Mon[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)

myhist_bg.TCL = hist(PCHiC$nCD4,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.TCL = hist(PCHiC$nCD4[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)




####* loop around files and label cell types

path = "/Users/dianaavalos/Programming/A_CRD_plots/trans_files"
file.names <- dir(path, pattern =".txt.gz")
for(i in 1:length(file.names)){
  cat(file.names[i], '  ')
  filename=file.names[i]
  name= substr(file.names[i],1,nchar(file.names[i])-11)
  
  if (str_detect(name, 'tcell')){
    PCHiC_cell=PCHiC_Tcl
  }
  if (str_detect(name, 'mono')){
    PCHiC_cell=PCHiC_Mon
  }
  if (str_detect(name, 'neut')){
    PCHiC_cell=PCHiC_Neu
  }
}


for (f in files){
  f=files[1]
  file=basename(f)
  name=paste0(str_sub(file, 1, - 27) ,str_sub(file, 34, - 5))
  
  TRH_signif = as.data.frame(data.table::fread(f, head=TRUE, stringsAsFactors=FALSE))
  trans_data_q001=TRH_signif
}


  