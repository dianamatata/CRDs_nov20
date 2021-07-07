# Section One ---------------------------------

# Section Two =================================

### Section Three #############################

# Clean environment ---------------------------------

rm(list=ls())
gc()


# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(corrplot)
color_palette=c("#00AFBB", "#E7B800", "#FC4E07","#972D15")


# Functions ---------------------------------

keep_shared_CRDs <- function(cell1_signif,cell1shared){
  crd_qtl_cell1_signif_shared=cell1_signif[(cell1_signif$V3 %in% cell1shared), ]
  crd_qtl_cell1_signif_shared
}


compute_overlap_nbr <- function(df_total, data){
  
  neut_vs_mono =(df_total %>% filter(data_type==data) %>% filter(cell1=='neut')  %>% filter(cell2=='mono'))$nbr_triplets_shared
  neut_vs_tcell=(df_total %>% filter(data_type==data) %>% filter(cell1=='neut')  %>% filter(cell2=='tcell'))$nbr_triplets_shared 
  mono_vs_neut=(df_total %>% filter(data_type==data) %>% filter(cell1=='mono')  %>% filter(cell2=='neut'))$nbr_triplets_shared
  mono_vs_tcell=(df_total %>% filter(data_type==data) %>% filter(cell1=='mono')  %>% filter(cell2=='tcell'))$nbr_triplets_shared
  tcell_vs_mono=(df_total %>% filter(data_type==data) %>% filter(cell1=='tcell')  %>% filter(cell2=='mono'))$nbr_triplets_shared
  tcell_vs_neut=(df_total %>% filter(data_type==data) %>% filter(cell1=='tcell')  %>% filter(cell2=='neut'))$nbr_triplets_shared
  M = matrix(c(1,as.integer(as.character(neut_vs_mono)[1]), as.integer(as.character(neut_vs_tcell)[1]),
               as.integer(as.character(mono_vs_neut)[1]),1,as.integer(as.character(mono_vs_tcell)[1]),
               as.integer(as.character(tcell_vs_neut)[1]),as.integer(as.character(tcell_vs_mono)[1]),1),
               ncol=3,byrow=T)
  M
}


compute_overlap_percent <- function(df_total, data){
  
  neut_vs_mono =(df_total %>% filter(data_type==data) %>% filter(cell1=='neut')  %>% filter(cell2=='mono'))$percent_triplets_shared 
  neut_vs_tcell=(df_total %>% filter(data_type==data) %>% filter(cell1=='neut')  %>% filter(cell2=='tcell'))$percent_triplets_shared 
  mono_vs_neut=(df_total %>% filter(data_type==data) %>% filter(cell1=='mono')  %>% filter(cell2=='neut'))$percent_triplets_shared
  mono_vs_tcell=(df_total %>% filter(data_type==data) %>% filter(cell1=='mono')  %>% filter(cell2=='tcell'))$percent_triplets_shared
  tcell_vs_mono=(df_total %>% filter(data_type==data) %>% filter(cell1=='tcell')  %>% filter(cell2=='mono'))$percent_triplets_shared
  tcell_vs_neut=(df_total %>% filter(data_type==data) %>% filter(cell1=='tcell')  %>% filter(cell2=='neut'))$percent_triplets_shared
  M = matrix(c(1,as.numeric(as.character(neut_vs_mono)[1]), as.numeric(as.character(neut_vs_tcell)[1]),
               as.numeric(as.character(mono_vs_neut)[1]),1,as.numeric(as.character(mono_vs_tcell)[1]),
               as.numeric(as.character(tcell_vs_neut)[1]),as.numeric(as.character(tcell_vs_mono)[1]),1),
             ncol=3,byrow=T)
  M
}


# Directories and Data ---------------------------------

dir_signif='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_signif' # here only triplets
dir_signif_pval='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/significants'
dir_all='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/all_triplets'
dir_shared='/Users/dianaavalos/Programming/A_CRD_plots/CRD_sharing/'
shared_crds=as.data.frame(data.table::fread(paste0(dir_shared,name,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE))
dir_BN='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/BN'
dir_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/plots'

cell_pairs=c(c('neut','mono'),c('neut','tcell'),c('mono','neut'),c('mono','tcell'),c('tcell','neut'),c('tcell','mono'))
module='mean'

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')

data_type='hist'
cell_type='neut'
pair=1
data='hist'

# Stats ---------------------------------

df_total=data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_total) <- c("data_type","cell1","cell2" ,"triplets_cell1_signif","triplets_cell1_signif_crd_shared",
                        "nbr_triplets_shared","percent_triplets_shared")

for(data_type in data_types){
  for (pair in c(1,3,5,7,9,11)){
    
    cell1=cell_pairs[pair]
    cell2=cell_pairs[pair+1]
    
    name=paste0(data_type,'_',module,'_',cell1,'_vs_',cell2)
    cat (name, '  ')
    
    # Load files =================================
    
    # load shared crd file  
    shared_crds=fread(paste0(dir_shared,name,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE)
    # load triplets from 2 cell types
    cell1_signif=fread(paste0(dir_signif,'/',data_type,'_',cell1,'_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)
    cell2_signif=fread(paste0(dir_signif,'/',data_type,'_',cell2,'_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)
    # load 2 BM
    BN_cell1=fread(paste0(dir_BN,'/','BN_',data_type,'_',cell1,'.txt'), head=FALSE, stringsAsFactors=FALSE)
    BN_cell2=fread(paste0(dir_BN,'/','BN_',data_type,'_',cell2,'.txt'), head=FALSE, stringsAsFactors=FALSE)
    colnames(BN_cell1)=colnames(BN_cell2)=c("name","nbr" ,"var","gene" ,"crd", "causal","reactive","independent","L1","L2","L3")
    BN_cell1$model <- colnames(BN_cell1[,6:8])[apply(BN_cell1[,6:8],1,which.max)]
    BN_cell2$model <- colnames(BN_cell2[,6:8])[apply(BN_cell2[,6:8],1,which.max)]
    
    # Algo =================================
    
    # filter triplets with CRD shared
    cell1_signif_shared=keep_shared_CRDs(cell1_signif,shared_crds$V1)
    cell2_signif_shared=keep_shared_CRDs(cell2_signif,shared_crds$V2)

    # compile all the shared significant triplets
    triplet_shared_cell1_with_cell2=data.frame()
    for (i in 1:dim(cell1_signif_shared)[1]) {
      var=cell1_signif_shared[i,]$V1
      gene=cell1_signif_shared[i,]$V2
      crd=cell1_signif_shared[i,]$V3
      
      cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == cell1_signif_shared[i,]$V3)] # many CRDs can correspond
      for (crds in cell2_equivalent_crd1){
        subset_triplets_cell2=cell2_signif_shared %>% filter(V1==var) %>% filter(V2==gene)  %>% filter(V3==crds) 
        subset_triplets_cell2$V4=crd
        if(dim(subset_triplets_cell2)[1]!=0) {
          triplet_shared_cell1_with_cell2=rbind(triplet_shared_cell1_with_cell2,subset_triplets_cell2)
        }
      }
    }
    
    colnames(triplet_shared_cell1_with_cell2)=c("var","gene","crd_c2","crd_c1")
    
    # write the shared significant triplets
    filename=paste0(dir_signif,'/',name,'_shared.txt')
    write.table(triplet_shared_cell1_with_cell2,file=filename)
    dim(triplet_shared_cell1_with_cell2)
    
    # check dimensions, and add to df_total
    df=cbind(data_type,
             cell1,
             cell2,
             dim(cell1_signif)[1],
             dim(cell1_signif_shared)[1],
             dim(triplet_shared_cell1_with_cell2)[1],
             dim(triplet_shared_cell1_with_cell2)[1]/dim(cell1_signif_shared)[1]*100
             )
    
    colnames(df)=c("data_type","cell1","cell2" ,"triplets_cell1_signif","triplets_cell1_signif_crd_shared",
                   "nbr_triplets_shared","percent_triplets_shared")
    df_total=rbind(df_total,df)
   
    
    # compare the 2 Bayesian models:
    triplet_shared_models=data.frame()
    for (t in 1:dim(triplet_shared_cell1_with_cell2)[1]) {
      triplet=triplet_shared_cell1_with_cell2[t,]
      print(triplet)
      triplet1=BN_cell1 %>% filter(var==triplet$var) %>% filter(gene==triplet$gene)  #%>% filter(crd==triplet$crd_c1) 
      triplet2=BN_cell2 %>% filter(var==triplet$var) %>% filter(gene==triplet$gene)  #%>% filter(crd==triplet$crd_c2) 
      model1=triplet1$model
      model2=triplet2$model
      # once that we have the 2 models, create list to check how many switch
      triplet_shared_models=rbind(triplet_shared_models,cbind(model1,model2))
    }
    # write in file
    filename=paste0(dir_signif,'/models_triplets_switch_',name,'.txt')
    write.table(triplet_shared_models,file=filename)

  }
}

# write the triplet sharing
filename=paste0(dir_signif,'/stats_all_triplets_shared.txt')
write.table(df_total,file=filename)



# Change of model in triplets shared, graphs ---------------------------------

for(data_type in data_types){
  all.files <- intersect(list.files(path=dir_signif, pattern='models_triplets_switch_', full.names=TRUE, recursive=FALSE),
                       list.files(path=dir_signif, pattern=data_type, full.names=TRUE, recursive=FALSE))
  
  for (file in all.files){
    f=paste(unlist(strsplit(unlist(strsplit(file, "/"))[8], "_"))[4:8], collapse = '_')
    BN1=fread(file, head=FALSE, stringsAsFactors=FALSE)[,2:3]
    colnames(BN1)=c("m1","m2")
    BN1 <- BN1 %>% arrange(m1, m2)
    BNsummary=table(BN1)
    BNmolten = melt(BNsummary,id.vars="m1")
    target <- c("causal","reactive","independent")
    if(!all(target %in% BNmolten$m1 == TRUE)){
      BNmolten=bind_rows(BNmolten, data.frame(m1=target[!target %in% BNmolten$m1],m2=target[!target %in% BNmolten$m1],value=0))
    }
    if(!all(target %in% BNmolten$m1 == TRUE)){
      BNmolten=bind_rows(BNmolten, data.frame(m1=target[!target %in% BNmolten$m2],m2=target[!target %in% BNmolten$m2],value=0))
    }
    BNmolten <- BNmolten %>% mutate(m1 = fct_relevel(m1, "causal","reactive","independent"))     # Reorder following a precise order
    
    # add percentage
    
    B_percent=data.frame()
    for(model in c("causal","reactive","independent")){
      cat(model)
      B_sub=BNmolten %>% filter (m1==model)
      B_sub$percent <- B_sub$value / sum(B_sub$value) * 100
      B_percent=rbind(B_percent,B_sub)
    }
    
    BNmolten=B_percent
    
    # stacked bar plot =================================
    
    pdf(paste0(dir_out,"/Switch_models_",f,".pdf"),paper='a4r',7,5)
    a  <- ggplot(BNmolten, aes(fill=m2, y=value, x=m1 ,label = scales::percent(round(value,digits=2))))+ 
      geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
      labs(x = "Model in cell 1",y = "Number of triplets switching models", title="Switch in model for shared triplets between cell types", fill='Model in cell 2') +
      geom_text(aes(label = paste0(round(value,digits=2))), position = position_stack(vjust = .5),size = 3 ) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(a)
    dev.off()
    
    # same with no title
    
    pdf(paste0(dir_out,"/Switch_models_notitles_",f,".pdf"),paper='a4r',7,5)
    
    a  <- ggplot(BNmolten, aes(fill=m2, y=value, x=m1 ,label = scales::percent(round(value,digits=2))))+ 
      geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
      labs(x = "Model in cell 1",y = "Number of triplets switching models", fill='Model in cell 2') +
      geom_text(aes(label = paste0(round(value,digits=2))), position = position_stack(vjust = .5),size = 3 ) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(a)
    
    dev.off()
    
    # stacked bar plot percent =================================
  
    
    pdf(paste0(dir_out,"/Switch_models_perc_",f,".pdf"),paper='a4r',7,5)
    a  <- ggplot(BNmolten, aes(fill=m2, y=percent, x=m1 ,label = scales::percent(round(percent,digits=2))))+ 
      geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
      labs(x = "Model in cell 1",y = "Percentage of triplets switching models", title="Switch in model for shared triplets between cell types", fill='Model in cell 2') +
      geom_text(aes(label = paste0('n:', round(value,digits=2), "\n",round(percent,digits=2),'%')), position = position_stack(vjust = .5),size = 3 ) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(a)
    dev.off()
        
  }
}


  

# Compare with intersection 3 ---------------------------------
# take neut intersect mono, add tcell, look at intersect
for(data_type in data_types){
    cell1='neut'
    cell2='mono'
    cell3='tcell'
    
    filename=paste0("/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_signif/",data_type,"_mean_",cell1,"_vs_",cell2,"_shared.txt")
    triplet_shared_c1_c2=fread(filename, head=FALSE, stringsAsFactors=FALSE)[,2:4]
    shared_crds=fread(paste0(dir_shared,data_type,'_',module,'_',cell1,'_vs_',cell3,'_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE)
    # load triplets from 2 cell types
    cell3_signif=fread(paste0(dir_signif,'/',data_type,'_',cell2,'_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)
    # filter triplets with CRD shared
    cell3_signif_shared=keep_shared_CRDs(cell3_signif,shared_crds$V2)
    dim(cell3_signif_shared)
    # compile all the shared significant triplets
    triplet_shared_c1_c2_c3=data.frame()
    
    for (i in 1:dim(cell3_signif_shared)[1]) {
      var=cell3_signif_shared[i,]$V1
      gene=cell3_signif_shared[i,]$V2
      crd=cell3_signif_shared[i,]$V3
      
      cell1_equivalent_crd3 = shared_crds$V1[which(shared_crds$V2 == crd)] # many CRDs can correspond
      for (crds in cell1_equivalent_crd3){
        print(length(crds))
        subset_triplets_cell3=triplet_shared_c1_c2 %>% filter(V2==var) %>% filter(V3==gene)  %>% filter(V4==crds) 
        if(dim(subset_triplets_cell3)[1]!=0) {
          triplet_shared_c1_c2_c3=rbind(triplet_shared_c1_c2_c3,subset_triplets_cell3)
        }
      }
    }
    print(triplet_shared_c1_c2_c3)
    
}


# Triplet sharing ---------------------------------

for(data_type in data_types){
  Mnbr=compute_overlap_nbr(df_total, data_type)
  Mpercent=compute_overlap_percent(df_total, data_type)

  pdf(paste0(dir_out,data_type,"_Mnbr.pdf"))
  colnames(Mnbr) = c("NEU","MON","TCL")
  rownames(Mnbr) = c("NEU","MON","TCL")
  corrplot(Mnbr,is.corr=F,cl.lim = c(0, 150),p.mat = Mnbr,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
  pdf(paste0(dir_out,data_type,"_Mpercent.pdf"))
  colnames(Mpercent) = c("NEU","MON","TCL")
  rownames(Mpercent) = c("NEU","MON","TCL")
  corrplot(Mpercent,is.corr=F,cl.lim = c(0, 100),p.mat = Mpercent,sig.level=-1,insig = "p-value",number.cex=1.5)
  dev.off()
  
}


# Sizes of Triplets per cell type ---------------------------------

# load triplets
NEU=fread(paste0(dir_signif,'/',data_type,'_','neut','_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)
MON=fread(paste0(dir_signif,'/',data_type,'_','mono','_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)
TCL=fread(paste0(dir_signif,'/',data_type,'_','tcell','_triplet.txt'), head=FALSE, stringsAsFactors=FALSE)

dim(NEU)[1]
dim(MON)[1]
dim(TCL)[1]

df_total=fread(paste0(dir_signif,'/stats_all_triplets_shared.txt'), head=FALSE, stringsAsFactors=FALSE)

# compute intersection sizes ---------------------------------

# filter triplets with CRD shared, naybe better to have a loop
shared_neu_mon=fread(paste0(dir_shared,data_type,'_',module,'_','neut','_vs_','mono','_sharedCRDs.txt'), head=FALSE, stringsAsFactors=FALSE)
dim(shared_neu_mon)
dim(shared_mon_neu)
cell1_signif_shared=keep_shared_CRDs(cell1_signif,shared_crds$V1)

cell2_signif_shared=keep_shared_CRDs(cell2_signif,shared_crds$V2)
signif_shared_neut=keep_shared_CRDs(cell1_signif,shared_crds$V1)
cell2_signif_shared=keep_shared_CRDs(cell2_signif,shared_crds$V2)

### to be written

N_M = length(intersect(NEU,MON))
N_T = length(intersect(NEU,TCL))
M_T = length(intersect(TCL,MON))
N_M_T = length(intersect(NEU,intersect(TCL,MON)))

NEU1 = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
MON1 = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
TCL1 = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T)

toplot2 = data.frame(NEU=c(length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+N_T-N_M_T+0,N_M_T),
                     MON=c(length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+0+M_T-N_M_T,N_M_T),
                     TCL=c(length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T),0+N_T-N_M_T+M_T-N_M_T,N_M_T),Group =c("1cell","2cells","3cells"))


toplot2.molten = melt(toplot2,id.vars="Group")
colnames(toplot2.molten)[2] = "CellType"
toplot2.molten$Group <- factor(toplot2.molten$Group ,levels = c("1cell","2cells","3cells"))

# plot intersection sizes ---------------------------------

expressionInput <- c(NEU = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), MON = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), 
                     TCL = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T), 
                     `NEU&MON` = N_M-N_M_T, `NEU&TCL` = N_T-N_M_T,`MON&TCL` = M_T-N_M_T, `NEU&MON&TCL` = N_M_T)

pdf(paste0(out_dir_shared,name,"_Overlap_chromatin_peaks.pdf"),paper='a4r')
#upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
dev.off()


    # check the q value, is triplet still significant
    head(cell1_signif_shared_equiv)
    
    # check if the triplet switched to another model
    # find the right triplet/line, substitute the crd equiv and check the winning model in the other one
    triplets_BN_cell1=fread(paste0(dir_BN,'/','BN_',data_type,'_',cell1,'.txt'), head=FALSE, stringsAsFactors=FALSE)
    triplets_BN_cell2=fread(paste0(dir_BN,'/','BN_',data_type,'_',cell2,'.txt'), head=FALSE, stringsAsFactors=FALSE)
    head(triplets_BN_cell1)
    colnames(triplets_BN_cell1)=colnames(triplets_BN_cell2)=c("name","nbr" ,"var","gene" ,"crd", "causal","reactive","independent","L1","L2","L3")


 ######
  for (i in 1:length(crd_qtl_cell1_signif_shared[,1])) {
    # cat(i,' ')

    # find equiv CRDs in cell2 + variant
    cell2_equivalent_crd1 = shared_crds$V2[which(shared_crds$V1 == crd_qtl_cell1_signif_shared[i,]$phenotype_ID)] # many CRDs can correspond
    variant = crd_qtl_cell1_signif_shared[i,]$var_ID
    # filter
    sub1=crd_qtl_cell2_all_shared %>% filter(phenotype_ID %in%  cell2_equivalent_crd1)   
    subset=sub1 %>% filter(var_ID %in%  variant) 
    subset$nom_pval
    
    # save values
    
    if (length(subset[,1]) != 0){
      line=paste0(subset$nom_pval)
      line_bis=paste0(subset$nom_pval, ' ',subset$var_ID, ' ',subset$phenotype_ID)
      write(line_bis,file=filename,append=TRUE)
    }
  }

