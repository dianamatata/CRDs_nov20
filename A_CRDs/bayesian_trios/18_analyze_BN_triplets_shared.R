#### plot the distribution of the probability of the most likely model
## this will let you know whether when you call a model causal, reactive or independent, whether the probability is High or not


## need to rerun with new data: cutoff of bayesian proba at 0,5

# Clean environment ---------------------------------

rm(list=ls())
gc()

# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)


# Directories and Data ---------------------------------

# directories for cluster =================================
dir_BN='/home/users/a/avalosma/scratch/12_TRIPLETS/BN'
dir_BN_union='/home/users/a/avalosma/scratch/12_TRIPLETS/BN_union'
dir_BN_union_out='/home/users/a/avalosma/scratch/12_TRIPLETS/BN_union_plots_stats'

# directories for Mac =================================
dir_BN='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/BN'
dir_BN_union='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/triplets_shared_mixed_CRDs/BN_union'
dir_BN_union_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/plots2'


# debug =================================
data_type='hist'
cell_type_ref1='neut'
cell_type_query2='mono'
cell_type_ref1='mono'
cell_type_query2='neut'
# Main ---------------------------------

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')
color_palette=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')

for(data_type in data_types){
  for(cell_type_ref1 in cell_types){
    for(cell_type_query2 in cell_types){
      if (cell_type_ref1 != cell_type_query2){
        
        triplet_shared_models=data.frame()
        
        name_condition=paste0(data_type,'_',cell_type_query2,'_vs_',cell_type_ref1,'_all' )
        cat(name_condition, '  \n')
        # cell_type_query2 is cell type 2
        
        BN_one_cell=fread(paste0(dir_BN,'/',data_type,'_',cell_type_ref1,'_1000000.txt'), head=FALSE, stringsAsFactors=FALSE)
        BN_two_cells=fread(paste0(dir_BN_union,'/',data_type,'_',cell_type_query2,'_vs_',cell_type_ref1,'_all.txt'), head=FALSE, stringsAsFactors=FALSE)
        colnames(BN_one_cell)=c("name","nbr" ,"var","gene" ,"crd", "dist","AF","causal","reactive","independent","L1","L2","L3")
        colnames(BN_two_cells)=c("name","nbr" ,"var","gene" ,"crd", "dist","AF","causal","reactive","independent","L1","L2","L3")
        
        BN_one_cell$max1 <- pmax(BN_one_cell$causal,BN_one_cell$reactive,BN_one_cell$independent)
        BN_one_cell=BN_one_cell %>% filter(max1>0.5)
        BN_one_cell$model <- colnames(BN_one_cell[,8:10])[apply(BN_one_cell[,8:10],1,which.max)]
        
        BN_two_cells$max1 <- pmax(BN_two_cells$causal,BN_two_cells$reactive,BN_two_cells$independent)
        BN_two_cells=BN_two_cells %>% filter(max1>0.5)
        BN_two_cells$model <- colnames(BN_two_cells[,8:10])[apply(BN_two_cells[,8:10],1,which.max)]

        # for all the shared triplets, look at what are the 2 models and save in triplet_shared_models
        for(i in 1:nrow(BN_one_cell)) {
          triplet_cell1 <- BN_one_cell[i,]
          triplet_cell2=BN_two_cells %>% filter(var==triplet_cell1$var) %>% filter(gene==triplet_cell1$gene) %>% filter(crd==triplet_cell1$crd)
          triplet_cell2=triplet_cell2[1]
          if(dim(triplet_cell2)[1]>0){
            triplet_shared_models=rbind(triplet_shared_models,cbind(triplet_cell1$model,triplet_cell2$model))
          }
        }
        cat(paste0(f, " ",dim(triplet_shared_models), '\n'))

        
        filename=paste0(dir_BN_union_out,'/',name_condition,'_models_triplets_shared.txt')
        write.table(triplet_shared_models, filename, quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
      }
    }
  }
}

# results where wrong because triplicates in the --- vs ---- so we redo the plots
### plot the % of model change and the % of remains the same`?`



# Change of model in triplets shared, graphs ---------------------------------
all.files <- list.files(path=dir_BN_union_out, pattern='txt', full.names=TRUE, recursive=FALSE)

for (file in all.files){
    f=unlist(strsplit(unlist(strsplit(file, "/")[1])[8], ".txt")[1])
    BN1=fread(file, head=FALSE, stringsAsFactors=FALSE)
    colnames(BN1)=c("m1","m2")
    BN1 <- BN1 %>% arrange(m1, m2)
    BNsummary=table(BN1)
    
    BNmolten= melt(BNsummary,id.vars="m1")
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
      B_sub=BNmolten %>% filter (m1==model)
      B_sub$percent <- B_sub$value / sum(B_sub$value) * 100
      B_percent=rbind(B_percent,B_sub)
    }
    
    BNmolten=B_percent
    cat(paste0(f, " ", sum(BNmolten$value), '\n'))
    
    # # stacked bar plot =================================
    # 
    # pdf(paste0(dir_BN_union_out,"/Switch_models_",f,".pdf"),paper='a4r',7,5)
    # a  <- ggplot(BNmolten, aes(fill=m2, y=value, x=m1 ,label = scales::percent(round(value,digits=2))))+ 
    #   geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
    #   labs(x = "Model in cell 1",y = "Number of triplets switching models", title="Switch in model for shared triplets between cell types", fill='Model in cell 2') +
    #   geom_text(aes(label = paste0(round(value,digits=2))), position = position_stack(vjust = .5),size = 3 ) +
    #   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    #   theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    # print(a)
    # dev.off()
    # 
    # # same with no title
    # 
    # pdf(paste0(dir_BN_union_out,"/Switch_models_notitles_",f,".pdf"),paper='a4r',7,5)
    # 
    # a  <- ggplot(BNmolten, aes(fill=m2, y=value, x=m1 ,label = scales::percent(round(value,digits=2))))+ 
    #   geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
    #   labs(x = "Model in cell 1",y = "Number of triplets switching models", fill='Model in cell 2') +
    #   geom_text(aes(label = paste0(round(value,digits=2))), position = position_stack(vjust = .5),size = 3 ) +
    #   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    #   theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    # print(a)
    # dev.off()
    
    # stacked bar plot percent =================================
    
    pdf(paste0(dir_BN_union_out,"/Switch_models_perc_nt_",f,".pdf"),paper='a4r',7,5)
    a  <- ggplot(BNmolten, aes(fill=m2, y=percent, x=m1 ,label = scales::percent(round(percent,digits=2))))+ 
      geom_bar(positio="stack", stat="identity",width=0.8) + scale_fill_manual(values = color_palette) + 
      labs(x = "Model in cell 1",y = "Percentage of triplets switching models", fill='Model in cell 2') +
      geom_text(aes(label = paste0('n:', round(value,digits=2), "\n",round(percent,digits=2),'%')), position = position_stack(vjust = .5),size = 3 ) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
    print(a)
    dev.off()
    sum(BNmolten$value)
    
}

# analyze the triplet change or not only
summary_change_models=data.frame()

for (file in all.files){
  f=unlist(strsplit(unlist(strsplit(file, "/")[1])[8], ".txt")[1])
  BN1=fread(file, head=FALSE, stringsAsFactors=FALSE)
  colnames(BN1)=c("m1","m2")
  BN1 <- BN1 %>% arrange(m1, m2)
  BN1 <- BN1 %>%  mutate(same_model = if_else(m1 == m2, "equal", "not_equal")) # equal =1
  BNsummary=table(BN1$same_model)
  pt=prop.table(table(BN1$same_model))*100
  pt=as.data.frame.matrix(rbind(pt))
  model=paste(unlist(strsplit(f, "_")[1])[1:4],collapse = '_')
  data_summary=cbind(model,pt)
  summary_change_models=bind_rows(summary_change_models,data_summary )
}

df=melt(summary_change_models) 
colnames(df)=c("model","switch","proportion")

pdf(paste0(df,"/Switch_models_1_0_",f,".pdf"),paper='a4r',7,5)
ggplot(df, aes(fill=switch, y=proportion, x=model)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100) + scale_fill_manual(values = c("#046C9A", "#E7B800" ))    + 
  geom_text(aes(label = paste0(round(proportion,digits=2),'%')), position = position_stack(vjust = .5),size = 3 ) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(text = element_text(size=15),axis.title = element_text(size = 15),axis.text = element_text(size = 10))
dev.off()


# corr plot with switch or not
# this is model swithc or not

library(corrplot)

summary_change_models$equal
histM=c(1,50.96, 47.62,49.46 ,1, 47.60 , 41.76,40.54,1 )
methylM=c(1, 57.76 ,44.51, 52.82, 1, 48.30,48.81 ,43.52,1)
round(summary_change_models$equal,2)
histM=as.data.frame.matrix(histM)

pdf(paste0(dir_BN_union_out,"triplets_not_changing_corrplot_hist.pdf"))
histM=data.frame(N=c(1,52.40, 51.14), M=c(53.62,1, 49.16 ),T=c( 46.66 ,43.01,1) )
colnames(histM) = rownames(histM) =c("NEU","MON","TCL")
histM=as.matrix(histM)
corrplot(histM,is.corr=FALSE, method="circle",cl.lim = c(0, 100),p.mat = histM,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()

pdf(paste0(dir_BN_union_out,"triplets_not_changing_corrplot_methyl.pdf"))
methylM=data.frame(N=c(1, 55.48, 52.32), M=c(60.28,1, 47.95),T=c(52.42, 47.41,1) )
colnames(methylM) = rownames(methylM) =c("NEU","MON","TCL")
methylM=as.matrix(methylM)
corrplot(methylM,is.corr=FALSE, method="circle",cl.lim = c(0, 100),p.mat = methylM,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()







