#### plot the distribution of the probability of the most likely model
## this will let you know whether when you call a model causal, reactive or independent, whether the probability is High or not
# last version
# Clean environment ---------------------------------

rm(list=ls())
gc()

# Packages ---------------------------------

library(qvalue)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)


# Directories and Data ---------------------------------

# directories for mac =================================
dir_BN='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/BN'
dir_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/plots'
# directories for cluster =================================
dir_BN='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/BN'
dir_out='/Users/dianaavalos/Programming/A_CRD_plots/12_bayesian_trios/plots'
# debug =================================
data_type='hist'
cell_type='neut'

# Main ---------------------------------

data_types = list('hist','methyl')
cell_types = list('neut','mono','tcell')

color_palette=c("#046C9A", "#00AFBB", "#E7B800", "#FC4E07","#972D15")


# get data for all categories ---------------------------------

df_total=data.frame()
data_all=data.frame()

for(data_type in data_types){
  for(cell_type in cell_types){
    
    ### loading files
    name_condition=paste0(data_type,'_',cell_type )
    cat(name_condition, '  \n')
    data=fread(paste0(dir_BN,'/',name_condition,'_1000000.txt'), head=FALSE, stringsAsFactors=FALSE)
    colnames(data)=c("name","nbr" ,"var","gene" ,"crd", "dist", "AF","causal","reactive","independent","L1","L2","L3")
    # define best model
    data1=data[,c("causal","reactive","independent")]
    # BN_one_cell$model <- colnames(BN_one_cell[,8:10])[apply(BN_one_cell[,8:10],1,which.max)]
    data$max1 <- pmax(data$causal,data$reactive,data$independent)
    data$model <- colnames(data1)[apply(data1[, c("causal","reactive","independent")],1,which.max)]
    head(data)
    summary(data$max1)
    # get stats for each model
    data_with_prob_sup_05=data %>% filter(max1>0.5)
    # data_with_prob_sup_05=data
    data_all=bind_rows(data_all,data_with_prob_sup_05 )
    t=table(data_with_prob_sup_05$model)
    pt=prop.table(table(data_with_prob_sup_05$model))*100
    # in data.frame
    t=as.data.frame.matrix(rbind(t))
    t$total=t$causal+t$independent+t$reactive
    pt=as.data.frame.matrix(rbind(pt))
    colnames(pt)=c("%causal","%independent","%reactive")
    data_summary=cbind(name_condition,t,pt)
    df_total=bind_rows(df_total,data_summary )
  }
}

data_all$AF_bins <- cut(data_all$AF, breaks=seq(0.0,1.0,0.2))
data_all$dist=abs(data_all$dist)
data_all$dist_bins <- cut(data_all$dist, breaks=c(1000,10000,100000,500000,1000000,10000000000))
head(data_all)
summary(data_all)


## melt
melted_df=melt(df_total[,1:4]) 
colnames(melted_df)=c("name","model","N")
melted_df_perc=melt(df_total[,c(1,6,7,8)]) 
melted_df_tot=cbind(melted_df,melted_df_perc[,3])
colnames(melted_df_tot)=c("name","model","N","percentage")


##### stacked barplot ---------------------------------
pdf(paste0(dir_out,"/models_likelihood_2_nbr.pdf"))

ggplot(melted_df_tot, aes(fill=model, y=percentage, x=name)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100) +
  scale_fill_manual(values=color_palette[c(2,3,4)]) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) +
  labs(x = "type",y ="proportion",  size=15) + ggtitle("Likelihood of the different scenarios") +
  geom_text(aes(label = paste0('n:', N, "\n",round(percentage,digits=2),'%')), position = position_stack(vjust = .5),size = 3 ) 
 
  #geom_text(aes(y=round(percentage,2), label=round(percentage,2)), color="black", vjust=0.5, size=4, position = position_stack(vjust = 0.5))
dev.off()


# histograms of AF and models in x---------------------------------


df_AF=data.frame()
for (model_idx in unique(data_all$model)){
  cat(model_idx)
  data_name= data_all %>% filter(model==model_idx)
  # cut AF in breaks, compute counts and percentage
  t=table(data_name$AF_bins)
  pt=prop.table(table(data_name$AF_bins))*100
  # in data.frame
  t=as.data.frame.matrix(rbind(t))
  pt=as.data.frame.matrix(rbind(pt))
  colnames(pt)=c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]" ,"(0.8,1]")
  data_summary=cbind(model_idx,pt)
  df_AF=bind_rows(df_AF,data_summary )
}

melted_AF_perc=melt(df_AF) 
colnames(melted_AF_perc)=c("model","AF_values","proportion")

pdf(paste0(dir_out,'/models_AF_',".pdf"))
ggplot(melted_AF_perc, aes(fill=AF_values, y=proportion, x=model)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100) + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) +
  labs(x = "model",y ="proportion",  size=15) + ggtitle("Likelihood of the different scenarios") 
dev.off()

pdf(paste0(dir_out,'/models_AF2_',".pdf"))
ggplot(melted_AF_perc, aes(fill=model, y=proportion, x=AF_values)) +  geom_bar(position="fill", stat="identity")+ 
  theme_minimal() + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) 
dev.off()


# histograms of AF and distance in x---------------------------------

head(data_all)
df1=melt(data_all[,c('AF_bins','dist_bins','model')])


df_dist_AF=data.frame()
for (idx in unique(data_all$dist_bins)){
  cat(idx)
  data_name= data_all %>% filter(dist_bins==idx)
  t=table(data_name$AF_bins)
  pt=prop.table(table(data_name$AF_bins))*100
  # in data.frame
  t=as.data.frame.matrix(rbind(t))
  pt=as.data.frame.matrix(rbind(pt))
  data_summary=cbind(idx,pt)
  df_dist_AF=bind_rows(df_dist_AF,data_summary )
}

melted_1=melt(df_dist_AF) 
colnames(melted_1)=c("dist_bins","AF_bins","proportion")

# plot stacked dist in x, AF in fill, percentage in y
ggplot(melted_1, aes(fill=AF_bins, y=proportion, x=dist_bins)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100) + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) +
  labs(x = "distance",y ="proportion",  size=15) + ggtitle(" ") 


ggplot(melted_1, aes(fill=dist_bins, y=proportion, x=AF_bins)) +  geom_bar(position="fill", stat="identity")+ 
  theme_minimal() + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) 

# plot line, problem dist not in the right order
ggplot(melted_1, aes(fill=AF_bins, y=proportion, x=dist_bins)) +  geom_point() + geom_line()  

ggplot(melted_1, aes(x=dist_bins, y=proportion, group=AF_bins, colour=AF_bins)) +
  geom_line() + geom_point() + scale_x_continuous()


####### do the same with dist
head(data_all)


df_dist=data.frame()
for (model_idx in unique(data_all$model)){
  cat(model_idx)
  data_name= data_all %>% filter(model==model_idx)
  t=table(data_name$dist_bins)
  pt=prop.table(table(data_name$dist_bins))*100
  # in data.frame
  t=as.data.frame.matrix(rbind(t))
  pt=as.data.frame.matrix(rbind(pt))
  data_summary=cbind(model_idx,pt)
  df_dist=bind_rows(df_dist,data_summary )
}

melted_dist_perc=melt(df_dist) 
melted_dist_perc=na.omit (melted_dist_perc)
colnames(melted_dist_perc)=c("model","dist","proportion")


pdf(paste0(dir_out,'/models_dist_',".pdf"))
ggplot(melted_dist_perc, aes(fill=dist, y=proportion, x=model)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100) + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) +
  labs(x = "model",y ="proportion",  size=15) + ggtitle("Likelihood of the different scenarios") 
dev.off()

pdf(paste0(dir_out,'/models_dist2_',".pdf"))
ggplot(melted_dist_perc, aes(fill=model, y=proportion, x=dist)) +  geom_bar(position="fill", stat="identity")+ 
  theme_minimal() + 
  scale_fill_viridis(discrete = T) +
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),) 
dev.off()

##### try again with dist
data_all=na.omit (data_all)

df_dist2=data.frame()
for (idx in unique(data_all$dist_bins)){
  data_name= data_all %>% filter(dist_bins==idx)
  t=table(data_name$model)
  pt=prop.table(table(data_name$model))*100
  # in data.frame
  t=as.data.frame.matrix(rbind(t))
  pt=as.data.frame.matrix(rbind(pt))
  data_summary=cbind(idx,pt)
  df_dist2=bind_rows(df_dist2,data_summary )
}

melted_dist2_perc=melt(df_dist2) 
colnames(melted_dist2_perc)=c("dist","model","proportion")


###
pt=prop.table(table(melted_dist_perc$model))*100

ggplot(melted_dist_perc, aes(fill=model, y=proportion, x=dist)) +  geom_bar(position="fill", stat="identity") 

  
  

ggplot(data=df_dist2, aes(x=idx, y=proportion, group=model)) +
  geom_line()+
  geom_point()



#### old

histinfo<-hist(data1$AF)
histinfo$counts

# cut AF in breaks
data1$AF_bins <- cut(data1$AF, breaks=seq(0.0,1.0,0.05))
#melt
melted_AF=melt(data1[,c(1,4,5)]) 

#stacked barplot
ggplot(melted_AF, aes(fill=model, y=percentage, x=AF_bins)) + 
  geom_bar(position="stack", stat="identity")+ theme_minimal() + ylim(0, 100)
ggplot(data1, aes(x=dist, fill=model)) +
geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')

ggplot(data1, aes(x=AF, fill=model)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity')
# melted_subset=melt(data_subset) 



###

### old

geom_text(aes(y=percentage, label=paste0(round(percentage,2)),'%\n n=',N), color="black", vjust=0.5, position = position_stack(vjust = 0.5))


ggplot(melted_df, aes(fill=perc_variable, y=percentage, x=name_condition)) + 
  geom_bar(stat = "stack", fill=color_palette[c(2,3,4)]) + theme_minimal() + ylim(0, 100)  + 
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),)+
  labs(x = "name_condition",y ="variable",vjust = -.5,  size=15) + ggtitle("Likelihood of the different scenarios") +
  geom_text(aes(label=paste0(round(data_summary$perc,2))), vjust=1.6, color="black")

#

ggplot(melted_df, aes(fill=perc_variable, y=percentage, x=name_condition)) + 
  geom_bar(stat = "stack", fill=color_palette[c(2,3,4)]) + theme_minimal() + ylim(0, 100)  + 
  theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),)+
  labs(x = "name_condition",y ="variable",vjust = -.5,  size=15) + ggtitle("Likelihood of the different scenarios") +
  geom_text(aes(label=paste0(round(data_summary$perc,2))), vjust=1.6, color="black")

####
colnames(data_with_prob_sup_05)=c("causal","reactive","independent", "max1", "model")
data_with_prob_sup_05 %>% count(model) %>% mutate(perc = n / nrow(data1)*100) -> data_summary


#############################################################################################
#
# plots per condtion
#
#############################################################################################


df_total=data.frame()
for(data_type in data_types){
  for(cell_type in cell_types){
    
    ### loading files
    name_condition=paste0(data_type,'_',cell_type )
    cat(name_condition, '  \n')
    data=fread(paste0(dir_BN,'/','BN_',name_condition,'.txt'), head=FALSE, stringsAsFactors=FALSE)
    colnames(data)=c("name","nbr" ,"var","gene" ,"crd", "causal","reactive","independent","L1","L2","L3")
    data1 <- data[, c("causal","reactive","independent")]
    data1$max1 <- pmax(data$causal,data$reactive,data$independent)
    data1$model <- colnames(data1)[apply(data1,1,which.max)]
    head(data1)
    data1 %>% count(model) %>% mutate(perc = n / nrow(data1)*100) -> data_summary
    data_summary$condition=name_condition
    df_total=bind_rows(df_total,data_summary )
    #############################################################################################
    #
    # PLOT models
    #
    #############################################################################################
    
    pdf(paste0(dir_out,'/models_likelihood_',name_condition,".pdf"))
    ggplot(data_summary, aes(x = model, y = perc)) + 
      geom_bar(stat = "identity", fill=color_palette[c(2,3,4)]) + theme_minimal() + ylim(0, 100)  + 
      theme(legend.title=element_text(size=15), axis.text= element_text(size=15),axis.title= element_text(size=15),plot.title = element_text(size = 15),)+
      labs(x = "Models",y = "Percentage",vjust = -.5,  size=15) + ggtitle("Likelihood of the different scenarios") +
      geom_text(aes(label=paste0(round(data_summary$perc,2))), vjust=1.6, color="black")
    dev.off()
    
    #############################################################################################
    #
    # PLOT distrib
    #
    #############################################################################################
    
    pdf(paste0(dir_out,'/distrib_best_proba_',name_condition,".pdf"))
    hist(data1$max1,
         main="Probability of the most likely model",
         xlab="bayesian probability",
         ylab="frequency")
    dev.off()
    
    ### stacked barplot
    # ggplot(data_summary, aes(fill=perc, y=perc, x=model)) + 
    #   geom_bar(position="fill", stat="identity")
  }
}
