

#10) R
## plot the distribution of the probability of the most likely model
## this will let you know whether when you call a model causal, reactive or independent, whether the probability is High or not

data <- as.data.frame(data.table::fread("cat BN_cases_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
data1 <- data[, c("causal","reactive","independent")]

data_cntr <- as.data.frame(data.table::fread("cat BN_controls_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
data_contr1 <- data_cntr[, c("causal","reactive","independent")]

overall <- as.data.frame(data.table::fread("cat BN_overall_cov_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
overall1 <- overall[, c("causal","reactive","independent")]

cntr59 <- as.data.frame(data.table::fread("cat BN_controls59_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
contr591 <- cntr59[, c("causal","reactive","independent")]

data1$max <- colnames(data1)[apply(data1,1,which.max)]
data_contr1$max <- colnames(data_contr1)[apply(data_contr1,1,which.max)]
overall1$max <- colnames(overall1)[apply(overall1,1,which.max)]
contr591$max <- colnames(contr591)[apply(contr591,1,which.max)]

barplot(prop.table(table(data1$max)))
barplot(prop.table(table(overall1$max)))
barplot(prop.table(table(data_contr1$max)))
barplot(prop.table(table(contr591$max)))

library(ggplot2)
ggplot(data.frame(data1), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack")

ggplot(data.frame(data_contr1), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack")


ggplot(data.frame(data1), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack",size = 7)+
  theme(legend.title=element_text(size=15), axis.text= element_text(size=30),axis.title.x=element_blank(),plot.title = element_text(size = 20),)+
  ggtitle("Most probable scenario per triplet in cases")

ggplot(data.frame(data_contr1), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack",size = 7)+
  theme(legend.title=element_text(size=15), axis.text= element_text(size=30),axis.title.x=element_blank(),plot.title = element_text(size = 20),)+
  ggtitle("Most probable scenario per triplet in controls")

ggplot(data.frame(overall1), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack",size = 7)+
  theme(legend.title=element_text(size=15), axis.text= element_text(size=30),axis.title.x=element_blank(),plot.title = element_text(size = 20),)+
  ggtitle("Most probable scenario per triplet in combined")






data1$max1 <- pmax(data$causal,data$reactive,data$independent)
data_contr1$max1 <- pmax(data_cntr$causal,data_cntr$reactive,data_cntr$independent)
overall1$max1 <- pmax(overall1$causal,overall1$reactive,overall1$independent)
contr591$max1 <- pmax(contr591$causal,contr591$reactive,contr591$independent)

hist(data_contr1$max1)
hist(overall1$max1)
hist(contr591$max1)

ggplot(data1, aes(x=max,y=max1))+
  geom_boxplot()

ggplot(data_contr1, aes(x=max,y=max1))+
  geom_boxplot()

ggplot(overall1, aes(x=max,y=max1))+
  geom_boxplot()

ggplot(contr591, aes(x=max,y=max1))+
  geom_boxplot()

data2 <- subset (data1, max1 >=0.7)
data2 <- subset (data_contr1, max1 >=0.7)
data2 <- subset (overall1, max1 >=0.8)
data2 <- subset (contr591, max1 >=0.7)

ggplot(data.frame(contr591), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack")

ggplot(data.frame(data2), aes(x=max)) +
  geom_bar()+
  theme_minimal()+
  geom_text(aes(label=..count..), vjust=-0.1,stat='count',position = "stack")


cases <- as.data.frame(data.table::fread("cat /home/users/a/alver/scratch/db_GaP/Analysis/CRD_genes/BN/cases_triplets/BN_cases_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
cases1 <- cases[, c("union_name","causal","reactive","independent")]
cases2 <- cases1
colnames(cases2) = c("union_name","cases_causal","cases_reactive","cases_independent")
colnames(cases2) = c("union_name","causal","reactive","independent") ## to get a nice plot over cases and controls
cases2$cases_max <- colnames(cases2)[apply(cases2,1,which.max)]
cases2$cases_max1 <- pmax(cases2$cases_causal,cases2$cases_reactive,cases2$cases_independent)
cases2$cases_max1 <- pmax(cases2$causal,cases2$reactive,cases2$independent) ## to get a nice plot over cases and controls

cases3 <- subset (cases2, cases_max1 >=0.7)
cases4 <- cases2[, c("union_name","cases_max","cases_max1")]
cases4$status <- 'cases'
colnames(cases4) = c("union_name","max","max1","status")


controls <- as.data.frame(data.table::fread("cat /home/users/a/alver/scratch/db_GaP/Analysis/CRD_genes/BN/controls_triplets/BN_controls_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
controls1 <- controls[, c("union_name","causal","reactive","independent")]
controls2 <- controls1
colnames(controls2) = c("union_name","controls_causal","controls_reactive","controls_independent")
colnames(controls2) = c("union_name","causal","reactive","independent") ## to get a nice plot over cases and controls
controls2$controls_max <- colnames(controls2)[apply(controls2,1,which.max)]
controls2$controls_max1 <- pmax(controls2$controls_causal,controls2$controls_reactive,controls2$controls_independent)
controls2$controls_max1 <- pmax(controls2$causal,controls2$reactive,controls2$independent) ## to get a nice plot over cases and controls

controls3 <- subset (controls2, controls_max1 >=0.7)
controls4 <- controls2[, c("union_name","controls_max","controls_max1")]
controls4$status <- 'controls'
colnames(controls4) = c("union_name","max","max1","status")


controls59 <- as.data.frame(data.table::fread("cat /home/users/a/alver/scratch/db_GaP/Analysis/CRD_genes/BN/controls59_triplets/BN_controls59_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
controls59_1 <- controls59[, c("union_name","causal","reactive","independent")]
controls59_2 <- controls59_1
colnames(controls59_2) = c("union_name","controls_causal","controls_reactive","controls_independent")
colnames(controls59_2) = c("union_name","causal","reactive","independent") ## to get a nice plot over cases and controls
controls59_2$controls_max <- colnames(controls59_2)[apply(controls59_2,1,which.max)]
controls59_2$controls_max1 <- pmax(controls59_2$controls_causal,controls59_2$controls_reactive,controls59_2$controls_independent)
controls59_2$controls_max1 <- pmax(controls59_2$causal,controls59_2$reactive,controls59_2$independent) ## to get a nice plot over cases and controls

controls59_3 <- subset (controls59_2, controls_max1 >=0.7)
controls59_4 <- controls59_2[, c("union_name","controls_max","controls_max1")]
controls59_4$status <- 'controls'
colnames(controls59_4) = c("union_name","max","max1","status")


overall <- as.data.frame(data.table::fread("cat /home/users/a/alver/scratch/db_GaP/Analysis/CRD_genes/BN/overall_cov_triplets/BN_overall_cov_triplet_all_modif2.txt",header=T,stringsAsFactors=F,sep="\t"))
overall1 <- overall[, c("union_name","causal","reactive","independent")]
overall2 <- overall1
colnames(overall2) = c("union_name","overall_causal","overall_reactive","overall_independent")
overall2$overall_max <- colnames(overall2)[apply(overall2,1,which.max)]
overall2$overall_max1 <- pmax(overall2$overall_causal,overall2$overall_reactive,overall2$overall_independent)

overall3 <- subset (overall2, overall_max1 >=0.8)






merge <- merge(cases2,controls2,by.x="union_name",by.y="union_name",all.x=TRUE)
merge1 <- na.omit(merge)

merge_modif <- merge(cases3,controls3,by.x="union_name",by.y="union_name",all.x=TRUE)
merge1_modif <- na.omit(merge_modif)

merge_new <- merge(overall2,merge,by.x="union_name",by.y="union_name",all.x=TRUE)
merge_new_modif <- na.omit(merge_new)

## 1:25561667:G:A;ENSG00000187010.14;1_internal_13807 ## in cases aCRDQTL results, but not cases specific, also in controls conditional top results

write.table(merge1,file="cases_controls_BN_triplet_overlap.txt",row.names=F, quote=F, col.names=T,sep="\t")
write.table(cases2,file="cases_BN_triplet.txt",row.names=F, quote=F, col.names=T,sep="\t")
write.table(controls2,file="controls_BN_triplet.txt",row.names=F, quote=F, col.names=T,sep="\t")


cases_controls <- rbind(cases4, controls4)
cases_controls1 <- subset (cases_controls, max1 >=0.7)

cases_controls_downsample <- rbind(cases4, controls59_4)
cases_controls_ds_1 <- subset (cases_controls_downsample, max1 >=0.7)

causal independent reactive
cases        10          12       12
controls    174         173      134

ggplot(data=cases_controls, aes(x=status,y=max1,fill=max)) +
  geom_bar(stat="identity",position="fill")

ggplot(data=cases_controls_ds_1, aes(x=status,y=max1,fill=max)) +
  geom_bar(stat="identity",position="fill")





## cases_vector <- as.vector(cases['union_name'])
## control_vector <- as.vector(controls['union_name'])
## overall_vector <- as.vector(overall['union_name'])

test <- as.data.frame(setdiff(cases$union_name,controls$union_name))
colnames(test) <- c("a")
test2 <- setdiff(test$a,overall$union_name)


