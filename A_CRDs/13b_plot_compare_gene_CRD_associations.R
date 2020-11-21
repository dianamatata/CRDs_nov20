# Plot 3.5 and 3.4 from paper
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS
path_out="/Users/dianaavalos/Programming/A_CRD_plots/figs_7_Rfile"
name_condition="hist"

# 3.5 not anymore
library(ggplot2)
library(corrplot)
library(reshape2)
library(data.table)
library(ggmosaic)

get_common_fraction <- function(set,refset,mindist,maxdist){
    if(mindist == 0 & maxdist == 0){
        set.filtered = set[which(set$distance==0),c(1,8)]
        refset.filtered = refset[which(refset$distance==0),c(1,8)]
    } else if(maxdist > 0) {
        set.filtered = set[which(set$distance>mindist & set$distance<=maxdist),c(1,8)]
        refset.filtered = refset[which(refset$distance>mindist & refset$distance<=maxdist),c(1,8)]
    } else {
        set.filtered = set[,c(1,8)]
        refset.filtered = refset[,c(1,8)]
    }
    myfraction = sum(duplicated(rbind(set.filtered,refset.filtered)))/nrow(refset.filtered)
    myfraction
}

plot_fractions <-function(fractions,discovered,replicated) {
    filename = paste0("CRD_gene_association_discovered_in_",discovered,"_replicated_in_",replicated,".pdf")
    pdf(filename,paper="a4r")
    dat_assos = data.frame(Distance = c("0","1-\n1e03","1e03-\n1e04","1e04-\n1e05","1e05-\n1e06"),fraction = fractions)
    dat_assos$Distance <- factor(dat_assos$Distance,levels = c("0","1-\n1e03","1e03-\n1e04","1e04-\n1e05","1e05-\n1e06"))

    p <- ggplot(dat_assos, aes(x = Distance, y = fraction))+ ggtitle(paste0("Gene-CRD associations in ",discovered)) +
      geom_bar(stat = "identity",fill="steelblue") + ylim(0,0.5) +
      labs(x = "Distance between genes and CRDs (bp)",y = paste0("Fraction of significant associations in ",replicated)) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
    print(p)
    dev.off()
}

### main

directory='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/CRD_genes_scripts/analysis_files'

map70_vs_70 = read.table(paste0(directory,'/70_vs_70_mapping_gene_CRD_mean_ALL.txt'))
map70_vs_72 = read.table(paste0(directory,'/70_vs_72_mapping_gene_CRD_mean_ALL.txt'))
map70_vs_73 = read.table(paste0(directory,'/70_vs_73_mapping_gene_CRD_mean_ALL.txt'))
map72_vs_70 = read.table(paste0(directory,'/72_vs_70_mapping_gene_CRD_mean_ALL.txt'))
map72_vs_72 = read.table(paste0(directory,'/72_vs_72_mapping_gene_CRD_mean_ALL.txt'))
map72_vs_73 = read.table(paste0(directory,'/72_vs_73_mapping_gene_CRD_mean_ALL.txt'))
map73_vs_70 = read.table(paste0(directory,'/73_vs_70_mapping_gene_CRD_mean_ALL.txt'))
map73_vs_72 = read.table(paste0(directory,'/73_vs_72_mapping_gene_CRD_mean_ALL.txt'))
map73_vs_73 = read.table(paste0(directory,'/73_vs_73_mapping_gene_CRD_mean_ALL.txt'))
colnames(map70_vs_70) = colnames(map70_vs_72) = colnames(map73_vs_73) =colnames(map70_vs_73) = colnames(map72_vs_70) = colnames(map72_vs_72) = colnames(map72_vs_73) = colnames(map73_vs_70) = colnames(map73_vs_72) =  c("gene","gene_chr","gene_start","gene_end","gene_strand","dummy","distance","CRD",
"CRD_chr","CRD_start","CRD_end","pval","slope","rank")


# set = map72_vs_70
# refset = map70_vs_70
# mindist = -1
# maxdist = 1
fraction_neut_replicated_in_mono = get_common_fraction(map72_vs_70,map70_vs_70,-1,-1)
fraction_neut_replicated_in_mono_bin1 = get_common_fraction(map72_vs_70,map70_vs_70,0,0)
fraction_neut_replicated_in_mono_bin2 = get_common_fraction(map72_vs_70,map70_vs_70,1,1e03)
fraction_neut_replicated_in_mono_bin3 = get_common_fraction(map72_vs_70,map70_vs_70,1e03,1e04)
fraction_neut_replicated_in_mono_bin4 = get_common_fraction(map72_vs_70,map70_vs_70,1e04,1e05)
fraction_neut_replicated_in_mono_bin5 = get_common_fraction(map72_vs_70,map70_vs_70,1e05,1e06)
fraction_neut_replicated_in_mono_all = c(fraction_neut_replicated_in_mono_bin1,fraction_neut_replicated_in_mono_bin2,fraction_neut_replicated_in_mono_bin3,
    fraction_neut_replicated_in_mono_bin4,fraction_neut_replicated_in_mono_bin5)

# bar plot of the distance between genes and CRDs associations
plot_fractions(fraction_neut_replicated_in_mono_all,"neutrophils","monocytes")
# fractions = fraction_neut_replicated_in_mono_all
# discovered ="neutrophils"
# replicated ="monocytes"

fraction_neut_replicated_in_tcel = get_common_fraction(map73_vs_70,map70_vs_70,-1,-1)
fraction_neut_replicated_in_tcel_all = c(get_common_fraction(map73_vs_70,map70_vs_70,0,0),
                                         get_common_fraction(map73_vs_70,map70_vs_70,1,1e03),
                                         get_common_fraction(map73_vs_70,map70_vs_70,1e03,1e04),
                                         get_common_fraction(map73_vs_70,map70_vs_70,1e04,1e05),
                                         get_common_fraction(map73_vs_70,map70_vs_70,1e05,1e06))

fractions=fraction_neut_replicated_in_tcel_all
discovered="neutrophils"
replicated="T cells"
plot_fractions(fraction_neut_replicated_in_tcel_all,"neutrophils","T cells")

fraction_mono_replicated_in_neut = get_common_fraction(map70_vs_72,map72_vs_72,-1,-1)
fraction_mono_replicated_in_neut_all = c(get_common_fraction(map70_vs_72,map72_vs_72,0,0),
                                         get_common_fraction(map70_vs_72,map72_vs_72,1,1e03),
                                         get_common_fraction(map70_vs_72,map72_vs_72,1e03,1e04),
                                         get_common_fraction(map70_vs_72,map72_vs_72,1e04,1e05),
                                         get_common_fraction(map70_vs_72,map72_vs_72,1e05,1e06))

plot_fractions(fraction_mono_replicated_in_neut_all,"monocytes","neutrophils")

fraction_mono_replicated_in_tcel = get_common_fraction(map73_vs_72,map72_vs_72,-1,-1)
fraction_mono_replicated_in_tcel_bin1 = get_common_fraction(map73_vs_72,map72_vs_72,0,0)
fraction_mono_replicated_in_tcel_bin2 = get_common_fraction(map73_vs_72,map72_vs_72,1,1e03)
fraction_mono_replicated_in_tcel_bin3 = get_common_fraction(map73_vs_72,map72_vs_72,1e03,1e04)
fraction_mono_replicated_in_tcel_bin4 = get_common_fraction(map73_vs_72,map72_vs_72,1e04,1e05)
fraction_mono_replicated_in_tcel_bin5 = get_common_fraction(map73_vs_72,map72_vs_72,1e05,1e06)
fraction_mono_replicated_in_tcel_all = c(fraction_mono_replicated_in_tcel_bin1,fraction_mono_replicated_in_tcel_bin2,fraction_mono_replicated_in_tcel_bin3,
    fraction_mono_replicated_in_tcel_bin4,fraction_mono_replicated_in_tcel_bin5)

plot_fractions(fraction_mono_replicated_in_tcel_all,"monocytes","T cells")

fraction_tcel_replicated_in_neut = get_common_fraction(map70_vs_73,map73_vs_73,-1,-1)
fraction_tcel_replicated_in_neut_bin1 = get_common_fraction(map70_vs_73,map73_vs_73,0,0)
fraction_tcel_replicated_in_neut_bin2 = get_common_fraction(map70_vs_73,map73_vs_73,1,1e03)
fraction_tcel_replicated_in_neut_bin3 = get_common_fraction(map70_vs_73,map73_vs_73,1e03,1e04)
fraction_tcel_replicated_in_neut_bin4 = get_common_fraction(map70_vs_73,map73_vs_73,1e04,1e05)
fraction_tcel_replicated_in_neut_bin5 = get_common_fraction(map70_vs_73,map73_vs_73,1e05,1e06)
fraction_tcel_replicated_in_neut_all = c(fraction_tcel_replicated_in_neut_bin1,fraction_tcel_replicated_in_neut_bin2,fraction_tcel_replicated_in_neut_bin3,
    fraction_tcel_replicated_in_neut_bin4,fraction_tcel_replicated_in_neut_bin5)

plot_fractions(fraction_tcel_replicated_in_neut_all,"T cells","neutrophils")

fraction_tcel_replicated_in_mono = get_common_fraction(map72_vs_73,map73_vs_73,-1,-1)
fraction_tcel_replicated_in_mono_bin1 = get_common_fraction(map72_vs_73,map73_vs_73,0,0)
fraction_tcel_replicated_in_mono_bin2 = get_common_fraction(map72_vs_73,map73_vs_73,1,1e03)
fraction_tcel_replicated_in_mono_bin3 = get_common_fraction(map72_vs_73,map73_vs_73,1e03,1e04)
fraction_tcel_replicated_in_mono_bin4 = get_common_fraction(map72_vs_73,map73_vs_73,1e04,1e05)
fraction_tcel_replicated_in_mono_bin5 = get_common_fraction(map72_vs_73,map73_vs_73,1e05,1e06)
fraction_tcel_replicated_in_mono_all = c(fraction_tcel_replicated_in_mono_bin1,fraction_tcel_replicated_in_mono_bin2,fraction_tcel_replicated_in_mono_bin3,
    fraction_tcel_replicated_in_mono_bin4,fraction_tcel_replicated_in_mono_bin5)

plot_fractions(fraction_tcel_replicated_in_mono_all,"T cells","monocytes")


#### plot 3.4 CRD-gene tissue sharing


path_out="/Users/dianaavalos/Programming/A_CRD_plots/figs_7_Rfile"
pdf(paste0(path_out,'/',name_condition,"_3.4_Pairwise_comparisons_between_cell_types.pdf"))

M = matrix(c(1,fraction_neut_replicated_in_mono,fraction_neut_replicated_in_tcel,fraction_mono_replicated_in_neut,1,fraction_mono_replicated_in_tcel,
    fraction_tcel_replicated_in_neut,fraction_tcel_replicated_in_mono,1),ncol=3,byrow=T)
colnames(M) = c("Neutrophils","Monocytes","T cells")
rownames(M) = c("Neutrophils","Monocytes","T cells")
#corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()

## Fig 3.5 , check in plots paper
# how did i find this data? >fraction_of_genes_shared_btw_tissues.py
