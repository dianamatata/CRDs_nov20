import glob
import pandas as pd
import csv


def write_to_file(dataframe, outputfile):
    """
    write dataframe in cvs file
    :param dataframe:
    :param outputfile: cvs filename
    """
    w = csv.writer(open(outputfile, "w"))
    for key, val in dataframe.items():
        w.writerow([key, val])


path1='/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5/significants'
path2='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:TRHs'

# with community file from igraph community algo

for cell in ('neut', 'mono', 'tcell'):
    for data_type in ('methyl','hist'):

        CRD_genes_list = pd.read_csv("%s/FDR_0.05_%s_%s_mean_mapping_CRD_gene_ALL.significant.txt" % (path1, data_type, cell), header=None, sep=' ')
        CRD_cluster_file = pd.read_csv('%s/TRH_%s_%s_mean_trans.significant_0.01.txt' % (path2,data_type, cell), header=0, sep='\t')

        CRD_genes_list.columns = ["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end", "phenotype_ID_strand",
          "nb_variants", "distance", "CRD_ID", "CRD_ID_chr", "CRD_ID_start", "CRD_ID_end", "rank",
          "fwd_pval", "fwd_r_squared", "fwd_slope", "fwd_best_hit", "fwd_sig", "bwd_pval", "bwd_r_squared", "bwd_slope",
          "bwd_best_hit", "bwd_sig"]

        # loop on clusters
        for cluster in range(1, max(CRD_cluster_file['communities.membership'])+1):
            CRDS_in_TRH=CRD_cluster_file['communities.names'][CRD_cluster_file['communities.membership'] == cluster]

            genes = pd.DataFrame(columns=["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end",
                                          "phenotype_ID_strand", "CRD_ID"])
            for CRD in CRDS_in_TRH:
                #print(CRD)  # CRD = '11_internal_10924'
                subdf = CRD_genes_list.loc[CRD_genes_list['CRD_ID'] == CRD].iloc[:, [0, 1, 2, 3, 4, 7]]
                genes = genes.append(subdf)

            print(len(genes))

            if (len(genes) > 15):

                    outname = "%s/trans_hubs_genes/genes_%s_%s_cluster%s.csv" % (path2, data_type, cell, cluster)
                    genes.to_csv(outname, index=False, sep='\t')
                    genes.phenotype_ID = genes.phenotype_ID.apply(lambda x: x.split(".")[0])
                    genes.phenotype_ID.to_csv("%s/trans_hubs_genes/genesonly_%s_%s_cluster%s.csv" % (path2, data_type, cell, cluster), index=False, sep='\t', header=False)

                    with open("%s/trans_hubs_genes/genesonly_%s_%s_cluster%s.csv" % (path2, data_type, cell, cluster),'w') as outfile:
                        genes.phenotype_ID.to_string(outfile) #pb it keeps the index
