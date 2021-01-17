### Goal: plot fig 3B: Overlap between genes in CRD-gene associations across cell types.

# libs and functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def calculate_gene_sharing_percentage(a, b, c):
    """
    :param a = set(dict_df['map70_vs_70'].gene)
    :param b = set(dict_df['map72_vs_70'].gene)
    :param c = set(dict_df['map73_vs_70'].gene)
    :return: fractions of genes shared by 1 2 or 3 cell types
    """
    shared3 = len(a.intersection(b).intersection(c)) / len(a) * 100
    shared2 = (len(a.intersection(b)) + len(a.intersection(c)) - 2 * len(a.intersection(b).intersection(c))) / len(
        a) * 100
    shared1 = (100 - shared3 - shared2)  # a.difference(b).difference(d)
    return round(shared1, 2), round(shared2, 2), round(shared3, 2)


def calculate_gene_sharing_number(a, b, c):
    """
    :param a = set(dict_df['map70_vs_70'].gene)
    :param b = set(dict_df['map72_vs_70'].gene)
    :param c = set(dict_df['map73_vs_70'].gene)
    :return: fractions of genes shared by 1 2 or 3 cell types
    """
    shared3 = len(a.intersection(b).intersection(c))
    shared2 = (len(a.intersection(b)) + len(a.intersection(c)) - 2 * len(a.intersection(b).intersection(c)))
    shared1 = len(a.difference(b).difference(c))
    return round(shared1, 2), round(shared2, 2), round(shared3, 2)


def survey(results, category_names, category_colors=None):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    if category_colors is None:
        category_colors = plt.get_cmap('Dark2')(
            np.linspace(0.15, 0.85, data.shape[1]))
    # category_colors=['teal', 'tomato', 'silver']

    fig, ax = plt.subplots(figsize=(9.2, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        # r, g, b, _ = color
        # text_color = 'black' if colname=='3 cell types' else 'white'
        text_color = 'white'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax.text(x, y, str(int(c)), ha='center', va='center',
                    color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return fig, ax


# create dict of df
column_names = ["gene", "gene_chr", "gene_start", "gene_end", "gene_strand", "dummy", "distance", "CRD",
                "CRD_chr", "CRD_start", "CRD_end", "pval", "slope", "rank"]


dict_df = {}
for i in [70, 72, 73]:
    for j in [70, 72, 73]:
        name = 'map' + str(i) + '_vs_' + str(j)
        # hist_mono_loom_mapping_CRD_gene_ALL.txt.gz
        # UPDATE PATH and load files
        path = '/Users/dianaavalos/Programming/A_CRD_plots/CRD_genes_5//merged/hist' + str(i) + '_vs_' + str(j) \
               + '_mapping_CRD_gene_ALL.txt.gz'
        dict_df[name] = pd.read_csv(path, sep=' ', header=None)
        dict_df[name].columns = column_names
        exec(name + " = pd.read_csv(path, sep=' ', header=None) ")

# for neutrophils CRDs (2nd 70), genes shared with 1 cell type, same for mono and tcells then
neut1, neut2, neut3 = calculate_gene_sharing_percentage(
    set(dict_df['map70_vs_70'].gene), set(dict_df['map72_vs_70'].gene), set(dict_df['map73_vs_70'].gene))
mono1, mono2, mono3 = calculate_gene_sharing_percentage(
    set(dict_df['map70_vs_72'].gene), set(dict_df['map72_vs_72'].gene), set(dict_df['map73_vs_72'].gene))
tcell1, tcell2, tcell3 = calculate_gene_sharing_percentage(
    set(dict_df['map70_vs_73'].gene), set(dict_df['map72_vs_73'].gene), set(dict_df['map73_vs_73'].gene))

# for neutrophils CRDs (2nd 70), genes shared with 1 cell type, same for mono and tcells then
neu1, neu2, neu3 = calculate_gene_sharing_number(
    set(dict_df['map70_vs_70'].gene), set(dict_df['map72_vs_70'].gene), set(dict_df['map73_vs_70'].gene))
mon1, mon2, mon3 = calculate_gene_sharing_number(
    set(dict_df['map70_vs_72'].gene), set(dict_df['map72_vs_72'].gene), set(dict_df['map73_vs_72'].gene))
tcel1, tcel2, tcel3 = calculate_gene_sharing_number(
    set(dict_df['map70_vs_73'].gene), set(dict_df['map72_vs_73'].gene), set(dict_df['map73_vs_73'].gene))

# plot

category_names = ['1 cell type', '2 cell types', '3 cell types']
results = {
    'Neutrophils': [neut1, neut2, neut3],
    'Monocytes': [mono1, mono2, mono3],
    'T-cells': [tcell1, tcell2, tcell3],
}

results_not_normalized = {
    'Neutrophils': [neu1, neu2, neu3],
    'Monocytes': [mon1, mon2, mon3],
    'T-cells': [tcel1, tcel2, tcel3],
}

# colors  c("#00AFBB", "#E7B800", "#FC4E07") c(NEU="#000066",MON="#663399",TCL="#339999") GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
# c("#000066", "#663399", "#339999"), greyscale= ['dimgrey','darkgrey', 'gainsboro']

fig, ax = survey(results, category_names, ["#F1BB7B", "#FD6467", "#5B1A18"])
# fig, ax = survey(results, category_names, ['navy','teal', 'turquoise'])
plt.box(False)
ax.set_xlabel("Fraction of overlap in %")
plt.xlabel("2013", fontsize=12)
plt.ylabel('Cell Types')
plt.show()
fig.savefig('overlap3b.svg', format='svg', dpi=1200)

overlap_trans_egenes_result3 = {
    'Neutrophils': [97.101, 1.449, 1.449],
    'Monocytes': [95.1219, 2.439, 2.439],
    'T-cells': [75, 0, 25],
}

overlap_trans_egenes_result3_not_norm = {
    'Neutrophils': [67, 1, 1],
    'Monocytes': [39, 1, 1],
    'T-cells': [3, 0, 1],
}

fig, ax = survey(overlap_trans_egenes_result3_not_norm, category_names, ["#F1BB7B", "#FD6467", "#5B1A18"])
plt.box(False)
ax.set_xlabel("Fraction of overlap in %")
plt.xlabel("2013", fontsize=12)
plt.ylabel('Cell Types')
plt.show()
fig.savefig('overlap_trans_egenes_result3notnorm.svg', format='svg', dpi=1200)

fig, ax = survey(overlap_trans_egenes_result3, category_names, ["#F1BB7B", "#FD6467", "#5B1A18"])
plt.box(False)
ax.set_xlabel("Fraction of overlap in %")
plt.xlabel("2013", fontsize=12)
plt.ylabel('Cell Types')
plt.show()
fig.savefig('overlap_trans_egenes_result3.svg', format='svg', dpi=1200)

# https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/horizontal_barchart_distribution.html#sphx-glr-gallery-lines-bars-and-markers-horizontal-barchart-distribution-py
# defining an array of color
#  not finished yet
bin0 = len(dict_df['map72_vs_70'].distance >= 1000)

plt.hist(abs(dict_df['map70_vs_70'].distance), bins=[0, 1, 1000, 10e4, 10e5, 10e6])
hist, bins, _ = plt.hist(abs(dict_df['map70_vs_70'].distance), bins=8)
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
plt.show()
plt.hist(abs(dict_df['map70_vs_70'].distance), bins=logbins)
plt.xscale('log')
plt.show()
