import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# file list
fileList = {
    "i_70_72": "CRD_scripts/files_for_python/70_72.ALLchr.module.merged.intersect.bed",
    "i_70_73": "CRD_scripts/files_for_python/70_73.ALLchr.module.merged.intersect.bed",
    "i_72_73": "CRD_scripts/files_for_python/72_73.ALLchr.module.merged.intersect.bed",
    "i_70_72_73": "CRD_scripts/files_for_python/70_72_73.ALLchr.module.merged.intersect.bed",

    "m_70": "CRD_scripts/files_for_python/EGAD00001002670.ALLchr.module.merged.bed",
    "m_72": "CRD_scripts/files_for_python/EGAD00001002672.ALLchr.module.merged.bed",
    "m_73": "CRD_scripts/files_for_python/EGAD00001002673.ALLchr.module.merged.bed"
}
# get nbr of bp per file and append in total_ bp
total_bp = pd.DataFrame(columns=['file', 'total_bp'])
T = {}
for file, filename in fileList.items():
    print(file)
    data = pd.read_csv(filename, sep="\t", header=None)
    data["bp"] = data[2] - data[1]
    total_bp_file = data["bp"].sum(axis=0)
    total_bp = total_bp.append({'file': file, 'total_bp': total_bp_file}, ignore_index=True)
    T[file] = total_bp_file

# matrix of overlap

overlap_array = ([[[1, T['i_70_72'] / T['m_70'], T['i_70_73'] / T['m_70']],
                   [T['i_70_72'] / T['m_72'], 1, T['i_72_73'] / T['m_72']],
                   [T['i_70_73'] / T['m_73'], T['i_72_73'] / T['m_73'], 1]]])

overlap_array = np.around(overlap_array, decimals=2)

# sns.heatmap(overlap_array)
# add empty diag and names and axes and title and round and size of rounds)

sns.set()
uniform_data = np.random.rand(3,3)
#ax = sns.heatmap(uniform_data, vmin=0, vmax=1,xticklabels=['NEUT','MONO','TCELL'], yticklabels=['NEUT','MONO','T-CELL'] )

ax = sns.heatmap(uniform_data, rasterized=False, annot=True, cmap="Blues", linewidths=4)
ax.set(xlabel='Query', ylabel='Reference', title='CRD Tissue Sharing')
plt.show()