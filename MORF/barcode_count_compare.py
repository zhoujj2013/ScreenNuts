#%%
import os
import argparse
import csv
from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
# 添加参数
parser.add_argument('-i', type=str, required=True,help='sample file',)
parser.add_argument('-b',type=str,help='barcode list,csv file,example /mnt/dfc_data2/project/linyusen/project/81_MORF/barcode_info.csv',required=True)
parser.add_argument('-o', type=str, required=True,help='Output directory to save plots and results')

args = parser.parse_args()

input = args.i
output_dir = args.o
barcode_file = args.b
os.makedirs(output_dir,exist_ok=True)
# input = '/mnt/dfc_data2/project/linyusen/project/81_MORF/samples.lst'
# output_dir = '/mnt/dfc_data2/project/linyusen/project/81_MORF/data/compare'



print('Precessing ... ')
barcode_dict = {}
with open(barcode_file, 'r') as infile:
    reader = csv.reader(infile)
    reader.__next__()
    for rows in reader:
        barcode_dict[rows[3]] = 0

barcode_tf_dict = {}
with open(barcode_file, 'r') as infile:
    reader = csv.reader(infile)
    reader.__next__()
    for rows in reader:
        barcode_tf_dict[rows[3]] = {}
        barcode_tf_dict[rows[3]]['Name'] = rows[0]
        barcode_tf_dict[rows[3]]['RefSeqGeneName'] = rows[1]
        barcode_tf_dict[rows[3]]['RefSeqandGencodeID'] = rows[2]
        barcode_tf_dict[rows[3]]['BarcodeSequence'] = rows[3]


f = open(input)
input_dir1_list = []
input_dir2_list = []
for i in f.readlines():
    i = i.strip('\n').split(' ')
    if i[0] == '1':
        input_dir1_list.append(i[-1])
        sample1 = i[1]
    if i[0] == '2':
        input_dir2_list.append(i[-1])
        sample2 = i[1]


TF_Count_file1_list = []
TF_Count_file2_list = []

for path in input_dir1_list:
    TF_Count_file1_list.append(os.path.join(path, 'TF_Count.txt'))
for path in input_dir2_list:
    TF_Count_file2_list.append(os.path.join(path, 'TF_Count.txt'))
#%%

# Function to read data from a CSV file and return a dictionary {barcode: cpm}
def read_data(file_path):
    f = open(file_path)
    r = csv.reader(f,delimiter='\t')
    r.__next__()  # Skip the first header line
    r.__next__()  # Skip the second header line
    data = {}
    for i in r:
        barcode,tf,gene,id,count, cpm = i
        cpm = float(cpm)          # Convert cpm to float
        data[barcode] = cpm       # Store barcode → cpm mapping
    return data

# Lists to store data dictionaries from multiple files
TF_Count1_list = []
TF_Count2_list = []

# Read data from all files in TF_Count_file1_list and TF_Count_file2_list
for file in TF_Count_file1_list:
    TF_Count1_list.append(read_data(file))
for file in TF_Count_file2_list:
    TF_Count2_list.append(read_data(file))

# Merge data from multiple files into one dictionary, grouped by barcode
# Each barcode maps to a list of CPMs (from different replicates/files)
TF_Count1 = {}
for temp in TF_Count1_list:
    for i in temp:
        if i not in TF_Count1:
            TF_Count1[i] = []
        TF_Count1[i].append(temp[i])

TF_Count2 = {}
for temp in TF_Count2_list:
    for i in temp:
        if i not in TF_Count2:
            TF_Count2[i] = []
        TF_Count2[i].append(temp[i])

# For each barcode, take the average CPM across all files
for i in TF_Count1:
    TF_Count1[i] = np.mean(TF_Count1[i])
for i in TF_Count2:
    TF_Count2[i] = np.mean(TF_Count2[i])
#%%


# 将 TF 的顺序固定好（交集）
common_keys = set(TF_Count1.keys()) & set(TF_Count2.keys())
x = [TF_Count1[k] for k in common_keys]
y = [TF_Count2[k] for k in common_keys]

# 计算皮尔森相关系数
r, p_value = pearsonr(x, y)

# 画散点图
plt.figure(figsize=(7, 7))
plt.scatter(x, y, alpha=0.6, color='steelblue', edgecolor='k')
plt.xlabel(f"TF CPM - {sample1}")
plt.ylabel(f"TF CPM - {sample2}")
plt.title(f"Pearson r = {r:.3f}, p = {p_value:.1e}")

# 添加对角线参考线（可选）
max_val = max(max(x), max(y))
plt.plot([0, max_val], [0, max_val], 'r--', lw=1)

plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'barcode_scatter_plot.png'))

#%%


# 伪计数，避免除以 0
log2fc = {
    k: np.log2((TF_Count1[k] + 0.1) / (TF_Count2[k] + 0.1))
    for k in set(TF_Count1) & set(TF_Count2)
}

# 排序（从大到小）
log2fc_sorted = dict(sorted(log2fc.items(), key=lambda x: x[1], reverse=True))

# 转换成 DataFrame（只有一列）
df = pd.DataFrame.from_dict(log2fc_sorted, orient='index', columns=['log2FC'])

vmax = -df['log2FC'].min()
vmin = df['log2FC'].min()

# 热图
fig, ax = plt.subplots(figsize=(4,8))  # 控制宽度和高度

fig.subplots_adjust(left=0.4, right=0.6)
sns.heatmap(
    df,
    cmap="coolwarm",
    cbar=True,
    yticklabels=False,
    ax=ax,
    vmax=vmax,
    vmin=vmin
)

# 美化
ax.set_title('Heatmap', fontsize=14)
ax.set_xticklabels([f'log2FC({sample1}/{sample2})'], rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'barcode_heatmap.png'))



#%%


# 假设有两个条件的计数字典
# TF_CountA, TF_CountB

# 统一 TF 集合
common_tfs = set(TF_Count1) & set(TF_Count2)

# 构建 DataFrame
df = pd.DataFrame({
    sample1: [TF_Count1[tf] for tf in common_tfs],
    sample2: [TF_Count2[tf] for tf in common_tfs],
}, index=list(common_tfs))

# 排名（rank 越小表示丰度越高）
df[f'{sample1}_rank'] = df[sample1].rank(ascending=False)
df[f'{sample2}_rank'] = df[sample2].rank(ascending=False)

# Rank-Rank plot
plt.figure(figsize=(6, 6))
sns.scatterplot(x=sample1, y=sample2, data=df, s=10)
plt.title('Rank-Rank Plot')
plt.xlabel(f'Condition {sample1} Rank')
plt.ylabel(f'Condition {sample2} Rank')
plt.plot([0, len(df)], [0, len(df)], color='gray', linestyle='--')  # y=x 参考线
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'Rank-Rank_Plot.png'))

#%%
# 构建字典：TF -> (rank1, rank2)
rank_dict = {
    tf: (df.loc[tf, f'{sample1}_rank'], df.loc[tf, f'{sample2}_rank'])
    for tf in df.index
}


#%%


# Dropout 判断逻辑
dropout_tfs_1 = [k for k in TF_Count1 if TF_Count1[k] == 0 and TF_Count2.get(k, 0) > 0]
dropout_tfs_2 = [k for k in TF_Count2 if TF_Count2[k] == 0 and TF_Count1.get(k, 0) > 0]

dropout_counts = [len(dropout_tfs_1), len(dropout_tfs_2)]
labels = [f'Dropout in {sample1}', f'Dropout in {sample2}']

# 绘图
plt.figure(figsize=(5, 4))
plt.bar(labels, dropout_counts, color=['salmon', 'skyblue'], width=0.5)
plt.title('Dropout TF Counts Between Samples')
plt.ylabel('Number of TFs')
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'Dropout_TF_Count.png'))


f = open(os.path.join(output_dir,'compare.info.txt'),'w')
w = csv.writer(f,delimiter='\t')
w.writerow(['Barcode','Name','RefSeqGeneName','RefSeqandGencodeID',f'{sample1}_count',f'{sample1}_count',f'log2FC({sample1}/{sample2})',f'{sample1}_rank',f'{sample1}_rank'])
for barcode in barcode_tf_dict:
    if barcode not in common_keys:
        continue
    w.writerow([barcode,barcode_tf_dict[barcode]['Name'],barcode_tf_dict[barcode]['RefSeqGeneName'],barcode_tf_dict[barcode]['RefSeqandGencodeID'],TF_Count1[barcode],TF_Count2[barcode],log2fc[barcode],rank_dict[barcode][0],rank_dict[barcode][1]])
f.close()
print('Done!')
