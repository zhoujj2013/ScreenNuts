#%%
import os
import csv
import argparse

parser = argparse.ArgumentParser()
# 添加参数
parser.add_argument('-input_dir1', type=str, required=True,
                    help='Directory containing barcode counts for sample1')
parser.add_argument('-input_dir2', type=str, required=True,
                    help='Directory containing barcode counts for sample2')
parser.add_argument('-o', type=str, required=True,
                    help='Output directory to save plots and results')
parser.add_argument('-sample1', type=str, default='sample1',
                    help='Name of the first sample (for labeling in plots)')
parser.add_argument('-sample2', type=str, default='sample2',
                    help='Name of the second sample (for labeling in plots)')
args = parser.parse_args()

input_dir1 = args.input_dir1
input_dir2 = args.input_dir2
output_dir = args.o
sample1 = args.sample1
sample2 = args.sample2


# input_dir1 = '/mnt/dfc_data2/project/linyusen/project/81_MORF/data/SRR22409540'
# input_dir2 = '/mnt/dfc_data2/project/linyusen/project/81_MORF/data/SRR22409541'
# output_dir = '/mnt/dfc_data2/project/linyusen/project/81_MORF/data/compare'
# sample1 = 'SRR22409540'
# sample2 = 'SRR22409541'


TF_Count_file1 = os.path.join(input_dir1,'TF_Count.csv')
TF_Count_file2 = os.path.join(input_dir2,'TF_Count.csv')

os.makedirs(output_dir,exist_ok=True)
def read_data(file_path):
    f = open(file_path)
    r = csv.reader(f)
    r.__next__()
    r.__next__()
    data = {}
    for i in r:
        barcode,count,cpm = i
        count = int(count)
        cpm = float(cpm)
        data[barcode] = cpm
    return data
TF_Count1 = read_data(TF_Count_file1)
TF_Count2 = read_data(TF_Count_file2)

#%%
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

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
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
import matplotlib.pyplot as plt

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


f = open(os.path.join(output_dir,'compare.info.csv'),'w')
w = csv.writer(f)
w.writerow([f'# Dropout TF Counts:{dropout_counts}'])
w.writerow(['barcode',f'{sample1}_count',f'{sample1}_count',f'log2FC({sample1}/{sample2})',f'{sample1}_rank',f'{sample1}_rank'])
for k in common_keys:
    w.writerow([k,TF_Count1[k],TF_Count2[k],log2fc[k],rank_dict[k][0],rank_dict[k][1]])
f.close()