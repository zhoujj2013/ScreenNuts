import os
import argparse
import csv
from Bio import SeqIO
from multiprocessing import Pool
from itertools import islice
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

parse = argparse.ArgumentParser()
parse.add_argument('-fq',type=str,help='fastq file',required=True)
parse.add_argument('-b',type=str,help='barcode list,csv file,example /mnt/dfc_data2/project/linyusen/project/81_MORF/barcode_info.csv',required=True)
parse.add_argument('-prefix',type=str,help='output name prefix',required=True)
parse.add_argument('-thread',type=int,help='number of threads to precess',required=True)
parse.add_argument('-key',type=str,default='GAAAGGACGA',help='a string before barcode,example GAAAGGACGA')
parse.add_argument('-KEY_REGION_START',type=int,default=25)
parse.add_argument('-KEY_REGION_END',type=int,default=50)
parse.add_argument('-BARCODE_LENGTH',type=int,default=24)
args = parse.parse_args()

fastq_file = args.fq
input_file = args.b
prefix = args.prefix
KEY = args.key
KEY_REGION_START = args.KEY_REGION_START
KEY_REGION_END = args.KEY_REGION_END
BARCODE_LENGTH = args.BARCODE_LENGTH
thread = args.thread




# KEY_REGION_START = 25  # start index of key region
# KEY_REGION_END = 50  # end index of key region
# BARCODE_LENGTH = 24
# KEY = "GAAAGGACGA"  # identifies sequence before barcode to determine barcode position
# input_file = '/mnt/dfc_data2/project/linyusen/project/81_MORF/barcode_info.csv'
# fastq_file = '/mnt/dfc_data2/project/linyusen/database/81_MORF/test_data/SRR22409537_1.fastq'
# output_dir = '/mnt/dfc_data2/project/linyusen/project/81_MORF/data/SRR22409537'
# thread = 32



#%%
cwd = os.getcwd()
output_dir = os.path.join(cwd,prefix)
os.makedirs(output_dir,exist_ok=True)
with open(input_file, 'r') as infile:
    barcode_dict = {}
    reader = csv.reader(infile)
    reader.__next__()
    for rows in reader:
        barcode_dict[rows[3]] = 0

barcode_tf_dict = {}
with open(input_file, 'r') as infile:
    reader = csv.reader(infile)
    reader.__next__()
    for rows in reader:
        barcode_tf_dict[rows[3]] = {}
        barcode_tf_dict[rows[3]]['Name'] = rows[0]
        barcode_tf_dict[rows[3]]['RefSeqGeneName'] = rows[1]
        barcode_tf_dict[rows[3]]['RefSeqandGencodeID'] = rows[2]
        barcode_tf_dict[rows[3]]['BarcodeSequence'] = rows[3]
#%%

def process_chunk(records_chunk):
    local_dict = {}
    local_stats = {'perfect': 0, 'nonperfect': 0, 'notfound': 0}

    for record in records_chunk:
        read_sequence = str(record.seq).upper()
        key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
        key_index = key_region.find(KEY)
        if key_index >= 0:
            start_index = key_index + KEY_REGION_START + len(KEY)
            barcode = read_sequence[start_index:(start_index + BARCODE_LENGTH)]
            if barcode in barcode_dict:
                local_dict[barcode] = local_dict.get(barcode, 0) + 1
                local_stats['perfect'] += 1
            else:
                local_stats['nonperfect'] += 1
        else:
            local_stats['notfound'] += 1

    return local_dict, local_stats

def chunk_iterable(iterable, chunk_size):
    """Yield chunks of `chunk_size` from iterable"""
    iterator = iter(iterable)
    while True:
        chunk = list(islice(iterator, chunk_size))
        if not chunk:
            break
        yield chunk



num_processes = thread
chunk_size = 1000  # 每个进程处理1万条 read

print('Precessing ... ')
results = []
with Pool(processes=num_processes) as pool:
    with open(fastq_file, "r") as handle:
        reader = SeqIO.parse(handle, "fastq")
        for chunk in chunk_iterable(reader, chunk_size):
            results.append(pool.apply_async(process_chunk, args=(chunk,)))
    pool.close()
    pool.join()
#%%

# 合并所有结果

final_dict = Counter()
total_stats = {'perfect': 0, 'nonperfect': 0, 'notfound': 0}

for r in results:
    partial_dict, stats = r.get()
    final_dict.update(partial_dict)
    for k in total_stats:
        total_stats[k] += stats[k]
#%%
# 更新 barcode_dict
for k in barcode_dict:
    barcode_dict[k] = final_dict.get(k, 0)
total_count = np.sum(list(total_stats.values()))
#%%

barcode_count_list = []
for i in barcode_dict:
    barcode_count_list.append(barcode_dict[i]+1)
barcode_count_list = np.array(barcode_count_list)
log_abundance = np.log10(barcode_count_list)
plt.hist(log_abundance, bins=20, color='skyblue', edgecolor='black')
plt.xlabel("log10(Abundance)")
plt.ylabel("Number of Barcodes")
plt.title("Log-Scale Abundance Histogram")
plt.savefig(os.path.join(output_dir,'Log-Scale_Abundance_Histogram.png'))
#%%
sorted_data = dict(sorted(barcode_dict.items(), key=lambda item: item[1], reverse=True))
y = list(sorted_data.values())
x = range(len(y))
max_idx = y.index(max(y))
min_idx = y.index(min(y))
max_kmer = list(sorted_data.keys())[max_idx]
min_kmer = list(sorted_data.keys())[min_idx]
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='K-mer Count')
plt.fill_between(x, y, alpha=0.3)
plt.xlabel("Ranked K-mers")
plt.ylabel("Count")
plt.title(f"TF Count Distribution")
plt.xticks([])
plt.text(max_idx, y[max_idx], f"Max: {max_kmer} ({y[max_idx]})",
         ha='left', va='bottom', fontsize=9, color='red')
plt.text(min_idx, y[min_idx], f"Min: {min_kmer} ({y[min_idx]})",
         ha='right', va='top', fontsize=9, color='blue')
plt.tight_layout()
plt.savefig(os.path.join(output_dir,'TF_Count_Distribution.png'))
#%%
barcode_count_list = []
for i in barcode_dict:
    barcode_count_list.append(barcode_dict[i])
skewness = np.percentile(barcode_count_list, 90) / np.percentile(barcode_count_list, 10)
#%%
with open(os.path.join(output_dir,'TF_Count.txt'), 'w') as csvfile:
    mywriter = csv.writer(csvfile, delimiter='\t')
    mywriter.writerow(['Barcode','Name','RefSeqGeneName','RefSeqandGencodeID',f'{prefix}_count',f'{prefix}_CPM'])
    for barcode in barcode_dict:
        count = sorted_data[barcode]
        mywriter.writerow([barcode,barcode_tf_dict[barcode]['Name'],barcode_tf_dict[barcode]['RefSeqGeneName'],barcode_tf_dict[barcode]['RefSeqandGencodeID'],count,count/total_count*1000000])


with open(os.path.join(output_dir,'readcount_stat_header.txt'), 'w') as csvfile:
    mywriter = csv.writer(csvfile, delimiter='\t')
    mywriter.writerow(['sample','total', 'barcode count','nonperfect','notfound','skewness'])

with open(os.path.join(output_dir,'readcount_stat.txt'), 'w') as csvfile:
    mywriter = csv.writer(csvfile, delimiter='\t')
    mywriter.writerow([prefix,total_count, total_stats['perfect'],total_stats['nonperfect'],total_stats['notfound'],skewness])
print('Done!')
