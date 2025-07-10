# Introduction

The pipeline designed for processing MORF library screening dataset using bulk RNA-seq and scRNA-seq.

# Requirements

Requires Python ≥ 3.8 and the following libraries:
```bash
# create a new conda enviroment
pip install numpy scipy matplotlib pandas seaborn biopython
```

# Usage

## TF barcode counting

Get the barcode abundance profiles.

### 1. Input

The script accepts the following required command-line arguments:

| Argument            | Type     | Description                                                                                                  | Default    |
|---------------------|----------|--------------------------------------------------------------------------------------------------------------|------------|
| `-fq`               | `str`    | fastq file                                                                                                   |            |
| `-b`                | `str`    | Path to barcode list in CSV format ,example /mnt/dfc_data2/project/linyusen/project/81_MORF/data/barcode.csv | None       |
| `-o`                | `str`    | Output directory to save results                                                                             |            |
| `-thread`           | `str`    | Number of threads to use for processing                                                                      |            |
| `-key`              | `str`    | Anchor key sequence located before the barcode (e.g., GAAAGGACGA)                                            | GAAAGGACGA |
| `-KEY_REGION_START` | `int`    | Start position (0-based) of the region where the key sequence is expected                                    | 25         |
| `-KEY_REGION_END`   | `int`    | End position (0-based, exclusive) of the region where the key sequence is expected                           | 50         |
| `-BARCODE_LENGTH`   | `int`    | Length of the barcode sequence to extract following the key                                                  | 24         |

### 2. Output

After completion, the output directory will contain:

TF_Count.csv：A table containing barcode counts and CPM (counts per million) for each TF

TF_Count_Distribution.png：	A plot showing the overall barcode count distribution

Log-Scale_Abundance_Histogram.png：A histogram of barcode abundances on a log scale

### 3. Example

Follow these steps to quickly analyse library diversity:
```angular2html
python tf_find_multi.py \
  -fq xxxx.fastq \
  -b ./barcode.csv \
  -o ./xxxx \
  -thread 64 \
  -key GAAAGGACGA \
  -KEY_REGION_START 25 \
  -KEY_REGION_END 50 \
  -BARCODE_LENGTH 24 \
```

## Comparision analysis

This script is designed for comparing two groups of MORF libraries based on their barcode abundance profiles.

### 1. Input

The script accepts the following required command-line arguments:

| Argument            | Type     | Description                      | Default    |
|---------------------|----------|----------------------------------|------------|
| `-i`                | `str`    | sample file                      |            |
| `-o`                | `str`    | Output directory to save results |            |
| `-b`                | `str`    | tf informa                       |            |

### 2. Output

After completion, the output directory will contain:

barcode_scatter_plot.png：A scatter plot comparing barcode CPM values between the two groups

barcode_heatmap.png：A heatmap showing clustering of TFs based on barcode abundance

Rank-Rank_Plot.png：A rank-rank plot visualizing the TFs overall consistency and shifts in ranks

compare.info.csv: A summary table including mean CPM, log2 fold changes

### 3. Example

Follow these steps to quickly analyse library diversity:
```angular2html
python tf_find_compare.py \
  -i samples.lst \
  -b xxx.txt \
  -o ./xxxx \
```
> ⚠️ Each path in samples.lst corresponds to the output directory of tf_find_multi.py


# Please Cite

If you use this script or parts of it in your research or project, please cite the repository or acknowledge the author appropriately. A suggested citation format:


# Contributors

Jiajian ZHOU, Yusen Lin, Tingting Deng and Meting.


# License

This project is licensed under the [MIT License](LICENSE.txt).  
You are free to use, modify, and distribute it with attribution.
