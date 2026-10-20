import pandas as pd
import os
from pathlib import Path
import yaml
import numpy as np
from collections import Counter
import random
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter
import matplotlib as mpl
from scipy.stats import mannwhitneyu
import seaborn as sns
import scifont

mpl.rcParams['font.family'] = 'Arimo'
mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.weight']      = 'bold'
# mpl.rcParams['axes.titleweight'] = 'bold'
# mpl.rcParams['axes.labelweight'] = 'bold'

# Analysis constants
nDictCountFilter1 = 100
nDictCountFilter2 = 50
nClvCountFilter   = 20
nMinInsertLen     = 15
nMaxInsertLen     = 25
nTopNumber        = 4000

def load_config(path="config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)

def load_manifest(cfg):
    df = pd.read_csv(cfg["sample_sheet"], dtype={"replicate": str})
    ppseq = Path(cfg["intermediate_dir"]) / cfg["preprocessed_seq_dir"]
    df["file"]  = df.apply(lambda r: str(ppseq /
                       f"{r.replicate}_{r.protein}-{r.type}.collapsed.tsv"), axis=1)
    df["label"] = df["protein"] + "-" + df["type"]
    return df

def display_name(sClvName):
    base = sClvName.split('/')[-1].removesuffix('.collapsed.tsv')   # '1_Lbu-NP'
    rep, rest = base.split('_', 1)                                  # '1', 'Lbu-NP'
    prot, typ = rest.split('-', 1)                                  # 'Lbu', 'NP'
    typ_disp = {'NP': 'noProtein', 'CR': 'withProtein', 'DT': 'dict'}[typ]
    return f'{prot}-{typ_disp}'

cfg = load_config("config.yaml")
m   = load_manifest(cfg)

sOutputDir = cfg["output_dir"]
os.makedirs(sOutputDir, exist_ok=True)
os.makedirs(os.path.join(sOutputDir, "supp_fig"), exist_ok=True)

clv1 = m[(m.replicate == "1") & (m.type != "DT")]
clv2 = m[(m.replicate == "2") & (m.type != "DT")]

lClvFile1  = clv1["file"].tolist()
lClvFile2  = clv2["file"].tolist()
lNames1    = clv1["label"].tolist()
lNames2    = clv2["label"].tolist()
sDictFile1 = m[(m.replicate == "1") & (m.type == "DT")]["file"].iloc[0]
sDictFile2 = m[(m.replicate == "2") & (m.type == "DT")]["file"].iloc[0]

### Process Dictionary; Supp-Fig-3 | Sequencing Depth of Dictionary ###

dfDict1 = pd.read_table(sDictFile1, names=['name', 'seq'], index_col=False)
dfDict2 = pd.read_table(sDictFile2, names=['name', 'seq'], index_col=False)

dfDict1['count'] = dfDict1['name'].apply(lambda x: int(x.split('-')[-1]))
dfDict2['count'] = dfDict2['name'].apply(lambda x: int(x.split('-')[-1]))

print("Printing: dfDict1")
print(dfDict1)
print("\n\n\n")
print("Printing: dfDict2")
print(dfDict2)
print("\n\n\n")

dfDict1 = dfDict1[dfDict1['count'] > nDictCountFilter1]
dfDict2 = dfDict2[dfDict2['count'] > nDictCountFilter2]

print("Printing: dfDict1")
print(dfDict1)
print("\n\n\n")
print("Printing: dfDict2")
print(dfDict2)
print("\n\n\n")

print(f'Dictionary filter rep1: {nDictCountFilter1} reads, '
      f'{len(dfDict1)} unique 35mers, total depth {dfDict1["count"].sum()}, '
      f'median {dfDict1["count"].median()}')
print(f'Dictionary filter rep2: {nDictCountFilter2} reads, '
      f'{len(dfDict2)} unique 35mers, total depth {dfDict2["count"].sum()}, '
      f'median {dfDict2["count"].median()}')

fig, axes = plt.subplots(1, 2, figsize=(5, 4))
for ax, df, filt, title in [
    (axes[0], dfDict1, nDictCountFilter1, 'Rep1'),
    (axes[1], dfDict2, nDictCountFilter2, 'Rep2'),
]:
    sns.boxplot(y='count', data=df, ax=ax, width=0.5, flierprops={"marker": "d", "markerfacecolor": "black", "markeredgecolor": "black"})
    ax.set_title(f'{title}, Dictionary\nunique sequences={len(df):,}\ntotal depth={df["count"].sum():,}',
                 fontsize=10)
    ax.set_ylabel('Number of reads per sequence')
    ax.axhline(y=filt, color='blue', linestyle='--')
    ax.set_ylim([16, 8000])
    ax.set_yscale('log', base=2)
    ax.yaxis.set_major_locator(FixedLocator([16, 64, 256, 1024, 4096]))
    ax.yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp3_dict_depth.png'), dpi=1200)
plt.close()

### Map Cleavage Fragments to Dictionary ###

def map_replicate(label_rep, clv_files, df_dict,
                  n_clv_count_filter=nClvCountFilter,
                  min_insert_len=nMinInsertLen, max_insert_len=nMaxInsertLen):
    mapped = []
    for sClvFile in clv_files:
        sClvName = sClvFile.split('/')[-1].replace('.collapsed.tsv', '')
        dfClv = pd.read_table(sClvFile, names=['name', 'seq'], index_col=False)
        dfClv['count'] = dfClv['name'].apply(lambda x: int(x.split('-')[-1]))
        dfClv = dfClv[dfClv['count'] > n_clv_count_filter]
        dfClv['insertlen'] = dfClv['seq'].apply(len)
        dfClv = dfClv[(dfClv['insertlen'] >= min_insert_len) &
                      (dfClv['insertlen'] <  max_insert_len)]

        lAllDf = []
        for L in range(min_insert_len, max_insert_len):
            dfClvOne  = dfClv[dfClv['insertlen'] == L]
            dfDictOne = df_dict.copy()
            dfDictOne['5pSeq']     = dfDictOne['seq'].apply(lambda x: x[:L])
            dfDictOne['8nt_motif'] = dfDictOne['seq'].apply(lambda x: x[L-4:L+4])
            dfDictOne['3nt_motif'] = dfDictOne['seq'].apply(lambda x: x[L-2:L+1])
            dfDictOne['2nt_motif'] = dfDictOne['seq'].apply(lambda x: x[L-1:L+1])

            dfClvOne  = dfClvOne.rename(columns={'name':'clv_name',  'seq':'clv_seq',  'count':'clv_count'})
            dfDictOne = dfDictOne.rename(columns={'name':'dict_name','seq':'dict_seq','count':'dict_count'})

            merged = pd.merge(dfClvOne, dfDictOne, how='left', left_on='clv_seq', right_on='5pSeq').dropna()
            merged['clvOverDict'] = merged['clv_count'] / merged['dict_count'].astype(float)
            lAllDf.append(merged)

        dfClvConcat = pd.concat(lAllDf)
        dfClvConcat = dfClvConcat.sort_values('clvOverDict', ascending=False).reset_index(drop=True)
        dfClvConcat['rank'] = dfClvConcat['clvOverDict'].rank(method='first', ascending=False).astype(int)
        mapped.append([label_rep, sClvName, dfClvConcat])
    return mapped

lMappedClv1 = map_replicate('Rep1', lClvFile1, dfDict1)
lMappedClv2 = map_replicate('Rep2', lClvFile2, dfDict2)

print("Printing: lMappedClv1")
print(lMappedClv1)
print("\n\n\n")
print("Printing: lMappedClv2")
print(lMappedClv2)
print("\n\n\n")

### Supp-Fig-4 | Sequencing Depth of Cleavage Libraries (distribution of reads per sequence; mapped reads on dictionary) — Rep1 ###

fig = plt.figure(figsize=(10, 8))
for i, (sLabelRep, sClvName, dfClvConcat) in enumerate(lMappedClv1):
    ax = fig.add_subplot(3, 4, i+1)
    sns.boxplot(y='clv_count', data=dfClvConcat, ax=ax, width=0.5)
    ax.set_title(f'{display_name(sClvName)}\nunique sequences={len(dfClvConcat):,}\n'
                 f'total depth={dfClvConcat["clv_count"].sum():,}',
                 fontsize=10)
    ax.set_ylabel('Number of reads per sequence')
    ax.axhline(y=nClvCountFilter, color='blue', linestyle='--')
    ax.set_ylim([16, 8000])
    ax.set_yscale('log', base=2)
    ax.yaxis.set_major_locator(FixedLocator([16, 64, 256, 1024, 4096]))
    ax.yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.suptitle('Rep1', y=0.99)
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp4_rep1_clv_depth.png'), dpi=1200)
plt.close()

### Supp-Fig-4 | Sequencing Depth of Cleavage Libraries (distribution of reads per sequence; mapped reads on dictionary) — Rep2 ###

fig = plt.figure(figsize=(10, 8))
for i, (sLabelRep, sClvName, dfClvConcat) in enumerate(lMappedClv2):
    ax = fig.add_subplot(3, 4, i+1)
    sns.boxplot(y='clv_count', data=dfClvConcat, ax=ax, width=0.5)
    ax.set_title(f'{display_name(sClvName)}\nunique sequences={len(dfClvConcat):,}\n'
                 f'total depth={dfClvConcat["clv_count"].sum():,}',
                 fontsize=10)
    ax.set_ylabel('Number of reads per sequence')
    ax.axhline(y=nClvCountFilter, color='blue', linestyle='--')
    ax.set_ylim([16, 8000])
    ax.set_yscale('log', base=2)
    ax.yaxis.set_major_locator(FixedLocator([16, 64, 256, 1024, 4096]))
    ax.yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.suptitle('Rep2', y=0.99)
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp4_rep2_clv_depth.png'), dpi=1200)
plt.close()

### Supp-Fig-5 | Distribution of relative cleavage efficiencies of cleaved sequences — Rep1 ###

fig = plt.figure(figsize=(10, 8))
for i, (sLabelRep, sClvName, dfClvConcat) in enumerate(lMappedClv1):
    ax = fig.add_subplot(3, 4, i+1)
    sns.scatterplot(x='rank', y='clvOverDict', data=dfClvConcat,
                    alpha=1, s=5, color='tab:blue',
                    edgecolor='tab:blue', linewidth=0.1, ax=ax)
    ax.set_title(display_name(sClvName))
    ax.set_ylabel('Relative cleavage efficiency')
    ax.axvline(x=nTopNumber, color='blue', linestyle='--')

plt.tight_layout()
plt.suptitle('Rep1', y=0.99)
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp5_rep1_clv_eff.png'), dpi=1200)
plt.close()

### Supp-Fig-5 | Distribution of relative cleavage efficiencies of cleaved sequences — Rep2 ###

fig = plt.figure(figsize=(10, 8))
for i, (sLabelRep, sClvName, dfClvConcat) in enumerate(lMappedClv2):
    ax = fig.add_subplot(3, 4, i+1)
    sns.scatterplot(x='rank', y='clvOverDict', data=dfClvConcat,
                    alpha=1, s=5, color='tab:blue',
                    edgecolor='tab:blue', linewidth=0.1, ax=ax)
    ax.set_title(display_name(sClvName))
    ax.set_ylabel('Relative cleavage efficiency')
    ax.axvline(x=nTopNumber, color='blue', linestyle='--')

plt.tight_layout()
plt.suptitle('Rep2', y=0.99)
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp5_rep2_clv_eff.png'), dpi=1200)
plt.close()
