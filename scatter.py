import pandas as pd
import os
import numpy as np
from collections import Counter
import random
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter
import matplotlib as mpl
from scipy.stats import mannwhitneyu
import seaborn as sns
from pathlib import Path
import yaml

import scifont
mpl.rcParams['font.family'] = 'Arimo'
mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['font.weight']      = 'bold'
# mpl.rcParams['axes.titleweight'] = 'bold'
# mpl.rcParams['axes.labelweight'] = 'bold'

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
    typ_disp = {'NP': 'no protein', 'CR': 'crRNA', 'DT': 'dict'}[typ]
    return f'{prot} ({typ_disp}, rep {rep})'

cfg = load_config("config.yaml")
m   = load_manifest(cfg)

sOutputDir = cfg["output_dir"]
os.makedirs(sOutputDir, exist_ok=True)

clv1 = m[(m.replicate == "1") & (m.type != "DT")]
clv2 = m[(m.replicate == "2") & (m.type != "DT")]

lClvFile1  = clv1["file"].tolist()
lClvFile2  = clv2["file"].tolist()
lNames1    = clv1["label"].tolist()
lNames2    = clv2["label"].tolist()
sDictFile1 = m[(m.replicate == "1") & (m.type == "DT")]["file"].iloc[0]
sDictFile2 = m[(m.replicate == "2") & (m.type == "DT")]["file"].iloc[0]

dfDict1 = pd.read_table(sDictFile1, names=['name', 'seq'], index_col=False)
dfDict2 = pd.read_table(sDictFile2, names=['name', 'seq'], index_col=False)

dfDict1['count'] = dfDict1['name'].apply(lambda x: int(x.split('-')[-1]))
dfDict2['count'] = dfDict2['name'].apply(lambda x: int(x.split('-')[-1]))

print("Printing: dfDict1")
print(dfDict2)
print("\n\n\n")
print("Printing: dfDict2")
print(dfDict2)
print("\n\n\n")

nDictCountFilter1 = 100
nDictCountFilter2 = 50
dfDict1 = dfDict1[dfDict1['count'] > nDictCountFilter1]
dfDict2 = dfDict2[dfDict2['count'] > nDictCountFilter2]

print("Printing: dfDict1")
print(dfDict2)
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
    ax.yaxis.set_major_locator(ticker.LogLocator(base=2))
    ax.yaxis.set_major_locator(FixedLocator([16, 64, 256, 1024, 4096]))
    ax.yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.savefig(os.path.join(sOutputDir, 'dict_depth.png'), dpi=1200)
plt.close()