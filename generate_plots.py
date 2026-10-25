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

# at the top of generate_plots.py
try:
    cfg_path = snakemake.params.config_path
except NameError:
    cfg_path = "config.yaml"                     # when run as `python generate_plots.py`
cfg = load_config(cfg_path)

m   = load_manifest(cfg)

sOutputDir = cfg["output_dir"]
os.makedirs(sOutputDir, exist_ok=True)
os.makedirs(os.path.join(sOutputDir, "main_fig"), exist_ok=True)
os.makedirs(os.path.join(sOutputDir, "supp_fig"), exist_ok=True)
logo_input_dir = os.path.join(cfg.get("intermediate_dir", "intermediate_files"), "logo_input")
os.makedirs(logo_input_dir, exist_ok=True)

### Manifest-derived replication structure ###
REPS = sorted(m["replicate"].unique())                 # e.g. ['1'] or ['1', '2']

# per-replicate cleavage sample lists (only replicates that actually have non-DT rows)
lClvFiles = {rep: m[(m.replicate == rep) & (m.type != "DT")]["file"].tolist()
             for rep in REPS}

# per-replicate dictionary file (only replicates that have a DT row)
dict_files = {rep: m[(m.replicate == rep) & (m.type == "DT")]["file"].iloc[0]
              for rep in REPS
              if ((m.replicate == rep) & (m.type == "DT")).any()}

print(f"replicates: {REPS}")
print(f"dictionary files: {dict_files}")

# Rep1 aliases (for backwards-compat with the code below)
lClvFile1 = lClvFiles.get("1", [])
lNames1   = [m[(m.replicate == "1") & (m.file == f)]["label"].iloc[0] for f in lClvFile1]

# Rep2 aliases — empty when only 1 replicate
lClvFile2 = lClvFiles.get("2", [])
lNames2   = [m[(m.replicate == "2") & (m.file == f)]["label"].iloc[0] for f in lClvFile2]

# Guarded dictionary files — None when the replicate has no dict
sDictFile1 = dict_files.get("1")
sDictFile2 = dict_files.get("2")

### Process Dictionary; Supp-Fig-3a,b | Sequencing Depth of Dictionary ###
def _load_dict(path, filt):
    if path is None:
        return None
    df = pd.read_table(path, names=['name', 'seq'], index_col=False)
    df['count'] = df['name'].apply(lambda x: int(x.split('-')[-1]))
    df = df[df['count'] > filt]
    return df

dfDict1 = _load_dict(sDictFile1, nDictCountFilter1)
dfDict2 = _load_dict(sDictFile2, nDictCountFilter2)

for label, df, filt in [("Rep1", dfDict1, nDictCountFilter1),
                        ("Rep2", dfDict2, nDictCountFilter2)]:
    if df is None:
        print(f"{label}: no dictionary")
    else:
        print(f"{label}: {filt} read filter → "
              f"{len(df)} unique 35mers, "
              f"total depth {df['count'].sum()}, "
              f"median {df['count'].median()}")

_panels = []
if dfDict1 is not None:
    _panels.append((dfDict1, nDictCountFilter1, "Rep1"))
if dfDict2 is not None:
    _panels.append((dfDict2, nDictCountFilter2, "Rep2"))

fig, axes = plt.subplots(1, len(_panels), figsize=(2.5 * len(_panels), 4),
                         squeeze=False)
for ax, (df, filt, title) in zip(axes[0], _panels):
    sns.boxplot(y='count', data=df, ax=ax, width=0.5,
                flierprops={"marker": "d", "markerfacecolor": "black",
                            "markeredgecolor": "black"})
    ax.set_title(f'{title}, Dictionary\nunique sequences={len(df):,}\n'
                 f'total depth={df["count"].sum():,}', fontsize=10)
    ax.set_ylabel('Number of reads per sequence')
    ax.axhline(y=filt, color='blue', linestyle='--')
    ax.set_ylim([16, 8000])
    ax.set_yscale('log', base=2)
    ax.yaxis.set_major_locator(FixedLocator([16, 64, 256, 1024, 4096]))
    ax.yaxis.set_major_formatter(ScalarFormatter())

plt.tight_layout()
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp3ab_dict_depth.png'), dpi=1200)
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

lMappedClv1 = map_replicate('Rep1', lClvFile1, dfDict1) if dfDict1 is not None and lClvFile1 else []
lMappedClv2 = map_replicate('Rep2', lClvFile2, dfDict2) if dfDict2 is not None and lClvFile2 else []
lAllMappedClv = lMappedClv1 + lMappedClv2

def _sample_key(sClvName):
    # '1_Lbu-NP' -> 'Lbu-NP'
    return sClvName.split('_', 1)[1]

by_key1 = {_sample_key(n): (rep, n, df) for (rep, n, df) in lMappedClv1}
by_key2 = {_sample_key(n): (rep, n, df) for (rep, n, df) in lMappedClv2}

missing1 = sorted(by_key1.keys() - by_key2.keys())
missing2 = sorted(by_key2.keys() - by_key1.keys())
if missing1: print(f'Rep1-only samples: {missing1}')
if missing2: print(f'missing in Rep1: {missing2}')

lZipMappedClv = [(by_key1[k], by_key2[k]) for k in sorted(by_key1.keys() & by_key2.keys())]

print("Printing: lMappedClv1")
print(lMappedClv1)
print("\n\n\n")
print("Printing: lMappedClv2")
print(lMappedClv2)
print("\n\n\n")

### Supp-Fig-4a | Sequencing Depth of Cleavage Libraries (distribution of reads per sequence; mapped reads on dictionary) — Rep1 ###

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
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp4a_rep1_clv_depth.png'), dpi=1200)
plt.close()

### Supp-Fig-4b | Sequencing Depth of Cleavage Libraries (distribution of reads per sequence; mapped reads on dictionary) — Rep2 ###

if lMappedClv2:
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
    plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp4b_rep2_clv_depth.png'), dpi=1200)
    plt.close()
else:
    print("no Rep2 — skipping Supp-Fig-4b")

### Supp-Fig-5a | Distribution of relative cleavage efficiencies of cleaved sequences — Rep1 ###

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
plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp5a_rep1_clv_eff.png'), dpi=1200)
plt.close()

### Supp-Fig-5b | Distribution of relative cleavage efficiencies of cleaved sequences — Rep2 ###

if lMappedClv2:
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
    plt.savefig(os.path.join(sOutputDir, 'supp_fig', 'supp5b_rep2_clv_eff.png'), dpi=1200)
    plt.close()
else:
    print("no Rep2 — skipping Supp-Fig-5b")
    
### Supp-Fig-6a | Cleavage efficiency by rank group — one figure per CR sample ###

S6A_BINS   = cfg["supp6a"]["bins"]
S6A_LABELS = cfg["supp6a"]["labels"]
S6A_PALETTE = dict(zip(S6A_LABELS, cfg["supp6a"]["colors"]))

# one figure per (protein) that has a CR sample — derive from the manifest
cr_proteins = sorted({row.protein for _, row in m.iterrows() if row.type == "CR"})

for prot in cr_proteins:
    target_label = f"{prot}-withProtein"

    fig = plt.figure(figsize=(4.5, 6))
    i = 0
    for sLabelRep, sClvName, dfClvConcat in lAllMappedClv:
        if display_name(sClvName) != target_label:
            continue
        i += 1

        dfClvConcat = dfClvConcat.copy()                    # don't mutate the cached df
        dfClvConcat['group'] = pd.cut(dfClvConcat['rank'],
                                       bins=S6A_BINS,
                                       labels=S6A_LABELS,
                                       right=True)
        counts = dfClvConcat['group'].value_counts()

        ax = fig.add_subplot(1, 2, i)
        sns.boxplot(
            data=dfClvConcat, x='group', y='clvOverDict',
            hue='group', palette=S6A_PALETTE, legend=False,
            dodge=False, width=0.5,
            flierprops={"marker": "d", "markerfacecolor": "black",
                        "markeredgecolor": "black"},
            ax=ax,
        )
        ax.set_xlabel('')
        xticklabels = [t.get_text() for t in ax.get_xticklabels()]
        ax.set_xticks(range(len(xticklabels)))
        ax.set_xticklabels([f'{lbl} (n={counts[lbl]})' for lbl in xticklabels], rotation=90)
        ax.set_title(sLabelRep)
        ax.set_ylabel('Relative cleavage efficiency')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(sOutputDir, 'supp_fig',
                             f'supp6a_{prot}_clv-eff.pdf'), dpi=1200)
    plt.close()

### Fig-2c,e | Motif-Cleavage Efficiency Distribution ###

def plot_motif_scatter(mapped_pairs, target, motif_type,
                       n_top, n_bottom, xlim, ylim, out_path,
                       highlight_motif=None):
    """Rep1 vs Rep2 mean clv/dict scatter for one (protein, type) sample.

    mapped_pairs : list of ((rep1, sClvName1, df1), (rep2, sClvName2, df2)) tuples
    target       : manifest label, e.g. 'Lbu-CR', 'Psp-CR'
    motif_type   : '2nt_motif' or '3nt_motif'
    highlight_motif : optional motif to force-highlight (e.g. 'CT' for Lwa)
    """
    for xRep1, xRep2 in mapped_pairs:
        sample_key = xRep1[1].split('_', 1)[1]          # 'Lbu-CR'
        if sample_key != target:
            continue

        sSampleName = display_name(xRep1[1])            # 'Lbu-withProtein' for the title
        dfClvRep1, dfClvRep2 = xRep1[2], xRep2[2]

        ave1 = dfClvRep1[[motif_type, 'clvOverDict']].groupby(motif_type)['clvOverDict'].mean().reset_index()
        ave1.columns = [motif_type, 'mean_clvOverDict']
        ave2 = dfClvRep2[[motif_type, 'clvOverDict']].groupby(motif_type)['clvOverDict'].mean().reset_index()
        ave2.columns = [motif_type, 'mean_clvOverDict']

        merged = ave1.merge(ave2, on=motif_type, suffixes=('_rep1', '_rep2')).dropna()
        merged['avg'] = merged[['mean_clvOverDict_rep1', 'mean_clvOverDict_rep2']].mean(axis=1)

        top    = merged.nlargest(n_top, 'avg')
        bottom = merged.nsmallest(n_bottom, 'avg')
        others = merged.drop(index=top.index.tolist() + bottom.index.tolist())

        fig, ax = plt.subplots(figsize=(3, 3))
        sns.scatterplot(data=others, x='mean_clvOverDict_rep1', y='mean_clvOverDict_rep2',
                        color='lightgray', s=40, alpha=0.8, edgecolor='white', linewidth=0.5, ax=ax)
        sns.scatterplot(data=bottom, x='mean_clvOverDict_rep1', y='mean_clvOverDict_rep2',
                        color='tab:blue', s=50, edgecolor='black', linewidth=0.5,
                        label=f'Bottom {n_bottom} (lowest)', ax=ax)

        label_groups = [(top, 'tab:red'), (bottom, 'tab:blue')]
        if highlight_motif is not None:
            df_interest = others[others[motif_type] == highlight_motif]
            if not df_interest.empty:
                sns.scatterplot(data=df_interest, x='mean_clvOverDict_rep1', y='mean_clvOverDict_rep2',
                                color='tab:blue', s=50, edgecolor='black', linewidth=0.5, ax=ax)
                label_groups.append((df_interest, 'tab:blue'))

        sns.scatterplot(data=top, x='mean_clvOverDict_rep1', y='mean_clvOverDict_rep2',
                        color='tab:red', s=50, edgecolor='black', linewidth=0.5,
                        label=f'Top {n_top} (highest)', ax=ax)

        for df, color in label_groups:
            for _, row in df.iterrows():
                ax.text(row['mean_clvOverDict_rep1'],
                        row['mean_clvOverDict_rep2'] + 0.02,
                        row[motif_type],
                        fontsize=10, fontweight='bold', color=color, ha='left', va='bottom')

        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xlabel('Rep 1', fontsize=8); ax.set_ylabel('Rep 2', fontsize=8)
        ax.set_title(f'{sSampleName} mean clv/dict', fontsize=8)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_path, dpi=1200)
        plt.close(fig)
        return
    raise ValueError(f'target {target!r} not found in mapped_pairs')

### Fig-2c / Fig-2e | Motif-Cleavage Efficiency Distribution ###

MOTIF_FIGS = cfg.get("motif_figs") or {}

for target, p in MOTIF_FIGS.items():
    plot_motif_scatter(
        lZipMappedClv, target,
        motif_type      = p["motif_type"],
        n_top           = p["n_top"],
        n_bottom        = p["n_bottom"],
        xlim            = p["xlim"],
        ylim            = p["ylim"],
        out_path        = os.path.join(sOutputDir, "main_fig", p["out"]),
        highlight_motif = p.get("highlight"),           # None if key absent or null
    )

### Fig-2f | Motif-cleavage efficiency vs reference, for one target sample ###

f2f = cfg.get("fig2f")

if f2f:
    sTargetSampleName = f2f["target"]
    lSelectedMotifs   = f2f["motifs"]
    reference_motif   = f2f["reference"]

    stats_rows = []

    fig = plt.figure(figsize=(4.5, 6))
    i = 0
    for sLabelRep, sClvName, dfClvConcat in lAllMappedClv:
        if sTargetSampleName != display_name(sClvName):
            continue
        i += 1
        sRepSampleName = f'{sLabelRep}_{display_name(sClvName)}'

        filtered_dfClvConcat = dfClvConcat[dfClvConcat['3nt_motif'].isin(lSelectedMotifs)]

        group_stats = filtered_dfClvConcat.groupby('3nt_motif')['clvOverDict'].agg(['mean', 'count']).reset_index()
        group_stats['label'] = group_stats.apply(
            lambda row: f"{row['3nt_motif']} (n={row['count']})", axis=1
        )
        group_stats_sorted = group_stats.sort_values('mean', ascending=False)
        ordered_motifs     = group_stats_sorted['3nt_motif'].tolist()
        label_map          = dict(zip(group_stats_sorted['3nt_motif'], group_stats_sorted['label']))

        ax = fig.add_subplot(1, 2, i)
        sns.boxplot(data=filtered_dfClvConcat, y='clvOverDict', x='3nt_motif',
                    order=ordered_motifs, color='#edc47e', width=0.5, showfliers=False, ax=ax)
        sns.despine(ax=ax)
        ax.set_xticks(range(len(ordered_motifs)))
        ax.set_xticklabels([label_map[m] for m in ordered_motifs], rotation=90)
        ax.set_title(sLabelRep)
        ax.set_ylabel('Relative cleavage efficiency')

        if sLabelRep == 'Rep1':
            ax.set_ylim(f2f["rep1_ylim"])
        elif sLabelRep == 'Rep2':
            ax.set_ylim(f2f["rep2_ylim"])

        groups = {motif: filtered_dfClvConcat[filtered_dfClvConcat['3nt_motif'] == motif]['clvOverDict']
                for motif in lSelectedMotifs}

        for motif in lSelectedMotifs:
            if motif == reference_motif:
                continue
            p = mannwhitneyu(groups[motif], groups[reference_motif],
                            alternative='two-sided').pvalue
            print(f'{sLabelRep} {motif} vs. {reference_motif}, '
                f'Mann-Whitney-U, p-value: {p:,}')
            stats_rows.append(dict(
                replicate   = sLabelRep,
                sample      = display_name(sClvName),
                comparison  = f'{motif}_vs_{reference_motif}',
                n_motif     = len(groups[motif]),
                n_reference = len(groups[reference_motif]),
                p_value     = p,
            ))

    ax.set_xlabel('')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(sOutputDir, 'main_fig', f2f["out"]), dpi=1200)
    plt.close()

    # export mann whitney test data

    stats_path = os.path.join(sOutputDir, 'main_fig', f2f["stats_out"])
    pd.DataFrame(stats_rows).to_csv(stats_path, sep='\t', index=False)
    print(f'wrote {stats_path}')
else:
    print("fig2f is null; skipping Fig-2f")
    
### Generate pLogo/kpLogo Foreground Input Sequences (Top 4000) ###
### Export Entire Mapped Table ###
### Write tables (pLogo input + full mapped table) ###

# sOutputTopSequences = os.path.join(logo_input_dir, 'Top_4000_Seq.xlsx')
sOutputMappedClv    = os.path.join(sOutputDir, 'mapped_clv_seq_full.xlsx')

with pd.ExcelWriter(sOutputMappedClv,    engine='xlsxwriter') as f2:
# pd.ExcelWriter(sOutputTopSequences, engine='xlsxwriter') as f1, \
    for sLabelRep, sClvName, dfClvConcat in lAllMappedClv:
        sample_key  = _sample_key(sClvName)                 # 'Lbu-CR' / 'Lbu-NP' / 'Dict-DT'
        
        if sample_key.endswith('-DT'):
            continue 
        
        sOutputName = f'{sLabelRep}-{display_name(sClvName)}'   # 'Rep1-Lbu-withProtein' (xlsx sheet name)

        # (1) per-sample .txt for kpLogo: top 4000 8-mers
        #     CR -> foreground for kpLogo; NP -> background for kpLogo; DT -> not used by kpLogo
        rep_num = sClvName.split('_', 1)[0]           # '1' or '2'
        out_txt = os.path.join(logo_input_dir, f'{rep_num}_{sample_key}.txt')
        dfClvConcat.iloc[:nTopNumber, ]['8nt_motif'].to_csv(
            out_txt, index=False, header=False)

        # # (2) Top 4000 8-mers per sheet
        # dfClvConcat.iloc[:nTopNumber, ]['8nt_motif'].to_excel(
        #     f1, sheet_name=sOutputName, index=False)

        # (3) Full mapped table per sheet
        dfClvConcat.to_excel(f2, sheet_name=sOutputName, index=False)

        print(sOutputName)