#!/usr/local/bin/snakemake
import os, time
import pandas as pd

# 
#
configfile: "config.yaml"

RAW_DATA_DIR = config.get("input_dir", "input_dir")
INTERMEDIATE_DIR = config.get("intermediate_dir", "intermediate_files")
PPSEQ_DIR = os.path.join(INTERMEDIATE_DIR, config.get("preprocessed_seq_dir", "preprocessed_sequences"))
OUTPUT_DIR = config["output_dir"]
os.makedirs(INTERMEDIATE_DIR, exist_ok=True)
os.makedirs(PPSEQ_DIR,        exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "main_fig"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "supp_fig"), exist_ok=True)

# meta = pd.read_csv(config.get("sample_sheet", "samples.csv")).set_index("sample")
# meta = pd.read_csv(config["sample_sheet"]).set_index(["replicate", "sample"])
# meta = pd.read_csv(config["sample_sheet"], dtype={"replicate": str}).set_index(["replicate", "sample"])
meta = pd.read_csv(config["sample_sheet"], dtype={"replicate": str}).set_index(["replicate", "protein", "type"])
PAIRS = list(meta.index)

# CLEAVAGE = meta.index[meta["class"] == "cleavage"].tolist()
# DICTIONARY = meta.index[meta["class"] == "dictionary"].tolist()
# SAMPLES = CLEAVAGE + DICTIONARY
# CLEAVAGE_PAIRS   = [(r, s) for (r, s) in PAIRS if meta.loc[(r, s), "class"] == "cleavage"]
# DICTIONARY_PAIRS = [(r, s) for (r, s) in PAIRS if meta.loc[(r, s), "class"] == "dictionary"]

MIN_QUAL = 30
MIN_QUALIFIED = 90

# write a log for each run
os.makedirs('logs', exist_ok=True)
DATE = time.strftime('%y%m%d')
LOG = 'logs/%s.log' % DATE
logf = open(LOG, 'a')
STAMP = '========== NEW RUN at %s ==========\n' % time.time()
logf.write(STAMP)
logf.close()

rule all:
    input:
        f'logs/{DATE}.log',
        expand(f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.collapsed.tsv',
               zip,
               replicate=[r for (r, _, _) in PAIRS],
               protein=[p for (_, p, _) in PAIRS],
               type=[t for (_, _, t) in PAIRS]),        
        # expand('qualitycheck/{sample}_L2_1.fastqc.zip', sample=SAMPLES),
        # expand(f'{PPSEQ_DIR}/{{sample}}.collapsed.tsv', sample=SAMPLES),
        # expand(f'{PPSEQ_DIR}/{{replicate}}/{{sample}}.collapsed.tsv', zip, replicate=[r for (r, _) in PAIRS], sample=[s for (_, s) in PAIRS])
    
        # plots + tables
        expand(f'{OUTPUT_DIR}/main_fig/fig2c_{{p}}_motif-clv-eff.pdf', p=['lbu','lwa','bt3','rfx']),
        f'{OUTPUT_DIR}/main_fig/fig2e_psp_motif-clv-eff.pdf',
        f'{OUTPUT_DIR}/main_fig/fig2f_psp_clv-eff.pdf',
        f'{OUTPUT_DIR}/main_fig/fig2_psp_mann-whitney.tsv',
        f'{OUTPUT_DIR}/supp_fig/supp3ab_dict_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp4a_rep1_clv_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp4b_rep2_clv_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp5a_rep1_clv_eff.png',
        f'{OUTPUT_DIR}/supp_fig/supp5b_rep2_clv_eff.png',
        f'{OUTPUT_DIR}/supp_fig/supp6a_psp_clv-eff.pdf',
        f'{INTERMEDIATE_DIR}/{config["logo_input_dir"]}/TopSeqForLogo_260129v1.xlsx',
        f'{OUTPUT_DIR}/mapped_clv_seq_full.xlsx',

#rule fastqc:
#    output: 'qualitycheck/{sample}_L2_1.fastqc.zip'
#    threads: 32
#    shell: 'fastqc -o qualitycheck -t {threads} {RAW_DATA_DIR}/{wildcards.sample}_L2_1.fq.gz'

rule trim_adaptor:
    input:  f'{RAW_DATA_DIR}/{{replicate}}_{{protein}}-{{type}}.fastq.gz'
    output: f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.fastq.gz'
    params:
        adapter_5p = lambda w: meta.loc[(w.replicate, w.protein, w.type), "adapter_5p"],
        adapter_3p = lambda w: meta.loc[(w.replicate, w.protein, w.type), "adapter_3p"],
        # adapter_5p = lambda w: meta.loc[(w.replicate, w.sample), "adapter_5p"],
        # adapter_3p = lambda w: meta.loc[(w.replicate, w.sample), "adapter_3p"],
    shell: 'cutadapt {input} --cores 30 -g {params.adapter_5p}...{params.adapter_3p} --trimmed-only 2>> {LOG} | gzip -c - > {output}'
    # shell: 'cutadapt {params.input_dir}/{wildcards.replicate}_{wildcards.protein}-{wildcards.type}_{READ}.fastq.gz \
            # --cores 30 -g {params.adapter_5p}...{params.adapter_3p} --trimmed-only 2>> {LOG} | \
			# gzip -c - > {output}'

rule trim_random_adaptor_sequences:
    input:  f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.fastq.gz'
    output: f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.trim3part.fastq.gz'
    params:
        random_at_3p = lambda w: meta.loc[(w.replicate, w.protein, w.type), "random_at_3p"],
        min_len      = lambda w: meta.loc[(w.replicate, w.protein, w.type), "min_len"],
        max_len      = lambda w: meta.loc[(w.replicate, w.protein, w.type), "max_len"],
    shell: 'cutadapt {input} -u {params.random_at_3p} -m {params.min_len} -M {params.max_len} --cores 30 2>> {LOG} | \
			gzip -c - > {output}'

rule quality_filter:
    input:  f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.trim3part.fastq.gz'
    output: f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.qc.fastq.gz'
    shell: 'zcat {input} | \
            fastq_quality_filter - -v -Q33 -q {MIN_QUAL} -p {MIN_QUALIFIED} 2>> {LOG} | \
            gzip -c - > {output}'

rule collapsing:
    input:  f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.trim.qc.fastq.gz'
    output: f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.collapsed.fa.gz'
    shell: 'zcat {input} | \
			fastx_collapser - -Q33 2>> {LOG} | \
			gzip -c - > {output}'

rule convert_to_tsv:
    input:  f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.collapsed.fa.gz'
    output: f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.collapsed.tsv'
    shell: 'zcat {input} | \
			seqkit fx2tab 2>> {LOG} > {output}'

# Call generate_plots.py
rule generate_plots:
    input:
        expand(f'{PPSEQ_DIR}/{{replicate}}_{{protein}}-{{type}}.collapsed.tsv',
               zip,
               replicate=[r for (r, _, _) in PAIRS],
               protein=[p for (_, p, _) in PAIRS],
               type=[t for (_, _, t) in PAIRS]),
        config["sample_sheet"],
    output:
        # main figures
        expand(f'{OUTPUT_DIR}/main_fig/fig2c_{{p}}_motif-clv-eff.pdf', p=['lbu','lwa','bt3','rfx']),
        f'{OUTPUT_DIR}/main_fig/fig2e_psp_motif-clv-eff.pdf',
        f'{OUTPUT_DIR}/main_fig/fig2f_psp_clv-eff.pdf',
        f'{OUTPUT_DIR}/main_fig/fig2_psp_mann-whitney.tsv',
        # supplementary figures
        f'{OUTPUT_DIR}/supp_fig/supp3ab_dict_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp4a_rep1_clv_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp4b_rep2_clv_depth.png',
        f'{OUTPUT_DIR}/supp_fig/supp5a_rep1_clv_eff.png',
        f'{OUTPUT_DIR}/supp_fig/supp5b_rep2_clv_eff.png',
        f'{OUTPUT_DIR}/supp_fig/supp6a_psp_clv-eff.pdf',
        # tables
        f'{INTERMEDIATE_DIR}/{config["logo_input_dir"]}/TopSeqForLogo_260129v1.xlsx',
        f'{OUTPUT_DIR}/mapped_clv_seq_full.xlsx',
    params:
        config_path = os.path.abspath("config.yaml"),
    threads: 1
    script: "generate_plots.py"