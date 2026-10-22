#!/usr/local/bin/snakemake
import os, time
import pandas as pd
import shutil, tempfile

configfile: "config.yaml"

RAW_DATA_DIR = config.get("input_dir", "input_dir")
INTERMEDIATE_DIR = config.get("intermediate_dir", "intermediate_files")
PPSEQ_DIR = os.path.join(INTERMEDIATE_DIR, config.get("preprocessed_seq_dir", "preprocessed_sequences"))
OUTPUT_DIR = config["output_dir"]
LOGO_OUTPUT_DIR = os.path.join(INTERMEDIATE_DIR, "logo_output")
os.makedirs(PPSEQ_DIR,        exist_ok=True)
os.makedirs(INTERMEDIATE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "main_fig"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "supp_fig"), exist_ok=True)
os.makedirs(os.path.join(INTERMEDIATE_DIR, "logo_input"), exist_ok=True)
os.makedirs(LOGO_OUTPUT_DIR, exist_ok=True)

KP_LOGO_ABS = shutil.which("kpLogo")
if not KP_LOGO_ABS:
    raise Exception("kpLogo not found. Please install according to the Installation Guide in README.")
print(f"Using kpLogo at: {KP_LOGO_ABS}")

# meta = pd.read_csv(config.get("sample_sheet", "samples.csv")).set_index("sample")
# meta = pd.read_csv(config["sample_sheet"]).set_index(["replicate", "sample"])
# meta = pd.read_csv(config["sample_sheet"], dtype={"replicate": str}).set_index(["replicate", "sample"])
meta = pd.read_csv(config["sample_sheet"], dtype={"replicate": str}).set_index(["replicate", "protein", "type"])
PAIRS = list(meta.index)
KPLOGO_PAIRS = [(r, p) for (r, p, t) in PAIRS if t == "CR"]
# print(f"kpLogo pairs ({len(KPLOGO_PAIRS)}): {KPLOGO_PAIRS}")

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
        # f'{INTERMEDIATE_DIR}/{config["logo_input_dir"]}/Top_4000_Seq.xlsx',
        f'{OUTPUT_DIR}/mapped_clv_seq_full.xlsx',

        #kpLogo
        expand(os.path.join(LOGO_OUTPUT_DIR,
                            "{replicate}_{protein}/image/{replicate}_{protein}.info.eps"),
               zip,
               replicate=[r for (r, p) in KPLOGO_PAIRS],
               protein=[p for (r, p) in KPLOGO_PAIRS]),

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
        # f'{INTERMEDIATE_DIR}/{config["logo_input_dir"]}/Top_4000_Seq.xlsx',
        f'{OUTPUT_DIR}/mapped_clv_seq_full.xlsx',
       # kpLogo input .txt files (top-4000 CR + all-8mer NP per sample)
        expand(f'{INTERMEDIATE_DIR}/logo_input/{{replicate}}_{{protein}}-{{type}}.txt',
               zip,
               replicate=[r for (r, _, t) in PAIRS if t != "DT"],
               protein=[p for (_, p, t) in PAIRS if t != "DT"],
               type=[t for (_, _, t) in PAIRS if t != "DT"]),
    params:
        config_path = os.path.abspath("config.yaml"),
    threads: 1
    script: "generate_plots.py"

rule run_kpLogo:
    input:
        fg = f'{INTERMEDIATE_DIR}/logo_input/{{replicate}}_{{protein}}-CR.txt',
        bg = f'{INTERMEDIATE_DIR}/logo_input/{{replicate}}_{{protein}}-NP.txt',
    output:
        # images
        freq_eps = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.freq.eps"),
        freq_pdf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.freq.pdf"),
        freq_png = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.freq.png"),
        info_eps = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.info.eps"),
        info_pdf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.info.pdf"),
        info_png = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.info.png"),
        kplogo_eps = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.kplogo.eps"),
        kplogo_pdf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.kplogo.pdf"),
        kplogo_png = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.kplogo.png"),
        plogo_eps = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.plogo.eps"),
        plogo_pdf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.plogo.pdf"),
        plogo_png = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.plogo.png"),
        all_pdf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/image/{replicate}_{protein}.all.pdf"),
        # data
        pwm_txt = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/data/{replicate}_{protein}.pwm.txt"),
        pass_cutoff_txt = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/data/{replicate}_{protein}.pass.p.cutoff.txt"),
        pass_cutoff_inf = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/data/{replicate}_{protein}.pass.p.cutoff.txt.inf"),
        most_sig_txt = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/data/{replicate}_{protein}.most.significant.each.position.txt"),
        bonferroni = os.path.join(LOGO_OUTPUT_DIR, "{replicate}_{protein}/data/{replicate}_{protein}.Bonferroni.significant.frequency.txt"),
    params:
        min_count  = config.get("min_count", 5),
        p_cutoff   = config.get("p_cutoff", 0.05),
        fdr        = config.get("fdr", False),
        markov     = config.get("markov", None),
        region     = config.get("region", None),
        weighted   = config.get("weighted", True),
        ranked     = config.get("ranked", False),
        max_k      = config.get("max_k", 4),
        extra      = config.get("extra_args", ""),
    threads: 1
    resources:
        kpLogo_lock = 1
    run:
        rep, prot = wildcards.replicate, wildcards.protein
        sample_tag = f"{rep}_{prot}"

        # 1. unique per-job temp dir
        job_root = os.path.join(tempfile.gettempdir(),
                                f"kpLogo_{sample_tag}_{os.getpid()}_{int(time.time())}")
        os.makedirs(job_root, exist_ok=True)
        for var in ("TMPDIR", "TEMP", "TMP", "R_TMPDIR", "GS_TMPDIR"):
            os.environ[var] = job_root
        print(f"Job temporary directory: {job_root}")

        # 2. destination dirs
        image_dir = os.path.join(LOGO_OUTPUT_DIR, sample_tag, "image")
        data_dir  = os.path.join(LOGO_OUTPUT_DIR, sample_tag, "data")
        os.makedirs(image_dir, exist_ok=True)
        os.makedirs(data_dir,  exist_ok=True)

        workflow_dir = os.getcwd()
        fg_abs = os.path.join(workflow_dir, input.fg)
        bg_abs = os.path.join(workflow_dir, input.bg)

        work_dir = os.path.join(job_root, "work")
        os.makedirs(work_dir, exist_ok=True)

        try:
            fg_name = os.path.basename(fg_abs)
            bg_name = os.path.basename(bg_abs)
            shutil.copy(fg_abs, os.path.join(work_dir, fg_name))
            shutil.copy(bg_abs, os.path.join(work_dir, bg_name))

            cmd = [KP_LOGO_ABS, os.path.join(work_dir, fg_name)]

            if params.weighted:
                cmd.append("-weighted")
            elif params.ranked:
                cmd.append("-ranked")

            cmd.extend(["-max_k", str(params.max_k)])
            cmd.extend(["-minCount", str(params.min_count)])

            if params.fdr:
                cmd.append("-FDR")
            else:
                cmd.extend(["-pc", str(params.p_cutoff)])

            if params.markov is not None:
                cmd.extend(["-markov", str(params.markov)])

            cmd.extend(["-bgfile", os.path.join(work_dir, bg_name)])

            if params.region:
                cmd.extend(["-region", params.region])

            cmd.extend(["-o", os.path.join(work_dir, "kpLogo")])

            if params.extra:
                cmd.extend(params.extra.split())

            print(f"Processing {sample_tag}")
            print("Running: " + " ".join(cmd))
            shell(" ".join(cmd))

            # 3. move kpLogo outputs to desired locations
            source_patterns = {
                "kpLogo.freq.eps": "freq_eps",
                "kpLogo.freq.pdf": "freq_pdf",
                "kpLogo.freq.png": "freq_png",
                "kpLogo.info.eps": "info_eps",
                "kpLogo.info.pdf": "info_pdf",
                "kpLogo.info.png": "info_png",
                "kpLogo.most.significant.each.position.eps": "kplogo_eps",
                "kpLogo.most.significant.each.position.pdf": "kplogo_pdf",
                "kpLogo.most.significant.each.position.png": "kplogo_png",
                "kpLogo.eps":     "plogo_eps",
                "kpLogo.pdf":     "plogo_pdf",
                "kpLogo.png":     "plogo_png",
                "kpLogo.all.pdf": "all_pdf",
                "kpLogo.pwm.txt": "pwm_txt",
                "kpLogo.pass.p.cutoff.txt":     "pass_cutoff_txt",
                "kpLogo.pass.p.cutoff.txt.inf": "pass_cutoff_inf",
                "kpLogo.most.significant.each.position.txt": "most_sig_txt",
                "kpLogo.Bonferroni.significant.frequency.txt": "bonferroni",
            }

            for src_name, out_attr in source_patterns.items():
                src_path = os.path.join(work_dir, src_name)
                dst_path = getattr(output, out_attr)
                if os.path.exists(src_path):
                    shutil.move(src_path, dst_path)
                    print(f"Moved: {src_name} -> {dst_path}")
                else:
                    print(f"Warning: {src_name} not found in {work_dir}")

            # save if additional barplot is generated
            barplot_src = os.path.join(work_dir, "kpLogo.barplot.pdf")
            if os.path.exists(barplot_src):
                shutil.move(barplot_src, os.path.join(image_dir, f"{sample_tag}.barplot.pdf"))
                print(f"Moved barplot -> {image_dir}")

            # touch all outputs so timestamps reflect this run
            for out_attr in source_patterns.values():
                p = getattr(output, out_attr)
                if os.path.exists(p):
                    os.utime(p, None)

            # verify that all expected output files exist
            missing = [getattr(output, a) for a in source_patterns.values()
                       if not os.path.exists(getattr(output, a))]
            if missing:
                print("WARNING: missing outputs after move:")
                for m in missing:
                    print(f"  {m}")

        finally:
            # cleanup
            shutil.rmtree(job_root, ignore_errors=True)
            print(f"Cleaned up: {job_root}")