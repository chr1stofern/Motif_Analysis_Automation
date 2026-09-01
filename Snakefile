# Snakefile

import os
import glob
import shutil
import tempfile
import subprocess
import fcntl
import shlex
import time
# ============================================
# 0. USAGE INSTRUCTIONS
# ============================================

# # Run the pipeline
# snakemake -j 4

# # Clean up logo images related outputs
# snakemake -j 4 clean_logo

# # Clean up tex related outputs
# snakemake -j 4 clean_tex

# # Clean everything (remove all outputs)
# snakemake -j 4 clean_all

# ============================================
# 1. CONFIGURATION
# ============================================
configfile: "config.yaml"

# ============================================
# 2. DEFINE PATHS AND FUNCTIONS
# ============================================
INPUT_DIR = config.get("input_dir", "fg")
OUTPUT_DIR = config.get("output_dir", "kplogo_output")
INPUT_SUFFIX = config.get("input_suffix", "_weighted_plogo_input_v2.txt")
KP_LOGO = config.get("kpLogo_path", "kpLogo")
BATCH_NAME = INPUT_DIR.split("_")[1]
BATCH_LOGO_DIR = os.path.join(OUTPUT_DIR, "logo", BATCH_NAME)

def get_samples():
    pattern = os.path.join(INPUT_DIR, "*" + INPUT_SUFFIX)
    files = glob.glob(pattern)
    samples = []
    for f in files:
        basename = os.path.basename(f)
        sample = basename.replace(INPUT_SUFFIX, "")
        samples.append(sample)
    
    # Sort by numeric prefix (e.g., "1", "2", "3", ..., "10", "11", ...)
    def extract_number(s):
        # Extract the leading number before the first dash
        parts = s.split("-")
        try:
            return int(parts[0])
        except ValueError:
            return 999  # Put samples without numbers at the end
    
    return sorted(samples, key=extract_number)

def get_input_file(wildcards):
    return os.path.join(INPUT_DIR, wildcards.sample + INPUT_SUFFIX)

def get_bg_file_path():
    bg_file = config.get("bg_file", None)
    if bg_file is None:
        return None
    if os.path.isabs(bg_file):
        return bg_file
    return os.path.abspath(bg_file)

BG_FILE_ABS = get_bg_file_path()
if BG_FILE_ABS:
    print("Using background file: " + BG_FILE_ABS)

def get_kpLogo_path():
    if os.path.isabs(KP_LOGO):
        return KP_LOGO
    if os.path.exists(KP_LOGO):
        return os.path.abspath(KP_LOGO)
    import shutil
    which_kpLogo = shutil.which(KP_LOGO)
    if which_kpLogo:
        return which_kpLogo
    raise Exception("kpLogo not found. Please install with: conda install -c bioconda kplogo")

SAMPLES = get_samples()
print("Found " + str(len(SAMPLES)) + " samples: " + str(SAMPLES))

KP_LOGO_ABS = get_kpLogo_path()
print("Using kpLogo at: " + KP_LOGO_ABS)

# LaTeX parameters (read from config or use defaults)
COLS = config.get("latex_cols", 7)
FIGURE_HEIGHT = config.get("latex_figure_height", "2cm")
FONT_SIZE = config.get("latex_font_size", "\\scriptsize")
PAGE_MARGIN = config.get("latex_page_margin", "0.4cm")
LOGO_TYPES = config.get("latex_logo_types", ["plogo", "kplogo"])
PAGE_TITLES = config.get("latex_page_titles", {
    "plogo": "Probability Logo (pLogo)",
    "kplogo": "K-Mers Probability Logo (kpLogo)"
})
MAX_K = config.get("max_k", 4)
LATEX_OUTPUT = config.get("latex_output", "Combine")

# ============================================
# 3. RULE: all
# ============================================
rule all:
    input:
        # Image files (now under logo/batch/sample/image/)
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.eps"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.pdf"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.png"), sample=SAMPLES),
        # ... similarly for all .info, .kplogo, .plogo, .all.pdf ...
        # Data files
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pwm.txt"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt.inf"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.most.significant.each.position.txt"), sample=SAMPLES),
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.Bonferroni.significant.frequency.txt"), sample=SAMPLES),
        # Combined PDF from LaTeX (still uses batch‑specific logo images)
        LATEX_OUTPUT + "_" + str(MAX_K) + "k.pdf",
        # Efficiency (unchanged)
        expand("eff/eff_{batch}.csv", batch=BATCH_NAME),
        # Heatmaps (unchanged)
        os.path.join(OUTPUT_DIR, "heatmap", "_all_heatmaps_generated.txt"),
        OUTPUT_DIR + "/heatmap/" + BATCH_NAME + "/" + config.get("combined_output", "combined_heatmap") + ".pdf",
        # Preference CSV (batch‑specific)
        f"preference_data_{BATCH_NAME}.csv"
# ============================================
# 4. RULE: run_kpLogo
# ============================================
rule run_kpLogo:
    input:
        seq = INPUT_DIR + "/{sample}" + INPUT_SUFFIX
    output:
        # Image files
        freq_eps = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.eps"),
        freq_pdf = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.pdf"),
        freq_png = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.freq.png"),
        info_eps = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.info.eps"),
        info_pdf = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.info.pdf"),
        info_png = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.info.png"),
        kplogo_eps = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.kplogo.eps"),
        kplogo_pdf = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.kplogo.pdf"),
        kplogo_png = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.kplogo.png"),
        plogo_eps = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.plogo.eps"),
        plogo_pdf = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.plogo.pdf"),
        plogo_png = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.plogo.png"),
        all_pdf = os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.all.pdf"),
        # Data files
        pwm_txt = os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pwm.txt"),
        pass_cutoff_txt = os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt"),
        pass_cutoff_inf = os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt.inf"),
        most_sig_txt = os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.most.significant.each.position.txt"),
        bonferroni = os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.Bonferroni.significant.frequency.txt")
    params:
        min_count = config.get("min_count", 5),
        p_cutoff = config.get("p_cutoff", 0.05),
        fdr = config.get("fdr", False),
        markov = config.get("markov", None),
        region = config.get("region", None),
        weighted = config.get("weighted", True),
        ranked = config.get("ranked", False),
        bg_file = BG_FILE_ABS,
        extra = config.get("extra_args", "")
    log:
        "logs/{sample}.log"
    threads:
        1
    resources:
        kpLogo_lock = 1
    run:
        sample = wildcards.sample

        # === 1. Create a unique job directory ===
        job_id = f"{sample}_{os.getpid()}_{int(time.time())}"
        job_root = os.path.join(tempfile.gettempdir(), f"kpLogo_job_{job_id}")
        os.makedirs(job_root, exist_ok=True)
        
        os.environ['TMPDIR'] = job_root
        os.environ['R_SESSION_TMPDIR'] = job_root
        os.environ['R_USER'] = job_root

        # === 2. Set ALL temporary directory environment variables ===
        os.environ['TMPDIR'] = job_root
        os.environ['TEMP'] = job_root
        os.environ['TMP'] = job_root
        os.environ['R_TMPDIR'] = job_root
        os.environ['GS_TMPDIR'] = job_root
        
        print(f"Job temporary directory: {job_root}")
        
        # Ghostscript options
        os.environ['GS_OPTIONS'] = '-dAutoRotatePages=/None'
        
        # Create destination directories inside BATCH_LOGO_DIR
        image_dir = os.path.join(BATCH_LOGO_DIR, sample, "image")
        data_dir = os.path.join(BATCH_LOGO_DIR, sample, "data")
        os.makedirs(image_dir, exist_ok=True)
        os.makedirs(data_dir, exist_ok=True)

        workflow_dir = os.getcwd()
        input_abs = os.path.join(workflow_dir, input.seq)
        input_name = os.path.basename(input_abs)

        # kpLogo's working directory
        work_dir = os.path.join(job_root, "work")
        os.makedirs(work_dir, exist_ok=True)
        
        try:
            shutil.copy(input_abs, os.path.join(work_dir, input_name))

            bg_name = None
            if params.bg_file and os.path.exists(os.path.join(workflow_dir, params.bg_file)):
                bg_name = os.path.basename(params.bg_file)
                shutil.copy2(os.path.join(workflow_dir, params.bg_file), os.path.join(work_dir, bg_name))

            cmd = [KP_LOGO_ABS, os.path.join(work_dir, input_name)]
            
            if params.weighted:
                cmd.append("-weighted")
            elif params.ranked:
                cmd.append("-ranked")

            cmd.extend(["-max_k", str(MAX_K)])
            cmd.extend(["-minCount", str(params.min_count)])

            if params.fdr:
                cmd.append("-FDR")
            else:
                cmd.extend(["-pc", str(params.p_cutoff)])

            if params.markov is not None:
                cmd.extend(["-markov", str(params.markov)])

            if bg_name:
                cmd.extend(["-bgfile", os.path.join(work_dir, bg_name)])

            if params.region:
                cmd.extend(["-region", params.region])

            cmd.extend(["-o", os.path.join(work_dir, "kpLogo")])

            if params.extra:
                cmd.extend(params.extra.split())

            print("Processing " + sample)
            print("Running: " + " ".join(cmd))
            print(f"Working in: {work_dir}")
            
            # Run kpLogo
            shell(" ".join(cmd), log=log[0])

            # ============================================================
            # Robust file moving: list all files in work_dir and move
            # ============================================================
            # Define a mapping from kpLogo source patterns to destination keys
            # The keys must match the output dictionary keys.
            source_patterns = {
                "kpLogo.freq.eps": "freq.eps",
                "kpLogo.freq.pdf": "freq.pdf",
                "kpLogo.freq.png": "freq.png",
                "kpLogo.info.eps": "info.eps",
                "kpLogo.info.pdf": "info.pdf",
                "kpLogo.info.png": "info.png",
                "kpLogo.most.significant.each.position.eps": "kplogo.eps",
                "kpLogo.most.significant.each.position.pdf": "kplogo.pdf",
                "kpLogo.most.significant.each.position.png": "kplogo.png",
                "kpLogo.eps": "plogo.eps",
                "kpLogo.pdf": "plogo.pdf",
                "kpLogo.png": "plogo.png",
                "kpLogo.all.pdf": "all.pdf",
                "kpLogo.pwm.txt": "pwm.txt",
                "kpLogo.pass.p.cutoff.txt": "pass.p.cutoff.txt",
                "kpLogo.pass.p.cutoff.txt.inf": "pass.p.cutoff.txt.inf",
                "kpLogo.most.significant.each.position.txt": "most.significant.each.position.txt",
                "kpLogo.Bonferroni.significant.frequency.txt": "Bonferroni.significant.frequency.txt",
            }

            # Define the destination dictionary (using the output objects)
            outputs = {
                "freq.eps": output.freq_eps,
                "freq.pdf": output.freq_pdf,
                "freq.png": output.freq_png,
                "info.eps": output.info_eps,
                "info.pdf": output.info_pdf,
                "info.png": output.info_png,
                "kplogo.eps": output.kplogo_eps,
                "kplogo.pdf": output.kplogo_pdf,
                "kplogo.png": output.kplogo_png,
                "plogo.eps": output.plogo_eps,
                "plogo.pdf": output.plogo_pdf,
                "plogo.png": output.plogo_png,
                "all.pdf": output.all_pdf,
                "pwm.txt": output.pwm_txt,
                "pass.p.cutoff.txt": output.pass_cutoff_txt,
                "pass.p.cutoff.txt.inf": output.pass_cutoff_inf,
                "most.significant.each.position.txt": output.most_sig_txt,
                "Bonferroni.significant.frequency.txt": output.bonferroni,
            }

            # Move each expected file
            for src_name, key in source_patterns.items():
                src_path = os.path.join(work_dir, src_name)
                dst_path = outputs[key]
                if os.path.exists(src_path):
                    shutil.move(src_path, dst_path)
                    print(f"Moved: {src_name} -> {dst_path}")
                else:
                    print(f"Warning: {src_name} not found in {work_dir}. Skipping.")

            # Also move any barplot file if present
            barplot_src = os.path.join(work_dir, "kpLogo.barplot.pdf")
            barplot_dst = os.path.join(image_dir, sample + ".barplot.pdf")
            if os.path.exists(barplot_src):
                shutil.move(barplot_src, barplot_dst)
                print(f"Moved barplot: {barplot_src} -> {barplot_dst}")

            # Touch all output files to ensure timestamps are updated
            for dst_path in outputs.values():
                if os.path.exists(dst_path):
                    os.utime(dst_path, None)

            # After moving, verify that all expected output files exist
            missing = [dst for dst in outputs.values() if not os.path.exists(dst)]
            if missing:
                print(f"WARNING: The following output files are still missing after moving:")
                for m in missing:
                    print(f"  {m}")
                # Optionally, list contents of work_dir for debugging
                print(f"Contents of {work_dir}:")
                for f in os.listdir(work_dir):
                    print(f"  {f}")

        finally:
            # Clean up the entire job directory
            shutil.rmtree(job_root, ignore_errors=True)
            print(f"Cleaned up: {job_root}")# ============================================
# 5. RULE: generate_latex
# ============================================
rule generate_latex:
    input:
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.{logo}.pdf"),
               sample=SAMPLES, logo=LOGO_TYPES)
    output:
        tex = LATEX_OUTPUT + "_" + str(MAX_K) + "k.tex"
    params:
        cols = COLS,
        figure_height = FIGURE_HEIGHT,
        font_size = FONT_SIZE,
        page_margin = PAGE_MARGIN,
        logo_types = LOGO_TYPES,
        page_titles = PAGE_TITLES
    run:
        samples = SAMPLES
        cols = params.cols
        rows = (len(samples) + cols - 1) // cols
        
        latex_lines = [
            "\\documentclass[a4paper,10pt]{article}",
            "\\usepackage[utf8]{inputenc}",
            "\\usepackage{adjustbox}",
            "\\usepackage{graphicx}",
            "\\usepackage{geometry}",
            "\\usepackage{caption}",
            "\\usepackage{subcaption}",
            "\\usepackage{float}",
            "\\usepackage{amsmath}",
            "",
            "\\DeclareGraphicsExtensions{.pdf,.PDF}",
            "",
            "\\geometry{",
            "    a4paper,",
            "    landscape,",
            "    top=" + params.page_margin + ",",
            "    bottom=" + params.page_margin + ",",
            "    left=" + params.page_margin + ",",
            "    right=" + params.page_margin + ",",
            "    headheight=0cm,",
            "    headsep=0cm,",
            "    footskip=0.3cm,",
            "    includeheadfoot",
            "}",
            "",
            "\\pagestyle{empty}",
            "\\captionsetup{font={small,sf}, labelfont={bf}, skip=2pt, justification=centering}",
            "",
            "\\begin{document}",
        ]
        
        for logo_type in params.logo_types:
            title = params.page_titles.get(logo_type, logo_type.upper() + " Logo")
            if logo_type == "kplogo":
                title = title + ": k=" + str(MAX_K)

            latex_lines.append("")
            latex_lines.append("% ============================================")
            latex_lines.append("% PAGE: " + title)
            latex_lines.append("% ============================================")
            latex_lines.append("\\begin{figure}[p]")
            latex_lines.append("    \\centering")
            latex_lines.append("    \\Large\\textbf{" + title + "}")
            latex_lines.append("    \\vspace{0.8cm}")
            latex_lines.append("")
            latex_lines.append("    \\begin{minipage}{0.97\\textwidth}")
            # latex_lines.append("    \\hspace{0.05\\textwidth}")
            latex_lines.append("        \\centering")

            for row_idx in range(rows):
                if row_idx > 0:
                    latex_lines.append("        \\vspace{0.4cm}")
                
                row_samples = samples[row_idx*cols:(row_idx+1)*cols]
                row_lines = []
                for sample in row_samples:
                    pdf_path = os.path.join(BATCH_LOGO_DIR, sample, "image", sample + "." + logo_type + ".pdf")
                    minipage = (
                        "    \\begin{minipage}[t][0.09\\textwidth][t]{0.12\\textwidth}\n"   # [t] = top-align
                        "        \\centering\n"
                        "        " + params.font_size + " \\textrm{" + sample + "}\n"
                        "        \\vspace{0.03cm}\n"
                        # "        \\adjustbox{width=\\textwidth, height=0.85\\textwidth, keepaspectratio}{\\includegraphics[width=\\textwidth]{" + pdf_path + "}}\n"
                        "        \\includegraphics[width=\\textwidth, keepaspectratio, trim=1.5cm 0cm 3cm 0cm, clip=true]{" + pdf_path + "}%\n"
                        # "        \\vspace{0.03cm}\n"
                        "    \\end{minipage}%\n"
                    )

                    row_lines.append(minipage)
                
                # === CHANGE IS HERE ===
                # Replace \hfill with \hspace{0.1cm} (or any value you want)
                latex_lines.append("    \\hspace*{0.3cm}\n".join(row_lines))
                # latex_lines.append("    \\hfill\n".join(row_lines))
                # latex_lines.append("    \\hspace{0.05\\textwidth}")
                latex_lines.append("    \\vspace{1.5cm}\n")
            
            latex_lines.append("    \\end{minipage}")
            latex_lines.append("\\end{figure}")
            latex_lines.append("")
            latex_lines.append("\\clearpage")
        
        latex_lines.append("")
        latex_lines.append("\\end{document}")
        
        with open(output.tex, "w") as f:
            f.write("\n".join(latex_lines))
        
        print("LaTeX file generated: " + output.tex)
# ============================================
# 6. RULE: compile_latex
# ============================================
rule compile_latex:
    input:
        tex = LATEX_OUTPUT + "_" + str(MAX_K) + "k.tex",
        logos = expand(OUTPUT_DIR + "/{sample}/image/{sample}.{logo}.pdf",
                       sample=SAMPLES, logo=LOGO_TYPES)
    params:
        aux_files = [
            LATEX_OUTPUT + ".aux",
            LATEX_OUTPUT + ".log",
            LATEX_OUTPUT + ".out"
        ]
    output:
        pdf = LATEX_OUTPUT + "_" + str(MAX_K) + "k.pdf"

    threads: 1
    run:
        # Run pdflatex twice to resolve cross-references
        for i in range(2):
            cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", input.tex]
            print("Running: " + " ".join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True)

            # Check for errors
            if result.returncode != 0:
                print("ERROR: pdflatex failed with return code " + str(result.returncode))
                print("STDOUT:")
                print(result.stdout)
                print("STDERR:")
                print(result.stderr)
                raise RuntimeError("pdflatex compilation failed")

            # Also print any warnings from the log
            if result.stderr:
                print("PDFLATEX output:")
                print(result.stderr)

        # Verify output exists
        if not os.path.exists(output.pdf):
            raise RuntimeError("PDF output not generated: " + output.pdf)

        # Clean up auxiliary files
        for aux_file in params.aux_files:
            if os.path.exists(aux_file):
                os.remove(aux_file)
                print("Removed: " + aux_file)

        print("PDF generated: " + output.pdf)

# ============================================
# 7 RULE: generate_efficiency_files
# ============================================
rule generate_efficiency_files:
    input:
        expand("fluor/fluor_{batch}.csv", batch=BATCH_NAME)
    output:
        expand("eff/eff_{batch}.csv", batch=BATCH_NAME)
    run:
        import glob
        
        # Find all fluor files
        fluor_files = glob.glob("fluor/fluor_*.csv")
        
        if not fluor_files:
            print("Warning: No fluorescence files found in fluor/ directory", file=sys.stderr)
            return
        
        # Run test-eff.py on all files
        cmd = ["python", "test-eff.py"]
        print("Running: " + " ".join(cmd))
        shell(" ".join(cmd))
        
        # Verify outputs were created
        for batch in ["d129", "n1060", "t136"]:
            expected = f"eff/eff_{batch}.csv"
            if not os.path.exists(expected):
                print(f"Warning: {expected} was not generated", file=sys.stderr)
# ============================================
# 7.5. RULE: generate_preference_csv
# ============================================
rule generate_preference_csv:
    input:
        pass_files = expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt"), sample=SAMPLES)
    output:
        csv = f"preference_data_{BATCH_NAME}.csv"
    run:
        import pandas as pd
        import os
        
        all_data = []
        base_order = ['A', 'C', 'G', 'T']
        
        for pass_file in input.pass_files:
            sample = os.path.basename(pass_file).replace(".pass.p.cutoff.txt", "")
            try:
                df = pd.read_csv(pass_file, sep='\t', comment='#', header=None)
                # Determine column structure
                try:
                    df.columns = ['kmer', 'pos', 'shift', 'statistics', 'pvalue', 
                                  'corrected_p', 'frac_obs', 'frac_exp', 'obs_exp', 'local_r']
                except ValueError:
                    df.columns = ['kmer', 'pos', 'shift', 'statistics', 'pvalue', 
                                  'corrected_p', 'n1', 'mean1', 'std1', 'n2', 'mean2', 'std2']
                # Keep only k=1 (monomers) and shift == 0
                kmer_df = df[df['kmer'].isin(base_order)]
                if 'shift' in kmer_df.columns:
                    kmer_df = kmer_df[kmer_df['shift'] == 0]
                # Drop duplicate rows (just in case)
                kmer_df = kmer_df.drop_duplicates(subset=['kmer', 'pos'])
                for _, row in kmer_df.iterrows():
                    all_data.append({
                        'sample': sample,
                        'pos': row['pos'],
                        'base': row['kmer'],
                        'value': row['statistics']   # test statistic
                    })
            except Exception as e:
                print(f"Warning: Could not parse {pass_file}: {e}", file=sys.stderr)
        
        if not all_data:
            raise RuntimeError("No data extracted from pass.cutoff files")
        
        df_all = pd.DataFrame(all_data)
        # Ensure uniqueness (group by sample, pos, base)
        df_all = df_all.groupby(['sample', 'pos', 'base'], as_index=False)['value'].first()
        df_all = df_all.sort_values(['sample', 'pos', 'base'])
        df_all.to_csv(output.csv, index=False)
        print(f"Preference CSV saved to: {output.csv}")
# ============================================
# 8. RULE: generate_heatmap
# ============================================
rule generate_heatmaps:
    input:
        pass_files = expand(os.path.join(BATCH_LOGO_DIR, "{sample}/data/{sample}.pass.p.cutoff.txt"), sample=SAMPLES),
        pref_csv = f"preference_data_{BATCH_NAME}.csv"
    output:
        done = os.path.join(OUTPUT_DIR, "heatmap", "_all_heatmaps_generated.txt")
    params:
        use_stats = config.get("heatmap_use_stats", True),
        cmap = config.get("heatmap_cmap", "inferno"),
        vmin = config.get("heatmap_vmin", None),
        vmax = config.get("heatmap_vmax", None),
        eff_vmin = config.get("eff_vmin", None),
        eff_vmax = config.get("eff_vmax", None),        
        figsize = config.get("heatmap_figsize", [6, 8]),
        fontsize = config.get("heatmap_fontsize", 10),
        annot = config.get("heatmap_annot", True),
        fmt = config.get("heatmap_fmt", ".2f"),
        square = config.get("heatmap_square", True),
        eff_dir = config.get("eff_dir", "eff"),
        cell_size = config.get("heatmap_cell_size", 0.4)
    run:
        import pandas as pd
        import os

        # Read positions from CSV
        df_pref = pd.read_csv(input.pref_csv)
        positions = sorted(df_pref['pos'].unique())
        
        heatmap_dir = os.path.join(OUTPUT_DIR, "heatmap")
        os.makedirs(heatmap_dir, exist_ok=True)
        
        batch_dir = os.path.join(heatmap_dir, BATCH_NAME)
        os.makedirs(batch_dir, exist_ok=True)
        
        # For each position, run heatmap.py with --pref-csv
        for pos in positions:
            output_pdf = os.path.join(batch_dir, f"Position-{pos}.pdf")
            cmd = [
                "python", "heatmap.py",
                "--pref-csv", input.pref_csv,
                output_pdf,
                "-p", str(pos),
                "--batch", BATCH_NAME,
                "--eff-dir", params.eff_dir,
                "--cell-size", str(params.cell_size)
            ]
            # Add global scales if set
            if params.vmin is not None:
                cmd.extend(["--global-vmin", str(params.vmin)])
            if params.vmax is not None:
                cmd.extend(["--global-vmax", str(params.vmax)])
            if params.eff_vmin is not None:
                cmd.extend(["--eff-global-vmin", str(params.eff_vmin)])
            if params.eff_vmax is not None:
                cmd.extend(["--eff-global-vmax", str(params.eff_vmax)])            
            if not params.use_stats:
                cmd.append("--no-stats")
            if params.cmap != "inferno":
                cmd.extend(["--cmap", params.cmap])
            if params.figsize != [6, 8]:
                cmd.extend(["--figsize", str(params.figsize[0]), str(params.figsize[1])])
            if params.fontsize != 10:
                cmd.extend(["--fontsize", str(params.fontsize)])
            if not params.annot:
                cmd.append("--no-annot")
            if params.fmt != ".2f":
                cmd.extend(["--fmt", params.fmt])
            if not params.square:
                cmd.append("--no-square")
            # Add skip samples
            skip_samples = config.get("skip_samples", [])
            if skip_samples:
                cmd.extend(["--skip-samples"] + skip_samples)
            
            cmd_str = " ".join(cmd)
            log_file = os.path.join("logs", f"heatmap_pos_{pos}.log")
            print(f"Generating heatmap for position {pos}: {output_pdf}")
            shell(cmd_str, log=log_file)
        
        # Create flag file
        with open(output.done, "w") as f:
            f.write(f"Generated heatmaps for positions: {positions}\n")
# ============================================
# 8.5 RULE: combine_heatmaps
# ============================================
rule combine_heatmaps:
    input:
        pref_csv = f"preference_data_{BATCH_NAME}.csv",
        eff_csv = expand("eff/eff_{batch}.csv", batch=BATCH_NAME)
    output:
        combined = OUTPUT_DIR + "/heatmap/{batch}/" + config.get("combined_output", "combined_heatmap") + ".pdf"
    params:
        positions = config.get("combine_positions", [3, 4, 5]),
        batch = BATCH_NAME,
        use_stats = config.get("heatmap_use_stats", True),
        cmap = config.get("heatmap_cmap", "inferno"),
        vmin = config.get("heatmap_vmin", None),
        vmax = config.get("heatmap_vmax", None),
        eff_vmin = config.get("eff_vmin", None),
        eff_vmax = config.get("eff_vmax", None),  
        figsize = config.get("heatmap_figsize", [6, 8]),
        fontsize = config.get("heatmap_fontsize", 10),
        annot = config.get("heatmap_annot", True),
        fmt = config.get("heatmap_fmt", ".2f"),
        square = config.get("heatmap_square", True),
        eff_dir = config.get("eff_dir", "eff"),
        cell_size = config.get("heatmap_cell_size", 0.4)
    run:
        cmd = [
            "python", "heatmap.py",
            "--pref-csv", input.pref_csv,
            output.combined,
            "--combine"
        ] + [str(p) for p in params.positions] + [
            "--batch", params.batch,
            "--eff-dir", params.eff_dir,
            "--cell-size", str(params.cell_size)
        ]
        # Add global scales
        if params.vmin is not None:
            cmd.extend(["--global-vmin", str(params.vmin)])
        if params.vmax is not None:
            cmd.extend(["--global-vmax", str(params.vmax)])
        if params.eff_vmin is not None:
            cmd.extend(["--eff-global-vmin", str(params.eff_vmin)])
        if params.eff_vmax is not None:
            cmd.extend(["--eff-global-vmax", str(params.eff_vmax)])     
        if not params.use_stats:
            cmd.append("--no-stats")
        if params.cmap != "inferno":
            cmd.extend(["--cmap", params.cmap])
        if params.figsize != [6, 8]:
            cmd.extend(["--figsize", str(params.figsize[0]), str(params.figsize[1])])
        if params.fontsize != 10:
            cmd.extend(["--fontsize", str(params.fontsize)])
        if not params.annot:
            cmd.append("--no-annot")
        if params.fmt != ".2f":
            cmd.extend(["--fmt", params.fmt])
        if not params.square:
            cmd.append("--no-square")
        skip_samples = config.get("skip_samples", [])
        if skip_samples:
            cmd.extend(["--skip-samples"] + skip_samples)
        
        cmd_str = " ".join(cmd)
        log_file = os.path.join("logs", f"heatmap_combined_{params.batch}.log")
        print(f"Generating combined heatmap: {output.combined}")
        shell(cmd_str, log=log_file)
# ============================================
# RULE: all_heatmaps (updated)
# ============================================
rule all_heatmaps:
    input:
        os.path.join(OUTPUT_DIR, "heatmap", "_all_heatmaps_generated.txt"),
        OUTPUT_DIR + "/heatmap/" + BATCH_NAME + "/" + config.get("combined_output", "combined_heatmap") + ".pdf"
# ============================================
# 9. RULE: clean
# ============================================
rule clean_logo:
    run:
        # Remove logs directory and its contents if present
        shell("rm -rf logs")
        # Remove the 'done' text file if it exists
        shell("rm -f done")
        # Remove the 'exit_with_error' file if it exists
        shell("rm -f exit_with_error")
        # Remove image output directory
        shell("rm -rf " + shlex.quote(OUTPUT_DIR))

rule clean_tex:
    run:
        # Remove LaTeX auxiliary files
        shell("rm -f " + "*.aux " + "*.log " + "*.out")
        # Remove generated LaTeX files
        shell("rm -f " + shlex.quote(LATEX_OUTPUT + "_" + str(MAX_K) + "k.tex"))
        shell("rm -f " + shlex.quote(LATEX_OUTPUT + "_" + str(MAX_K) + "k.pdf"))

rule clean_heatmap:
    run:
        # Remove heatmap output
        shell("rm -rf " + OUTPUT_DIR + "/heatmap")

# ============================================
# 10. RULE: clean_all (includes removing all output)
# ============================================
rule clean_all:
    run:
        # Remove logs directory and its contents if present
        shell("rm -rf logs")
        # Remove the 'done' text file if it exists
        shell("rm -f done")
        # Remove the 'exit_with_error' file if it exists
        shell("rm -f exit_with_error")
        # Remove image output directory
        shell("rm -rf " + shlex.quote(OUTPUT_DIR))
        # Remove temporary and LaTeX auxiliary files
        shell("rm -rf " + shlex.quote(OUTPUT_DIR + "/*/*.tmp"))
        shell("rm -f " + "*.aux " + "*.log " + "*.out")
        # Remove generated LaTeX files
        shell("rm -f *.tex *.pdf")