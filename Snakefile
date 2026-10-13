# Snakefile

import os
import glob
import shutil
import tempfile
import subprocess
import shlex
import time
# ============================================
# 0. USAGE INSTRUCTIONS
# ============================================

# # Run the pipeline
# (uv run) snakemake --core 4

# # Clean up logo images related outputs
# snakemake -j 4 clean_logo

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
    return os.path.abspath(bg_file) if bg_file else None

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
        expand(os.path.join(BATCH_LOGO_DIR, "{sample}/image/{sample}.info.eps"), sample=SAMPLES)
# ============================================
# 4. RULE: run_kpLogo #Add minimal mode where all output except .plogo is kept for ease of organisation
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
            print(f"Cleaned up: {job_root}")
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