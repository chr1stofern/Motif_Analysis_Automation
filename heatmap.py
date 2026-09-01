#!/usr/bin/env python3
"""
Generate heatmaps from kpLogo output files.

Usage:
    # Single sample from its pass.cutoff file
    python heatmap.py <pass.cutoff.txt> <output.pdf> [options]
    
    # Multiple samples from multiple pass.cutoff files
    python heatmap.py <pass1.cutoff.txt> <pass2.cutoff.txt> ... <output.pdf> [options]
    
    # Multiple samples from a samples file (one per line)
    python heatmap.py @samples.txt <output.pdf> [options]

Arguments:
    input_files   One or more .pass.p.cutoff.txt files, or a file with sample names (use @ prefix)
    output_pdf    Output PDF file path (e.g., kplogo_output/heatmap/Position-1.pdf)

Options:
    -p POSITION, --position POSITION   Position to extract (default: 1)
    -t TITLE, --title TITLE            Custom plot title (optional, auto-generated if not provided)
    --no-stats                         Use -log10(p-value) instead of test statistics
    --cmap CMAP                        Matplotlib colormap (default: inferno)
    --vmin VMIN                        Minimum color scale value (default: auto)
    --vmax VMAX                        Maximum color scale value (default: auto)
    --figsize W H                      Figure size in inches (default: 6 8)
    --fontsize FONTSIZE                Font size for annotations (default: 10)
    --no-annot                         Hide values in each cell
    --fmt FMT                          Number format for annotations (default: .2f)
    --no-square                        Don't force square cells
    -h, --help                         Show this help message
"""

import os
import sys
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import random
random.seed(42)
np.random.seed(42)

matplotlib.rcParams['font.family'] = 'Arial'

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate heatmaps from kpLogo output files.",
        usage="python heatmap.py <input_files...> <output_pdf> [options]"
    )

    # Positional: files and output PDF
    parser.add_argument("files", nargs="*",
                        help="One or more .pass.p.cutoff.txt files, or @samples.txt file (ignored if --pref-csv)")
    parser.add_argument("output_pdf", nargs="?", default=None,
                        help="Output PDF file path")

    # Plot options
    parser.add_argument("-p", "--position", type=int, default=1,
                        help="Position to extract (default: 1)")
    parser.add_argument("-t", "--title", default=None,
                        help="Custom plot title (auto-generated if not provided)")
    parser.add_argument("--no-stats", dest="use_stats", action="store_false", default=True,
                        help="Use -log10(p-value) instead of test statistics")
    parser.add_argument("--cmap", default="inferno",
                        help="Matplotlib colormap (default: inferno)")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--figsize", nargs=2, type=float, default=[8, 8],
                        help="Figure size in inches (default: 8 8)")
    parser.add_argument("--fontsize", type=int, default=10)
    parser.add_argument("--no-annot", dest="annot", action="store_false", default=False)
    parser.add_argument("--fmt", default=".2f")
    parser.add_argument("--no-square", dest="square", action="store_false", default=True)

    # Efficiency
    parser.add_argument("-e", "--eff-dir", default="eff",
                        help="Directory containing efficiency tables (default: eff)")
    parser.add_argument("--eff-file", default=None,
                        help="Explicit efficiency file (overrides auto-detection)")

    # Cell size
    parser.add_argument("--cell-size", type=float, default=0.4,
                        help="Cell size in inches (default: 0.4)")

    # Skipping samples
    parser.add_argument("--skip-samples", nargs="+", default=[],
                        help="List of sample names to skip")
    parser.add_argument("--skip-file", default=None,
                        help="File containing sample names to skip (one per line)")

    # Batch and global scales
    parser.add_argument("--batch", type=str, default=None,
                        help="Batch name (e.g., d129) for computing global scales")
    parser.add_argument("--global-vmin", type=float, default=None,
                        help="Global vmin for preference heatmap (overrides auto)")
    parser.add_argument("--global-vmax", type=float, default=None,
                        help="Global vmax for preference heatmap (overrides auto)")
    parser.add_argument("--eff-global-vmin", type=float, default=None,
                        help="Global vmin for efficiency strip (overrides auto)")
    parser.add_argument("--eff-global-vmax", type=float, default=None,
                        help="Global vmax for efficiency strip (overrides auto)")

    # Combine positions
    parser.add_argument("--combine", nargs="+", type=int, default=None,
                        help="Combine multiple positions (e.g., --combine 3 4 5)")

    # NEW: Pre‑extracted preference CSV
    parser.add_argument("--pref-csv", type=str, default=None,
                        help="Use pre‑extracted preference CSV instead of pass.cutoff files")

    args = parser.parse_args()

    # Validate arguments
    if args.pref_csv is None:
        # Need at least one input file and an output PDF
        if not args.files:
            parser.error("At least one input file is required when --pref-csv is not used")
        if args.output_pdf is None:
            parser.error("output_pdf is required when --pref-csv is not used")
    else:
        # With --pref-csv, the output PDF must be provided either as the last positional
        if args.files and args.output_pdf is None:
            args.output_pdf = args.files[-1]
            args.files = args.files[:-1] if len(args.files) > 1 else []
        elif args.output_pdf is None:
            parser.error("output_pdf is required")
    return args

def ensure_output_dir(output_path):
    """Create the output directory if it doesn't exist."""
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

def load_preference_csv(csv_path, skip_samples=None, samples=None):
    df = pd.read_csv(csv_path)
    if skip_samples:
        df = df[~df['sample'].isin(skip_samples)]
    if samples is not None:
        df = df[df['sample'].isin(samples)]
    positions = df['pos'].unique()
    pref_dfs = {}
    # Define desired order and naming (e.g., RNA)
    base_order = ['A', 'U', 'C', 'G']   # change as needed
    for pos in positions:
        pos_df = df[df['pos'] == pos].copy()
        pivoted = pos_df.pivot_table(index='sample', columns='base', values='value', aggfunc='first')
        # Rename T to U if desired
        if 'T' in pivoted.columns:
            pivoted.rename(columns={'T': 'U'}, inplace=True)
        # Ensure all bases are present
        for b in base_order:
            if b not in pivoted.columns:
                pivoted[b] = 0
        pivoted = pivoted[base_order]   # reorder columns
        # ===== CONTROL THE X-AXIS LABEL =====
        # pivoted.columns.name = ""   # or "" to remove, or "Base"
        pivoted.columns.name = None
        pref_dfs[pos] = pivoted
    return pref_dfs

def compute_global_scales_from_dfs(pref_dfs):
    """Compute global vmin/vmax from a dict of DataFrames."""
    all_values = []
    for df in pref_dfs.values():
        all_values.extend(df.values.flatten())
    if all_values:
        max_abs = max(abs(np.nanmin(all_values)), abs(np.nanmax(all_values))) or 1.0
        return -max_abs, max_abs
    return None, None

def is_pass_file(filename):
    """Check if the file is a .pass.p.cutoff.txt file."""
    return filename.endswith(".pass.p.cutoff.txt")

def expand_input_files(files):
    """
    Expand input files list.
    
    Handles:
    - .pass.p.cutoff.txt files directly
    - @samples.txt files (one sample per line)
    - Glob patterns
    """
    expanded = []
    
    for f in files:
        # If file starts with '@', it's a samples file
        if f.startswith('@'):
            samples_file = f[1:]  # Remove '@'
            if not os.path.exists(samples_file):
                raise FileNotFoundError(f"Samples file not found: {samples_file}")
            
            with open(samples_file, "r") as sf:
                for line in sf:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        expanded.append(line)
            continue
        
        # If it's a glob pattern (contains wildcards)
        if '*' in f or '?' in f or '[' in f or ']' in f:
            matched = glob.glob(f)
            if not matched:
                print(f"Warning: No files found matching pattern: {f}", file=sys.stderr)
            expanded.extend(matched)
            continue
        
        # Regular file
        expanded.append(f)
    
    return expanded

def get_sample_from_pass_file(pass_file):
    """Extract sample name from a .pass.p.cutoff.txt file."""
    basename = os.path.basename(pass_file)
    return basename.replace(".pass.p.cutoff.txt", "")

def infer_data_dir(pass_file):
    """Infer data directory from a pass.cutoff file path."""
    parts = Path(pass_file).parts
    # Find the 'data' directory and go up one level (to the sample directory)
    # Then go up one more level (to the root directory)
    for i, part in enumerate(parts):
        if part == "data" and i > 1:
            # Go up two levels: from 'data' to the root directory
            return os.path.join(*parts[:i-1])
    # Fallback: go up three levels from the file
    return os.path.dirname(os.path.dirname(os.path.dirname(pass_file)))

def load_samples(files):
    """
    Load samples from input files.
    
    Returns:
        tuple: (samples, data_dir)
            samples: list of sample names
            data_dir: inferred from the first pass.cutoff file
    """
    # Expand input files
    expanded = expand_input_files(files)
    
    if not expanded:
        raise RuntimeError("No input files found")
    
    samples = []
    data_dir = None
    
    for f in expanded:
        if is_pass_file(f):
            sample = get_sample_from_pass_file(f)
            samples.append(sample)
            if data_dir is None:
                data_dir = infer_data_dir(f)
        else:
            print(f"Warning: {f} is not a .pass.p.cutoff.txt file. Skipping.", file=sys.stderr)
    
    if not samples:
        raise RuntimeError("No valid .pass.p.cutoff.txt files found")
    
    return samples, data_dir

def extract_position_scores(sample, data_dir, position, use_stats=True):
    """
    Extract scores for a specific position from pass.p.cutoff.txt.
    
    Returns None if the position is not found (no fallback to other positions).
    """
    pass_file = os.path.join(data_dir, sample, "data", f"{sample}.pass.p.cutoff.txt")
    
    if not os.path.exists(pass_file):
        print(f"Warning: {pass_file} not found", file=sys.stderr)
        return None
    
    try:
        df = pd.read_csv(pass_file, sep='\t', comment='#', header=None)
    except pd.errors.EmptyDataError:
        print(f"Warning: {pass_file} is empty", file=sys.stderr)
        return None
    
    # Try different column formats
    try:
        df.columns = ['kmer', 'pos', 'shift', 'statistics', 'pvalue', 
                      'corrected_p', 'frac_obs', 'frac_exp', 'obs_exp', 'local_r']
    except ValueError:
        try:
            df.columns = ['kmer', 'pos', 'shift', 'statistics', 'pvalue', 
                          'corrected_p', 'n1', 'mean1', 'std1', 'n2', 'mean2', 'std2']
        except ValueError:
            print(f"Warning: Unexpected column count in {pass_file}", file=sys.stderr)
            return None
    
    # Filter to k=1 (monomers only)
    kmer_df = df[df['kmer'].isin(['A', 'T', 'C', 'G'])].copy()
    
    if kmer_df.empty:
        print(f"Warning: No k=1 entries in {pass_file}", file=sys.stderr)
        return None
    
    # Map position - EXACT match only, NO FALLBACK
    pos_df = kmer_df[kmer_df['pos'] == position]
    
    if pos_df.empty:
        # Position not found - return None (will be filled with 0 later)
        print(f"Warning: Position {position} not found for {sample}. No data available.", file=sys.stderr)
        return None
    
    # Extract scores
    scores = {}
    for _, row in pos_df.iterrows():
        base = row['kmer']
        if use_stats:
            scores[base] = row['statistics']
        else:
            scores[base] = row['pvalue']
    
    return scores

def infer_batch(samples):
    """
    Infer batch token from sample names.
    Examples:
        3-Psp-D129A  -> D129
        10-Psp-N1060K -> N1060
    """
    import re
    tokens = []
    for s in samples:
        # look for patterns like D129, N1060, T136, etc.
        m = re.search(r'([A-Za-z]\d{2,})', s)
        if m:
            tokens.append(m.group(1).upper())
    if not tokens:
        return None
    # most common token
    from collections import Counter
    return Counter(tokens).most_common(1)[0][0]

def find_eff_file(eff_dir, batch, explicit=None):
    """Return path to the efficiency file (now CSV format)."""
    if explicit:
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"Efficiency file not found: {explicit}")
        return explicit

    if batch is None:
        print("Warning: could not infer batch; skipping efficiency panel.", file=sys.stderr)
        return None

    # Try new naming convention first: eff_{batch}.csv
    candidates = [
        os.path.join(eff_dir, f"eff_{batch.lower()}.csv"),  # eff_d129.csv
        os.path.join(eff_dir, f"eff_{batch.upper()}.csv"),  # eff_D129.csv
        os.path.join(eff_dir, f"psp_{batch.lower()}_eff.txt"),  # psp_d129_eff.txt (old)
        os.path.join(eff_dir, f"psp_{batch.upper()}_eff.txt"),  # psp_D129_eff.txt (old)
        os.path.join(eff_dir, f"{batch.lower()}_eff.txt"),
        os.path.join(eff_dir, f"{batch.upper()}_eff.txt"),
        os.path.join(eff_dir, f"*{batch.lower()}*eff*.csv"),
        os.path.join(eff_dir, f"*{batch.upper()}*eff*.csv"),
        os.path.join(eff_dir, "*.csv"),  # Fallback: any CSV file
    ]
    
    import glob
    for pat in candidates:
        matches = glob.glob(pat)
        if matches:
            return matches[0]
    
    print(f"Warning: no efficiency file found for batch '{batch}' in '{eff_dir}'.",
          file=sys.stderr)
    return None

def pipeline_to_eff_key(sample):
    """
    Map pipeline sample name -> key used in efficiency table (CSV format).
    
    Examples:
        1-Psp-WT-D129NT   -> WT+crNT
        2-Psp-WT-D129EGFP7 -> WT+crEGFP7
        3-Psp-D129A        -> D129A
        10-Psp-D129K       -> D129K
    """
    # strip leading number and optional "Psp-"
    parts = sample.split("-", 1)
    rest = parts[1] if len(parts) > 1 else sample
    if rest.startswith("Psp-"):
        rest = rest[4:]

    # special WT cases - match CSV format EXACTLY
    if "WT" in rest:
        if "NT" in rest:
            return "WT+crNT"
        elif "EGFP7" in rest:
            return "WT+crEGFP7"

    # For mutants, return the last part
    # 3-Psp-D129A -> D129A
    # 14-Psp-N1060P -> N1060P
    # 15-Psp-T136P -> T136P
    # Also handle cases where the batch name is in the middle (e.g., WT-D129NT)
    # If it's something like "D129A", return it directly
    if rest in ["D129A", "D129C", "D129E", "D129F", "D129G", "D129H", "D129I", 
                "D129K", "D129L", "D129M", "D129N", "D129P", "D129Q", "D129R", 
                "D129S", "D129T", "D129V", "D129W", "D129Y"]:
        return rest
    
    # Fallback: return the last part after splitting by "-"
    return rest.split("-")[-1]

def load_efficiency(eff_path, samples):
    """
    Load efficiency table from CSV format (eff_d129.csv, etc.).
    Returns a tuple of (efficiency_series, pvalue_series) aligned to `samples`.
    """
    if eff_path is None or not os.path.exists(eff_path):
        print(f"Warning: Efficiency file not found: {eff_path}", file=sys.stderr)
        return None, None

    try:
        df = pd.read_csv(eff_path)
        
        print(f"Loaded efficiency file: {eff_path}")
        print(f"Columns: {df.columns.tolist()}")
        print(f"Rows: {len(df)}")
        print(f"Mutants in file: {df['Mutants'].tolist()}")
        
        # Create mapping from mutant name to efficiency and p-value
        eff_map = {}
        pval_map = {}
        
        for _, row in df.iterrows():
            mutant = str(row['Mutants']).strip()
            # eff_map[mutant] = row['Inverted_Log2_FC']
            eff_map[mutant] = row['Inverted_Log2_FC_Norm']
            if 'Dunnett_P_Value' in df.columns:
                pval_map[mutant] = row['Dunnett_P_Value']
            else:
                pval_map[mutant] = np.nan
        
    except Exception as e:
        print(f"Warning: Could not read {eff_path}: {e}", file=sys.stderr)
        return None, None

    if not eff_map:
        print(f"Warning: No efficiency data found in {eff_path}", file=sys.stderr)
        return None, None

    values = []
    pvalues = []
    missing_samples = []
    
    for s in samples:
        key = pipeline_to_eff_key(s)
        
        # Try different key formats
        eff_val = eff_map.get(key, np.nan)
        p_val = pval_map.get(key, np.nan)
        
        if np.isnan(eff_val):
            # Try using the last part of the sample name
            alt_key = s.split("-")[-1] if "-" in s else s
            eff_val = eff_map.get(alt_key, np.nan)
            p_val = pval_map.get(alt_key, np.nan)
            
            if np.isnan(eff_val):
                missing_samples.append(f"{s} (looked up as '{key}' and '{alt_key}')")
        
        values.append(eff_val)
        pvalues.append(p_val)
    
    if missing_samples:
        print(f"Warning: {len(missing_samples)} samples not found in efficiency file:")
        for m in missing_samples[:5]:  # Show first 5
            print(f"  - {m}")
        if len(missing_samples) > 5:
            print(f"  ... and {len(missing_samples) - 5} more")
    
    eff_series = pd.Series(values, index=samples, name="Efficiency")
    pval_series = pd.Series(pvalues, index=samples, name="PValue")
    
    print(f"Loaded efficiency for {eff_series.count()} out of {len(samples)} samples")
    
    return eff_series, pval_series

def extract_mutant_name(sample: str) -> str:
    parts = sample.split('-')
    if len(parts) >= 3 and parts[2] == "WT":
        return "WT"
    return parts[2] if len(parts) >= 3 else sample

def shorten_sample_name(sample: str) -> str:
    # Example Names
    # "1-Psp-WT-NT"
    # "2-Psp-WT-EGFP7"
    # "3-Psp-D129A"
    # "10-Psp-D129K"
    # "21-Psp-D129Y"

    # split on "-"
    parts = sample.split("-")          # ["3", "Psp", "D129A"]

    if parts[2] == "WT":
        # shortened = parts[2] + parts[3]   # e.g. "WT-NT"
        shortened = parts[2]  # e.g. "WT-NT"
    else:
        shortened = parts[2][-1] # "D129A"

    # everything after the first "-"
    # sample.split("-", 1)[-1] # "Psp-D129A"
    # short = rest.replace("Psp-", "")   # "D129A"

    return shortened

def create_heatmap_data(samples, data_dir, position, use_stats=True):
    """
    Create a DataFrame for the heatmap.
    
    Missing positions are filled with 0.
    """
    
    #change label here for DNA or RNA motifs
    all_bases = ['A', 'U', 'C', 'G']
    
    heatmap_data = {}
    missing_samples = []
    
    for sample in samples:
        scores = extract_position_scores(sample, data_dir, position, use_stats)
        if scores is not None:
            # Ensure all bases are present
            for base in all_bases:
                if base not in scores:
                    scores[base] = 0
            heatmap_data[sample] = scores
        else:
            # Sample has no data for this position
            # Add all bases with value 0
            heatmap_data[sample] = {base: 0 for base in all_bases}
            missing_samples.append(sample)
    
    if missing_samples:
        print(f"Warning: Missing position data for samples: {missing_samples}. Filled with 0.", file=sys.stderr)
    
    if not heatmap_data:
        raise RuntimeError("No data found for any samples")
    
    # Convert to DataFrame and reorder columns
    df = pd.DataFrame(heatmap_data).T
    df = df[all_bases]
    
    return df

def compute_global_scales(samples, data_dir, position=1, use_stats=True):
    """
    Compute global vmin/vmax for all positions.
    Returns (pref_vmin, pref_vmax)
    """
    import numpy as np
    
    print("Computing global color scales from all positions...")
    
    # Determine positions to scan (1-8)
    positions = list(range(1, 9))
    all_pref_values = []
    
    for pos in positions:
        try:
            df = create_heatmap_data(samples, data_dir, pos, use_stats)
            all_pref_values.extend(df.values.flatten())
        except Exception as e:
            print(f"Warning: Could not compute for position {pos}: {e}", file=sys.stderr)
    
    if all_pref_values:
        pref_vmin = np.nanmin(all_pref_values)
        pref_vmax = np.nanmax(all_pref_values)
        # Make symmetric around 0
        max_abs = max(abs(pref_vmin), abs(pref_vmax)) or 1.0
        pref_vmin, pref_vmax = -max_abs, max_abs
        print(f"Global preference scale: vmin={pref_vmin:.3f}, vmax={pref_vmax:.3f}")
    else:
        pref_vmin, pref_vmax = None, None
        print("Warning: No preference values found for global scale.", file=sys.stderr)
    
    return pref_vmin, pref_vmax

def plot_combined(pref_df, eff_series, pval_series, output_file, title=None, position=None, **kwargs):
    """
    Layout:
        Title (centered, spanning both columns)
        Preference Heatmap (left) | Efficiency Strip (right)
        Preference Colorbar (full width, below both)
        Efficiency Colorbar (full width, below preference colorbar)
    """
    ensure_output_dir(output_file)

    if title is None:
        title = f"Base Preference at Position {position}" if position is not None else "Base Preference Heatmap"

    # Sort rows
    def extract_number(s):
        try:
            return int(s.split("-")[0])
        except Exception:
            return 999

    # order = sorted(pref_df.index, key=extract_number)
    order = sorted(pref_df.index, key=extract_mutant_name)
    pref_df = pref_df.loc[order]
    
    # Handle None case for efficiency
    has_efficiency = eff_series is not None and not eff_series.isna().all()
    if has_efficiency:
        eff_series = eff_series.reindex(order)
        if pval_series is not None:
            pval_series = pval_series.reindex(order)
    else:
        eff_series = None
        pval_series = None

    ylabels = [shorten_sample_name(s) for s in pref_df.index]

    # ===== COLOR LIMITS (USE GLOBAL VALUES IF PROVIDED) =====
    vmin = kwargs.get("vmin", None)
    vmax = kwargs.get("vmax", None)
    if vmin is None or vmax is None:
        # Fallback to per-position computation
        max_abs = max(abs(pref_df.min().min()), abs(pref_df.max().max())) or 1.0
        vmin, vmax = -max_abs, max_abs
        print(f"Using per-position scale: vmin={vmin:.3f}, vmax={vmax:.3f}")

    # ===== EFFICIENCY COLOR LIMITS =====
    if has_efficiency:
        eff_vmin = kwargs.get("eff_vmin", None)
        eff_vmax = kwargs.get("eff_vmax", None)
        if eff_vmin is None or eff_vmax is None:
            eff_max = eff_series.max()
            eff_min = eff_series.min()
            eff_abs_max = max(abs(eff_min), abs(eff_max)) or 0.1
            eff_vmin, eff_vmax = -eff_abs_max, eff_abs_max
            print(f"Using per-position efficiency scale: vmin={eff_vmin:.3f}, vmax={eff_vmax:.3f}")

    # ===== CALCULATE FIGURE SIZE =====
    n_rows = len(pref_df)
    n_cols = len(pref_df.columns)
    
    cell_size = kwargs.get("cell_size", 0.4)
    
    fig_width = (n_cols + 1) * cell_size + 0.3
    fig_height = n_rows * cell_size + 0.9
    
    print(f"Figure size: {fig_width:.2f} x {fig_height:.2f} inches (cell_size={cell_size})")
    
    fig = plt.figure(figsize=(fig_width, fig_height))

    # ===== LAYOUT WITH TOP SPACER =====
    if has_efficiency:
        # Full layout with efficiency strip
        axd = fig.subplot_mosaic(
            [
                ["top_spacer", "top_spacer"],
                ["heat", "eff"],
                ["empty", "empty"],
                ["cbar", "cbar"],
                ["cbar_spacer", "cbar_spacer"],
                ["eff_cbar", "eff_cbar"]
            ],
            width_ratios=[n_cols, 1],
            height_ratios=[0.01, 1.0, 0.04, 0.04, 0.06, 0.04],
            gridspec_kw={"wspace": 0.02, "hspace": 0.04}
        )
        # Hide spacer axes
        axd["top_spacer"].axis("off")
        axd["empty"].axis("off")
        axd["cbar_spacer"].axis("off")
    else:
        # Layout without efficiency strip (only heatmap)
        axd = fig.subplot_mosaic(
            [
                ["top_spacer"],
                ["heat"],
                ["empty"],
                ["cbar"]
            ],
            width_ratios=[1.0],
            height_ratios=[0.01, 1.0, 0.04, 0.06],
            gridspec_kw={"hspace": 0.04}
        )
        axd["top_spacer"].axis("off")
        axd["empty"].axis("off")

    # ===== TITLE =====
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.95)

    # ===== PREFERENCE HEATMAP =====
    sns.heatmap(
        pref_df,
        ax=axd["heat"],
        cmap=kwargs.get("cmap", "inferno"),
        center=0,
        vmin=vmin,
        vmax=vmax,
        annot=kwargs.get("annot", False),
        fmt=kwargs.get("fmt", ".2f"),
        square=True,
        yticklabels=ylabels,
        cbar=True,
        cbar_ax=axd["cbar"] if has_efficiency else axd["cbar"],
        cbar_kws={
            "label": "Residue Preference (z-score)",
            "orientation": "horizontal"
        },
        annot_kws={"fontsize": kwargs.get("fontsize", 10)}
    )
    axd["heat"].set_ylabel("")

    # ===== EFFICIENCY STRIP WITH ASTERISKS =====
    if has_efficiency and "eff" in axd:
        eff_df = eff_series.to_frame(name="")

        sns.heatmap(
            eff_df,
            ax=axd["eff"],
            cmap="inferno",
            center=0,
            vmin=eff_vmin,
            vmax=eff_vmax,
            annot=False,
            fmt=".2f",
            square=True,
            xticklabels=False,
            yticklabels=False,
            cbar=True,
            cbar_ax=axd["eff_cbar"],
            cbar_kws={
                "label": "Knockdown Efficiency (-log Fold Change in Fluorescence)",
                "orientation": "horizontal"
            }
        )
        axd["eff"].set_xlabel("Eff.", fontsize=11)
        axd["eff"].set_ylabel("")

        # ===== ADD ASTERISKS FOR SIGNIFICANT MUTANTS =====
        if pval_series is not None:
            for i, (idx, val) in enumerate(eff_series.items()):
                p_val = pval_series[idx] if idx in pval_series.index else np.nan
                
                # Check if significant (p < 0.05)
                is_significant = not np.isnan(p_val) and p_val < 0.05
                
                # Get the color of the cell to determine text color
                eff_val = eff_series[idx]
                if not np.isnan(eff_val):
                    # norm_val = (eff_val - eff_vmin) / (eff_vmax - eff_vmin) if eff_vmax != eff_vmin else 0.5
                    # text_color = 'white' if norm_val > 0.6 or norm_val < 0.4 else 'black'
                    text_color = 'white' if eff_val < 40 else 'black'
                    
                    # Add asterisk
                    if is_significant:
                        axd["eff"].text(
                            0.5, i + 0.65,
                            "*",
                            ha='center', va='center',
                            fontsize=12,
                            fontweight='bold',
                            color=text_color
                            # color='white'
                        )

    # ===== PADDING =====
    if has_efficiency:
        plt.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.06)
    else:
        plt.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.06)

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Combined heatmap saved to: {output_file}")
    plt.close(fig)

def plot_combined_positions(pref_dfs, eff_series, pval_series, output_file, 
                            positions, title=None, **kwargs):
    """
    Combine multiple positions into a single figure.
    Efficiency strip appears only once on the right.
    Single combined colorbar with two labels (top: z-score, bottom: efficiency).
    """
    ensure_output_dir(output_file)
    
    if title is None:
        title = f"Positions {', '.join(map(str, positions))}"

    # Sort rows
    def extract_number(s):
        try:
            return int(s.split("-")[0])
        except Exception:
            return 999

    first_pos = positions[0]
    pref_df_first = pref_dfs[first_pos]
    # order = sorted(pref_df_first.index, key=extract_number)
    order = sorted(pref_df_first.index, key=extract_mutant_name)
    
    pref_dfs_ordered = {}
    for pos, df in pref_dfs.items():
        pref_dfs_ordered[pos] = df.loc[order]
    
    has_efficiency = eff_series is not None and not eff_series.isna().all()
    if has_efficiency:
        eff_series = eff_series.reindex(order)
        if pval_series is not None:
            pval_series = pval_series.reindex(order)
    
    ylabels = [shorten_sample_name(s) for s in order]

    # ===== COLOR LIMITS =====
    vmin = kwargs.get("vmin", None)
    vmax = kwargs.get("vmax", None)
    if vmin is None or vmax is None:
        all_values = []
        for df in pref_dfs_ordered.values():
            all_values.extend(df.values.flatten())
        max_abs = max(abs(np.nanmin(all_values)), abs(np.nanmax(all_values))) or 1.0
        vmin, vmax = -max_abs, max_abs
        print(f"Combined preference scale: vmin={vmin:.3f}, vmax={vmax:.3f}")

    if has_efficiency:
        eff_vmin = kwargs.get("eff_vmin", None)
        eff_vmax = kwargs.get("eff_vmax", None)
        if eff_vmin is None or eff_vmax is None:
            eff_max = eff_series.max()
            eff_min = eff_series.min()
            eff_abs_max = max(abs(eff_min), abs(eff_max)) or 0.1
            eff_vmin, eff_vmax = -eff_abs_max, eff_abs_max
            print(f"Combined efficiency scale: vmin={eff_vmin:.3f}, vmax={eff_vmax:.3f}")

    # ===== CALCULATE FIGURE SIZE =====
    n_rows = len(pref_df_first)
    n_cols = len(pref_df_first.columns)
    n_positions = len(positions)
    
    cell_size = kwargs.get("cell_size", 0.4)
    spacer_width = kwargs.get("spacer_width", 0.15)
    n_cols_total = n_positions * n_cols + 2
    
    fig_width = (n_positions * n_cols + 1) * cell_size + spacer_width + 0.5
    fig_height = n_rows * cell_size + 1.2
    
    print(f"Figure size: {fig_width:.2f} x {fig_height:.2f} inches (cell_size={cell_size})")
    print(f"Total columns: {n_cols_total}")
    
    fig = plt.figure(figsize=(fig_width, fig_height))

    # ===== SPACING PARAMETERS =====
    wspace = kwargs.get("wspace", 0.03)
    hspace = kwargs.get("hspace", 0.04)
    pad_left = kwargs.get("pad_left", 0.05)
    pad_right = kwargs.get("pad_right", 0.95)
    pad_top = kwargs.get("pad_top", 0.89)
    
    # ===== HEIGHT RATIOS =====
    spacer_height = kwargs.get("spacer_height", 0.01 )
    cbar_spacer_height = kwargs.get("cbar_spacer_height", 0.1)
    cbar_height = kwargs.get("cbar_height", 0.05)  # Moderate height
    
    height_ratios = [
        spacer_height,              # Top spacer
        1.0,                        # Heatmap + efficiency row
        cbar_spacer_height,         # Spacer between heatmaps and colorbars
        cbar_height                 # Single combined colorbar
    ]
    
    # ===== LAYOUT WITH SPACER =====
    row_names = []
    for pos in positions:
        for _ in range(n_cols):
            row_names.append("heat_" + str(pos))
    row_names.append("spacer_between")
    row_names.append("eff")
    
    width_ratios = []
    for _ in positions:
        for _ in range(n_cols):
            width_ratios.append(1)
    width_ratios.append(spacer_width)
    width_ratios.append(1)
    
    mosaic = [
        ["top_spacer"] * n_cols_total,
        row_names,
        ["cbar_spacer"] * n_cols_total,
        ["cbar"] * n_cols_total
    ]
    
    axd = fig.subplot_mosaic(
        mosaic,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
        gridspec_kw={"wspace": wspace, "hspace": hspace}
    )
    
    # Hide spacer axes
    axd["top_spacer"].axis("off")
    axd["cbar_spacer"].axis("off")
    axd["spacer_between"].axis("off")

    # ===== TITLE =====
    fig.suptitle(title, fontsize=18, fontweight='bold', y=0.97)

    # ===== PREFERENCE HEATMAPS =====
    show_spines = kwargs.get("show_spines", True)
    spine_color = kwargs.get("spine_color", "black")
    spine_width = kwargs.get("spine_width", 0.5)
    
    for i, pos in enumerate(positions):
        heat_ax = axd["heat_" + str(pos)]
        heat_data = pref_dfs_ordered[pos]
        
        sns.heatmap(
            heat_data,
            ax=heat_ax,
            cmap=kwargs.get("cmap", "inferno"),
            center=0,
            vmin=vmin,
            vmax=vmax,
            annot=kwargs.get("annot", False),
            fmt=kwargs.get("fmt", ".2f"),
            square=True,
            yticklabels=ylabels if i == 0 else False,
            xticklabels=True,
            cbar=False,
            annot_kws={"fontsize": kwargs.get("fontsize", 10)}
        )
        # ===== TITLE WITH RELATIVE POSITION (Dictionary Mapping) =====
        position_map = {
            3: "-2 Position",
            4: "-1 Position",
            5: "+1 Position"
        }
        heat_ax.set_title(position_map.get(pos, f"Position {pos}"), fontsize=12)
        heat_ax.set_ylabel("")
        
        # ===== HEATMAP SPINES =====
        if show_spines:
            heat_ax.spines["left"].set_visible(True)
            heat_ax.spines["bottom"].set_visible(True)
            heat_ax.spines["top"].set_visible(True)
            heat_ax.spines["right"].set_visible(True)
            

    # ===== SINGLE EFFICIENCY STRIP =====
    if has_efficiency:
        eff_df = eff_series.to_frame(name="")
        eff_ax = axd["eff"]
        
        sns.heatmap(
            eff_df,
            ax=eff_ax,
            cmap=kwargs.get("cmap", "inferno"),
            center=0,
            vmin=eff_vmin,
            vmax=eff_vmax,
            annot=False,
            fmt=".2f",
            square=True,
            xticklabels=False,
            yticklabels=False,
            cbar=False
        )
        eff_ax.set_xlabel("Eff.", fontsize=11)
        eff_ax.set_ylabel("")
        
        # ===== EFFICIENCY STRIP SPINES =====
        if show_spines:
            eff_ax.spines["left"].set_visible(True)
            eff_ax.spines["right"].set_visible(True)
            eff_ax.spines["bottom"].set_visible(True)
            eff_ax.spines["top"].set_visible(True)
            
            for spine in eff_ax.spines.values():
                if spine.get_visible():
                    spine.set_color(spine_color)
                    spine.set_linewidth(spine_width)
        
        # ===== ASTERISKS =====
        if pval_series is not None:
            for j, (idx, val) in enumerate(eff_series.items()):
                p_val = pval_series[idx] if idx in pval_series.index else np.nan
                is_significant = not np.isnan(p_val) and p_val < 0.05
                
                if is_significant:
                    eff_val = eff_series[idx]
                    if not np.isnan(eff_val):
                        # norm_val = (eff_val - eff_vmin) / (eff_vmax - eff_vmin) if eff_vmax != eff_vmin else 0.5
                        # text_color = 'white' if abs(norm_val - 0.5) > 0.2 else 'black'
                        text_color = 'white' if eff_val < 40 else 'black'
                        
                        eff_ax.text(
                            0.5, j + 0.65,
                            "*",
                            ha='center', va='center',
                            fontsize=12,
                            fontweight='bold',
                            color=text_color
                        )

    # ===== SINGLE SHARED COLORBAR WITH TWO LABELS =====
    # 1. Place cbar for eff map at the bottom
    sm = plt.cm.ScalarMappable(cmap=kwargs.get("cmap", "inferno"), norm=plt.Normalize(eff_vmin, eff_vmax))
    sm.set_array([])
    
    cbar = fig.colorbar(sm, cax=axd["cbar"], orientation="horizontal")
    
    # Configure Bottom Scale
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.set_label("Knockdown Efficiency (-log Fold Change in Fluorescence)", fontsize=10, labelpad=6)

    # 2. Place cbar for pref map at the top using twiny() (twinx doesnt work for horizontal bars)
    ax_twin = cbar.ax.twiny()
    
    # Match the geometric bounds of the colorbar, but map to efficiency data ranges
    ax_twin.set_xlim(vmin, vmax)
    
    # Configure Top Scale (Z-score or Preference scale)
    ax_twin.xaxis.set_label_position('top')
    ax_twin.xaxis.set_ticks_position('top')
    ax_twin.set_xlabel("Residue Preference (z-score)", fontsize=10, labelpad=6)

    # CRITICAL FIX: Explicitly ensure upper axis marks are visible 
    # (Seaborn often strips these settings globally on subplots)
    ax_twin.tick_params(axis='x', which='both', top=True, labeltop=True, direction='out')

    # ===== COLORBAR SPINES =====
    if show_spines:
        for spine in axd["cbar"].spines.values():
            spine.set_visible(True)
            spine.set_color(spine_color)
            spine.set_linewidth(spine_width)

    # ===== PADDING =====
    plt.subplots_adjust(
        left=pad_left,
        right=pad_right,
        top=pad_top,
        bottom=0.06
    )

    # ===== SAVE =====
    plt.savefig(output_file, dpi=300, bbox_inches="tight", pad_inches=0.1)
    print(f"Combined heatmap saved to: {output_file}")
    plt.close(fig)
          
def main():
    args = parse_args()

    # Collect skip samples from file if provided
    skip_set = set(args.skip_samples)
    if args.skip_file and os.path.exists(args.skip_file):
        with open(args.skip_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    skip_set.add(line)

    # ======================================================================
    # We require --pref-csv; no fallback to scanning pass files.
    # ======================================================================
    if args.pref_csv is None:
        raise RuntimeError("The --pref-csv argument is required. "
                           "Please generate preference_data.csv via the Snakefile rule.")
    if not os.path.exists(args.pref_csv):
        raise FileNotFoundError(f"Preference CSV not found: {args.pref_csv}")

    # Load all positions from CSV
    pref_dfs = load_preference_csv(args.pref_csv, skip_samples=skip_set)
    if not pref_dfs:
        raise RuntimeError("No data loaded from preference CSV")
    # Get samples from the first position (all positions share the same sample list)
    first_pos = sorted(pref_dfs.keys())[0]
    samples = list(pref_dfs[first_pos].index)
    print(f"Loaded preference data for {len(samples)} samples, positions {sorted(pref_dfs.keys())}")

    # Determine batch
    batch = args.batch
    if batch is None:
        batch = infer_batch(samples)
    print(f"Inferred batch: {batch}")

    # Load efficiency
    eff_path = args.eff_file
    if eff_path is None and batch is not None:
        import glob
        candidates = [
            os.path.join(args.eff_dir, f"eff_{batch.lower()}.csv"),
            os.path.join(args.eff_dir, f"eff_{batch.upper()}.csv"),
            os.path.join(args.eff_dir, f"psp_{batch.lower()}_eff.txt"),
            os.path.join(args.eff_dir, f"psp_{batch.upper()}_eff.txt"),
        ]
        for c in candidates:
            if os.path.exists(c):
                eff_path = c
                break
    if eff_path is None:
        print("Warning: No efficiency file found. Continuing without efficiency strip.", file=sys.stderr)
    eff_series, pval_series = load_efficiency(eff_path, samples) if eff_path else (None, None)

    # Compute global preference scale from all positions
    if args.global_vmin is None or args.global_vmax is None:
        pref_vmin, pref_vmax = compute_global_scales_from_dfs(pref_dfs)
        if pref_vmin is None:
            pref_vmin, pref_vmax = -1, 1  # fallback
        print(f"Global preference scale: vmin={pref_vmin:.3f}, vmax={pref_vmax:.3f}")
    else:
        pref_vmin, pref_vmax = args.global_vmin, args.global_vmax
        print(f"User specified preference scale: vmin={pref_vmin:.3f}, vmax={pref_vmax:.3f}")

    # Compute global efficiency scale
    if args.eff_global_vmin is None or args.eff_global_vmax is None:
        if eff_series is not None and not eff_series.isna().all():
            # eff_max = max(abs(eff_series.min()), abs(eff_series.max())) or 0.1
            eff_vmin = eff_series.min()
            eff_vmax = eff_series.max()
            # eff_vmin, eff_vmax = -eff_max, eff_max
            print(f"Global efficiency scale: vmin={eff_vmin:.3f}, vmax={eff_vmax:.3f}")
        else:
            eff_vmin, eff_vmax = None, None
    else:
        eff_vmin, eff_vmax = args.eff_global_vmin, args.eff_global_vmax
        print(f"User specified efficiency scale: vmin={eff_vmin:.3f}, vmax={eff_vmax:.3f}")

    # Decide: combine positions or single position?
    if args.combine and len(args.combine) > 0:
        positions = sorted(args.combine)
        # Ensure we have data for all requested positions
        for pos in positions:
            if pos not in pref_dfs:
                raise RuntimeError(f"Position {pos} not found in preference CSV")
        title = args.title or f"Cleavage Preference & \nEfficiency of {batch.upper()} Mutants"
        plot_combined_positions(
            pref_dfs, eff_series, pval_series,
            args.output_pdf,
            positions,
            title=title,
            vmin=pref_vmin, vmax=pref_vmax,
            eff_vmin=eff_vmin, eff_vmax=eff_vmax,
            cmap=args.cmap, annot=args.annot, fmt=args.fmt,
            square=args.square, fontsize=args.fontsize,
            cell_size=args.cell_size
        )
    else:
        # Single position
        position = args.position
        if position not in pref_dfs:
            raise RuntimeError(f"Position {position} not found in preference CSV")
        pref_df = pref_dfs[position]
        title = args.title or f"Base Preference at Position {position}"
        plot_combined(
            pref_df, eff_series, pval_series,
            args.output_pdf,
            title=title,
            position=position,
            vmin=pref_vmin, vmax=pref_vmax,
            eff_vmin=eff_vmin, eff_vmax=eff_vmax,
            cmap=args.cmap, annot=args.annot, fmt=args.fmt,
            square=args.square, fontsize=args.fontsize,
            cell_size=args.cell_size
        )
    print("Done!")
       
if __name__ == "__main__":
    main()