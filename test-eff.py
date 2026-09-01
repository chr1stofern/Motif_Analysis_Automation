import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import dunnett

def process_fluorescence_file(file_path, output_dir):
    print(f"Processing: {file_path}")
    
    # 1. Read input CSV
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    file_name = os.path.basename(file_path)
    batch_name = file_name.split('_')[1].split('.')[0]
    # Derive WT name (adjust if your WT label is different)
    wt_name = batch_name.capitalize() + batch_name[0].capitalize()
    print(f"WT identifier: {wt_name}")

    # Clean up columns
    rep_cols = ['Rep1', 'Rep2', 'Rep3']
    for col in rep_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=rep_cols)

    # ============================================================
    # Part A: Summary statistics from raw (original) data
    # ============================================================
    df['Mean_Fluorescence'] = df[rep_cols].mean(axis=1)
    df['SD_Fluorescence'] = df[rep_cols].std(axis=1)
    wt_mean_raw = df[df['Mutants'] == wt_name]['Mean_Fluorescence'].values[0]
    df['Inverted_Log2_FC'] = -np.log2(df['Mean_Fluorescence'] / wt_mean_raw)
    df['Knockdown_Percentage'] = (1 - (df['Mean_Fluorescence'] / wt_mean_raw)) * 100

    # ============================================================
    # Part B: Normalise each replicate column by its median
    # ============================================================
    medians = df[rep_cols].median(axis=0)
    norm_df = df[rep_cols].div(medians, axis=1)
    norm_df['Mutants'] = df['Mutants']

    # Compute the mean of the normalised replicates for each mutant
    # Step 1: group by Mutants and compute the mean of each replicate column
    norm_means_df = norm_df.groupby('Mutants')[rep_cols].mean()
    # Step 2: average across the three replicates to get a single mean per mutant
    norm_means = norm_means_df.mean(axis=1)   # now a Series with index = Mutants

    # WT mean from normalised data
    wt_norm_mean = norm_means[wt_name]

    # Compute normalised Inverted_Log2_FC (using the normalised means)
    df['Inverted_Log2_FC_Norm'] = -np.log2(df['Mutants'].map(norm_means) / wt_norm_mean)

    # ============================================================
    # Part C: Statistical test on per‑replicate log₂ fold‑changes
    #         (from normalised data)
    # ============================================================
    # Melt normalised data to long format
    norm_long = norm_df.melt(id_vars=['Mutants'], value_vars=rep_cols,
                             var_name='Replicate', value_name='Norm_Fluorescence')
    # Per-replicate log2 FC
    norm_long['Log2_FC'] = np.log2(norm_long['Norm_Fluorescence'] / wt_norm_mean)

    # Group by mutant and collect the Log2_FC values
    mutant_groups = norm_long.groupby('Mutants')['Log2_FC'].apply(list).to_dict()

    control_data = mutant_groups[wt_name]
    alternative_mutants = [m for m in mutant_groups.keys() if m != wt_name]
    alternative_data = [mutant_groups[m] for m in alternative_mutants]

    p_values_dict = {wt_name: 1.0}
    if alternative_data:
        try:
            res = dunnett(*alternative_data, control=control_data, alternative='two-sided')
            for idx, mutant in enumerate(alternative_mutants):
                p_val = res.pvalue[idx]
                p_values_dict[mutant] = max(0.0, min(p_val, 1.0))
        except Exception as e:
            print(f"Warning: Dunnett's test failed: {e}")
            for mutant in alternative_mutants:
                p_values_dict[mutant] = np.nan

    df['Dunnett_P_Value'] = df['Mutants'].map(p_values_dict)
    df['Significant_0.05'] = df['Dunnett_P_Value'] < 0.05

    # ============================================================
    # Part D: Save output with both raw and normalised FC
    # ============================================================
    base_name = os.path.basename(file_path)
    output_name = base_name.replace("fluor_", "eff_")
    output_path = os.path.join(output_dir, output_name)

    ordered_cols = [
        'Mutants', 'Rep1', 'Rep2', 'Rep3',
        'Mean_Fluorescence', 'SD_Fluorescence',
        'Inverted_Log2_FC',      # raw effect size
        'Inverted_Log2_FC_Norm', # normalised effect size
        'Knockdown_Percentage',
        'Dunnett_P_Value', 'Significant_0.05'
    ]
    df[ordered_cols].to_csv(output_path, index=False)
    print(f"Saved compiled data file to: {output_path}")

def main():
    input_dir = 'fluor'
    output_dir = 'eff'

    if not os.path.exists(input_dir):
        print(f"Error: Target directory '{input_dir}' not found.")
        return

    os.makedirs(output_dir, exist_ok=True)
    csv_files = glob.glob(os.path.join(input_dir, 'fluor_*.csv'))

    if not csv_files:
        print(f"No source files found.")
        return

    print(f"Found {len(csv_files)} batch data file(s).")
    for file_path in csv_files:
        process_fluorescence_file(file_path, output_dir)
    print("\nAll operations completed successfully!")

if __name__ == '__main__':
    main()