# High-Throughput in vitro Cleavage Assay for Cas13 Motif Discovery

## Overview
### Experiment
In our previous publication, we mapped the in vitro cleavage sites of various Cas13 proteins using random substrates. In total, 22 samples/sequencing files (11 for each of two replicate) were obtained for downstream analysis. The samples for each replicate are as follows:

Lbu: No Protein\
Lbu: Targetting crRNA\
Lwa: No Protein\
Lwa: Targetting crRNA\
Psp: No Protein\
Psp: Targetting crRNA\
bt3: No Protein\
bt3: Targetting crRNA\
Rfx: No Protein\
Rfx: Targetting crRNA\
Dictionary: oKWW2182

### What the Code Does
1. **Raw NGS Reads Preprocessing:**\
Each sequencing read is reduced to the insert/substrate fragment it contains, i.e. the sequence the Cas13 protein actually recognises. The preprocessing is conducted with various command line packages, including cutadapt, fastq_quality_filter, fastx_collapser and seqkit, all of which are called by a fully automated `snakemake` script (`Snakefile`).
    - **Adapter trimming:** the 5' and 3' sequencing adapters are removed from both cleavage and dictionary samples.
    - **Trimming the last 4nt random region:** for cleavage libraries, the last 4 nt of the insert were synthesised as a random region and carry no motif information and shall be trimmed. Dictionary inserts are fixed and thus unmodified.
    - **Length filtering:** reads are kept only if their insert falls within the expected size range: 1–35 nt for cleavage libraries, and exactly 35 nt for the dictionary. This removes empty inserts, adapter dimers, and reads longer than the library design.
    - **Quality filtering:** reads are discarded unless at least 90% of their bases have a Phred quality score of 30 or higher. This removes reads with sequencing errors that would otherwise corrupt the insert.
    - **Collapsing reads:** identical inserts are collapsed into unique sequences with their read counts, and the result is written as a tab-separated table (*.collapsed.tsv). This table will be read by `generate_plots.py` which conducts the downstream motif analysis.

2. **Motif Analysis:**\
The collapsed reads are filtered for a final time. Then, the inserts are mapped back to the dictionary to measure how efficiently each sequence was cleaved, and the enriched motifs are then visualised.
    - **Dictionary filtering:** dictionary sequences are kept only if they were observed above a minimum read count (100 in replicate 1, 50 in replicate 2).
    - **Insert filtering:** cleavage fragments are kept only if they were observed above a minimum read count of 20 and have an insert length between 15 and 24 nt.
    - **Mapping to dictionary:** each cleavage fragment is matched to the dictionary prefix of the same length. The position of the fragment's 3' end in the dictionary is taken as the cleavage site, and the 2-nt and 3-nt motifs around that position are extracted.
    - **Computing relative cleavage efficiency:** for each fragment, the ratio of its count in the cleavage library to its count in the dictionary gives a relative cleavage efficiency, which is used to rank all fragments.
    - **Replicate comparison:** the mean efficiency of each motif is compared between the two replicates and plotted as a scatterplot. Motifs that are consistently enriched or depleted across replicates can then be identified.
    - **Logo generation:** the top-ranked motifs are written as input sequences for pLogo, which visualises the significantly enriched bases at the cleavage site.

### Intermediate Output Files
The intermediate output directory is specified as `intermediate_dir` in the configuration file (more details in [Configuring Input & Output Directories/Files](#configuring-input--output-directoriesfiles)). The `intermediate_dir` is utilised by the script for depositing three products, 1., partially or fully processed sequencing reads (but not RAW reads) from `Snakefile`; 2., filtered and mapped sequences from `generate_plots.py` that serve as inputs for kpLogo motif visualisation; and 3, the full output of the kpLogo package for every sample. The first is stored in `intermediate_dir/preprocessed_seq_dir`, the second is stored in `intermediate_dir/logo_input`, and the third is stored in `intermediate_dir/logo_output`.

Note that the processed example reads are stored in `intermediate_files/sampled_tsv` inside the cloned repository. Refer to [Initial Configuration](#initial-configuration) on how to use the example sequences and configure the directory paths.

### Final Output Files:
A total of 3 subdirectories, namely `main_fig`, `supp_fig` and `motif_plogo`, and one xlsx sheet will be produced in the `output_dir`. If you choose to use our dataset from GSE313507, there will be 5 files in `main_fig` by default, each corresponding to parts of the main figures in the publication; 10 files in `supp_fig`, each corresponding to parts of the supplementary figures in the publication; and 10 files in `motif_plogo` corresponding to the motifs of the 5 Cas13 protein tested in two replicates. For details on choosing parameters, please refer to the [Initial Configuration](#initial-configuration) section

By default, the following files with be produced in each subdirectory:
1. **`supp_fig`:**
    - **Supp. fig. 3a, b:** sequencing depth of dictionary libraries, plotted as the distribution of the number of reads per unique sequence in each replicate. Fig. 3a corresponds to replicate 1, fig. 3b corresponds to replicate 2. The dashed blue line indicates the minimum read threshold used for filtering the dictionary reads before generating the boxplot.
    - **Supp. fig. 4a, b:** sequencing depth of cleavage libraries, plotted as the distribution of the number of reads per unique sequence each replicate. Fig. 4a corresponds to replicate 1, fig. 4b corresponds to replicate 2. The dashed blue line indicates the minimum read threshold used for filtering the cleavage reads before generating the boxplot. Note that the plots correspond to the subset of reads after mapping onto the dictionary.
    - **Supp. fig. 5a, b:** distribution of relative cleavage efficiencies of cleaved sequences, plotted against the relative cleavage efficiency rank. Fig. 5a corresponds to replicate 1, fig. 5b corresponds to replicate 2. The dashed blue line indicates the rank threshold, used for filtering in only the top sequences to the left of the line for downstream pLogo motif visualisation.
    - **Supp. fig. 6a:** relative cleavage efficiency of various subsets among all ranked motifs from cleaved samples.

Note that for supp. fig. 6a, only one pair of duo-replicate plot for the protein PspCas13b is shown in the original publication. Here, the script automatically generates one paired plot (or single plot if only 1 replicate is provided) for each protein sample involved. Thus, when using the full dataset, a total of of 5 plots will be generated, each corresponding to one of five protein samples.

2. **`mapped_clv_seq_full.xlsx`:**
A full record of every cleavage fragment that was successfully mapped to the dictionary. Each sample is recorded on a separate sheet. Data obtained from the analysis, including 2/3nt motifs and cleavage efficiency ranks are also included along with the fragment sequence itself.

3. **`main_fig`:**
    - **Fig. 2c:** relative cleavage efficiency of dinucleotide motifs in both replicates of non-PspCas13b proteins. Most efficient (top few) motifs are shown in red, inefficient motifs are shown in gray
    - **Fig. 2e:** relative cleavage efficiency of trinucleotide motifs in both replicates of PspCas13b. Most efficient (top few) motifs are shown in red, inefficient motifs are shown in gray
    - **Fig. 2f:** Mann–Whitney U test of relative cleavage efficiency of representative trinucleotide motifs of PspCas13b. Note that the significant values are exported in `fig2_psp_mann-whitney.tsv` by default.

4. **`motif_plogo`:**
This directory contains a minimal version of kpLogo's output, i.e., the .pdf file of the pLogo generated for each CR-NP sample pair. The full output generated by kpLogo can be found in `intermediate_dir/logo_output` as mentioned in [Intermediate Output Files](#intermediate-output-files).

Note that the file names of the main and supplementary figures can be modified. This can be done via the `out:` and `out_template:` parameters in config.yaml.

## Prerequisites
### Conda
Note that this package is designed to be run on a Linux distribution as the kpLogo conda package is installed for linux-64 on the genomedk channel. You may also choose to install kpLogo with other methods, but this is not recommended. Alternatively, if you wish to obtain the logos manually using [pLogo's web server](https://plogo.uconn.edu/), you may refer to the [pLogo Generation with pLogo Web Page](#plogo-generation-with-plogo-web-page) section.

We use conda to manage the installation of kpLogo and other packages. Please refer to [Conda's official installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html) if you do not have Conda installed.
            
### UV
We also recommend using [UV](https://docs.astral.sh/uv/) for managing python-related enviroinments. It can be installed by simply running on the terminal:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Installation Guide
Please clone the repository from GitHub to a desired location on your device or server:
```bash
git clone <repository-url>
```

### Setting Up Conda Environment:
Please move into the cloned directory before proceeding:
```bash
cd <repository-folder>
```

On older versions of Conda where libmamba may not be installed, install it by running the following:
```bash
conda install conda-forge::conda-libmamba-solver
```
Run the following command to create an environment with the necessary dependencies as declared by environment.yaml, including kpLogo, via conda:
```bash
conda env create -f environment.yaml --solver=libmamba
```

Please activate the conda environment every time before running the script:
```bash
conda activate motif_analysis_env
```

### Setting Up UV Environment:
Simply run:
```bash
uv sync
```

## Initial Configuration
### Configuring Input & Output Directories/Files:
If you are testing the example subset or the full dataset (more in the next section), need not modify any files. Simply move on to the next section. However, for the full dataset, please Please make sure the cloned directory has sufficient storage (at least 100 GB).

If you are using your own data, please locate the configuration file `config.yaml` and modify the parameters `input_dir`, `intermediate_dir` and `output_dir` to point to directories of your liking. Again, make sure these locations have sufficient storage. You may also modify the names of `preprocessed_seq_dir` and `output_dir` if you wish. But they will still be produced under the same parent directory. Note that both relative paths (without "./", see the file for examples) and full paths are allowed, although full paths are recommended.

### Setting Up Samples:
This package can be used in one of 3 ways: with an example subset, with the full dataset from which we produced the figures for our previous publication, and with your own data.

The example sequences are included as a compressed .zip file under <program_repo>/intermediate_files/preprocessed_sequences.zip. These sequences make up a subset of the preprocessed intermediates of the full dataset. In particular, the preprocessed `.collapsed.tsv` of LwaCas13a cleavage samples and the dictionary samples for both replicates. Before proceeding, unzip the file in the same directory (intermediate_files) and check that the folder produced is named as `preprocessed_example_sequences`. Do not place the folder elsewhere for this use case.

To download the full dataset (deposited in the NCBI Gene Expression Omnibus, GEO, under the accession number GSE313507) necessary for reproducing figures in our publication, please add execution rights to the download script by running: `chmod +x ./fetch_raw.sh`, then run the bash script fetch_raw.sh with the syntax:
```bash
./fetch_raw.sh <full/path/to/samples.csv> <full/path/to/input_dir>
```

The download process will take considerable time. If the process crashed mid-run or some samples failed to download, simply run it again. Files that are already downloaded will be skipped. Please check that a total of 22 files of the extension `.fastq.gz` are downloaded.

Alternately, you may also download the raw data via the [SRA Explorer](https://sra-explorer.info/#) and search for the accession number GSE313507. Include all 22 files within the series and click "Add to collection". Click on "22 saved datasets" and run the automatically generated download script of your choice. A total of 44 sequence files will be downloaded, since each sample produces two reads (read 1/forward read, and read 2/reverse read). Note that only read 1 (forward read) is required for our analysis. 

Note that the script supports 1-2 replicates. To process more than 2 replicates, you may either edit the internal code, or simply repeat the process and specify a different set of samples for the script to run (save the previous outputs first before re-running).

### Naming Samples
If you choose to use the example sequences or the .sh download script, the downloaded files will already be in the correct naming format, and you may skip this and the next two section. For such use cases, you may skip the next section few sections on configuring parameters. The sample information sheet and config.yaml file are already ready.

If you choose to use your own data, please make sure your raw sequence files match the following naming format:
`{rep#}_{sample/protein name}-{type: CR (crRNA and protein added), or NP (no protein), or DT (dictionary)}.fastq.gz`. For example, for a crRNA added PspCas13b conducted in replicate 1, it shall be named 1_Psp-CR.fastq.gz ({sample/protein name} can be changed to any string of your liking, but must not include "-" or "\_"); for a dictionary conducted in replicate 2, it shall be named 2_Dict-DT.fastq.gz (again, "Dict" can be changed to other strings except "-" and "\_").

The full list of sequencing samples should also be documented in the sample information sheet. A template `my_samples.csv` is provided in `<program_repo>/my_config/`.
- Column 1 `accession`: specifies the GEO assession number of each sample. This column can be left empty if samples are not downloaded from GEO and are already in your local storage.
- Column 2 `replicate`: specifies the replicate number --- `1`, or `2`.
- Column 3 `protein`: specifies the protein or sample name (any string without "-" or "_").
- Column 4 `type`: specifies the sample type --- `CR` (crRNA and protein added), or `NP` (no protein), or `DT` (dictionary).
- The remaining columns are described in the next section.

### Configuring Parameters for Raw Sequence Preprocessing:
Please make sure filenames of your raw NGS sequences follow the naming pattern as aforementioned. After that, obtain the appropriate 5' and 3' adaptor sequences used in your experiment's sequencing setup, and paste them into column 5 and 6 respectively of the sample sheet. 

The remaining 3 columns of `sample_sheet` allows for additional control on the trimming of sequences after the initial adaptor sequence is removed.
- Column 7 `random_at_3p`: trims a number of bases from one end of the sequence, e.g., if set to 3, 3 bases will be removed from the 5' end of a sequence, if set to -4, 4 bases will be removed from the 3' end of a sequence, if set to 0, nothing will be removed.
- Column 8 `min_len`: the minimum length for the sequence to not be filtered out (i.e. any sequence with length < `min_len` will be discarded).
- Column 9 `max_len`:the maximum length for the sequence to not be filtered out (i.e. any sequence with length > `max_len` will be discarded).

### Dictionary Mapping & Scatter Plot Generation:
After the adaptor trimming, the adaptor-trimmed sequences are further processed with a count and length filter. The relevant parameters can be found and modified in the template `<program_repo>/my_config/my_config.yaml`.

**Count Filters:**
- `dict_min_count_rep1`: minimum count of a particular sequence in the dictionary sample of replicate 1 for it to be filtered in
- `dict_min_count_rep2`: minimum count of a particular sequence in the dictionary sample of replicate 2 for it to not be filtered in
- `clv_min_count_rep1`: minimum count of a particular sequence in the cleavage sample of replicate 1 for it to be filtered in
- `clv_min_count_rep2`: minimum count of a particular sequence in the cleavage sample of replicate 1 for it to be filtered in

Note that we use different values for `dict_min_count_rep1` and `dict_min_count_rep2`, since the sequencing depth of the dictionary from replicate 1 is higher than that of replicate 2. For processing our example dataset, or full dataset, keep the parameters as is. For processing your own data, please adjust them accordingly.

**Length Filters:**
- `min_insert_len`: minimum length of a particular sequence in the cleavage sample for it to be filtered in
- `max_insert_len`: maximum length of a particular sequence in the cleavage sample for it to be filtered in\

The final processed sequence (insert) to be mapped to the dictionary must be equal to or longer than the `min_insert_len`, and shorter than the `max_insert_len` before proceeding.
For processing our example dataset, or full dataset, keep the parameters as is. For processing your own data, please adjust them accordingly.

**Rank Filter:**
After all sequences are mapped to the dictionary, the sequences are ranked based on their relative cleavage efficiency, and only the top 4000 sequences are taken for pLogo visualisation for our full data. This can be controlled by the parameter `top_n_seq` in `my_config.yaml`, and it also changes the line drawn in the plot of relative cleavage efficiency against their ranking (supp. fig. 5) for ease of reference. Please adjust the value accordingly for your own data. You may first run the program with an arbitrary value, then inspect the resulting Supplementary Figure 5 to see where efficiency starts to level off. Choose a value near this point. 

**Parameters for Supp. Fig. 6a:**
The grouping of ranked motifs into separate bins can be modified to your liking. In the publication, we separated the motifs into the top 0-400 (top 400) ranked motifs, the top 401-1600, the top 1601-4000, and the remaining (4001-infinity). The actual grouping of data can be modified at `bins:` in the `supp6a:` section of `my_config.yaml`. For the full data, we write `bins:   [0, 400, 1600, 4000, .inf]` where .inf denotes infinity. To change the labels along the x-axis of the figure, modify the second parameter `labels:`. To change the colour of individual boxes for each group, modify `colors:`.

**Parameters for Fig. 2c,e:**
Fig. 2c,e plots the cleavage efficiencies of both replicates 1 and 2. If you are using your own data, please change the first indented line ({protein}-{type}:) to your sample names following the specified format (replicate number is not necessary). You may refer to the pre-built `example_config.yaml` and `full_config.yaml` as examples. 

Next, specify whether the protein recognises a 2nt or 3nt motif in `motif_type:` as `2nt_motif` or `3nt_motif` respectively. 

n_top specifies the number of top ranked motif(s) you wish to highlight on the plot, and conversely, n_bottom specifies the number of the lowest ranked motif(s) you wish to highlight on the plot. If there are other motifs you wish to highlight but are not among the top and bottom few motifs, you can specify them with the format `highlight: ["motif_string_1", "motif_string_2"]`, e.g., `highlight: ["CT", "GG"]`. 

`xlim` and `ylim` controls the visible range of the two axes. Modify them accordingly if some data points clip outside of the figure. You may do so by checking the range of the data to be plotted, which is printed in the "Producing Fig-2c,e" section during run-time, then increasing the parameters to cover the entire range.

**Parameters for Fig. 2f:**
This figure is only generated for PspCas13b when using the full dataset. If you wish to construct this plot for your own data, please modify `target:` to select the desired sample for analysis. 

To select other motifs relevant for your protein, please modify `motifs:`. The `reference:` motif specifies the motif that the target `motifs:` will be tested against for significance. 

`rep1_ylim:` and `rep2_ylim:` controls the height of the y-axis of the plot for replicate 1 and replicate 2 respectively. A .tsv file containing the significance values (obtained from the Mann–Whitney U test of the motifs) will also be exported. You may refer to the relevant part in the pre-built “example_config.yaml” and “full_config.yaml” as examples.

### pLogo Generation with kpLogo Package:
The generation of motif logos for all samples is automated with the kpLogo package. You may simply run the script without changing any parameters for kpLogo for all data use cases. If you do wish to modify certain parameters, you may do so at the kpLogo block in `config.yaml`. For detailed instructions on using kpLogo, please consult the package manual by running `kpLogo` on the terminal (while the [Conda environment](#setting-up-conda-environment) is activated). 

```yaml
max_k: 1
min_count: 0
p_cutoff: 0.05
fdr: false
region: null
weighted: false
ranked: false
extra_args: "-plot s -stack_order -1 -fix 1.1 -pseudo 2.5"
```
Please keep the parameter `-pseudo 2.5` unmodified. This adds a "pseudocount" (default = 1.0) during the kpLogo's computation of p-values and prevents it from being 0. Since the NGS sequences contains an extremely large number of sample (reads), the pseudocount is necessary to compensate for the highly deflated p-values.

### pLogo Generation with pLogo Web Page:
First open [pLogo generation website](https://plogo.uconn.edu/) on a browser of your choice. Keep the option `remove duplicate sequences` checked for both foreground (left box with the label "1. Paste or upload aligned sequences" on top) and background (the right box with the label "2. Select, upload, or paste a background"). Then, drag and drop the processed and mapped sequences into the appropriate input boxes, or upload them by clicking `choose file`. For example, 1_RfxExample-CR.collapsed.tsv (stored in `intermediate_files/logo_input`) shall be used as the foreground, while 1_RfxExample-NP.collapsed.tsv shall be used as the background. When both foreground and background sequences are approved, click the `generate pLogo` button. You may directly export the generated logo by clicking `Export Image` on the `export` tab. Or, you may also modify other customisation and export file parameters before doing so. If you wish to obtain the raw data of the plot, you may also obtain it by clicking the `Export raw position weight matrix values` button on the `export tab`. For other samples/replicates, repeat the process with the relevant foreground and background files.

## Running the Script
To run the entire pipeline for the example dataset, execute the following:
````bash
uv run snakemake --configfile ./example_config/example_config.yaml --cores 30

````

To run the entire pipeline for the full dataset, execute the following:
````bash
uv run snakemake --configfile ./full_config/full_config.yaml --cores 30 

````

To run the entire pipeline for your own dataset, execute the following:
````bash
uv run snakemake --configfile ./my_config/my_config.yaml --cores 30 

````

## Cleanup Guide
### Removing Intermediate Output & Output Files:
Simply remove the directories you specified in `config.yaml` by running `rm -rf <dir>`. Alternatively, you may use the automatic clean up rule in snakemake by running:
```bash
uv run snakemake --cores 4 clean_all
``` 

### Removing Environments
For the Conda environment, deactivate it first:
```bash
conda deactivate
```

Then, remove the conda environment:
````conda env remove -n motif_analysis_env````

and clear the package cache:
```bash
conda clean --packages --tarballs
```

For the UV environment, deactivate with:
```bash
deactivate
```

Then, remove the environment with:
```bash
rm -rf .venv
```