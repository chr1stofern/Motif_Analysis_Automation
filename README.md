# Title Here

## Overview
### Experiment
In vitro cleavage site mapping using random substrates

Samples involved (22 in total, 11 for each replicate):

Lbu: No Protein  
Lbu: Targetting crRNA 
Lwa: No Protein
Lwa: Targetting crRNA  
Psp: No Protein
Psp: Targetting crRNA  
bt3: No Protein  
bt3: Targetting crRNA  
Rfx: No Protein
Rfx: Targetting crRNA  
Dictionary: oKWW2182  

### What the Code Does
1. Raw NGS Reads Preprocessing
- adaptor trimming: use different 3' sequences for dictionary (dict) and cleavage (clv) sequences
- clv: trim the last 4 nt
- assign min and max length after cutadapt (MIN: 35 for dict, 1 for clv; MAX: both 35)
2. Motif Analysis
- both: draw a histogram to check the cutoff
- dict: Get 5' sequences with the selected length
- pandas: pd.merge. Use the length number to get the sequences near the cleavage position
- Calculate the mean of each 2/3nt motif
- Make a scatterplot to compare rep1 and 2 to see which motifs are more frequently cleaved
- generate plogo to visualise significant motifs

## Prerequisites
Note that this package is designed to be run on a Linux distribution as the kpLogo conda package is installed for linux-64 on the genomedk channel. You may also choose to install kpLogo with other methods, but this is not recommended. Alternatively, if you wish to obtain the logos manually using pLogo's [web server](https://plogo.uconn.edu/), you may skip "Setting Up kpLogo" in the "Installation Guide" and refer to the "How to Use" section.

We use conda to manage the installation of kpLogo. If you wish to opt for this route, please refer to [Conda's official installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html):
            

We also recommend using [uv](https://docs.astral.sh/uv/) for managing python-related enviroinments. It can be installed by simply running on the terminal:
````curl -LsSf https://astral.sh/uv/install.sh | sh```` 

## Installation Guide
Please clone the repository from GitHub to a desired location on your device or server:
````git clone <repository-url>````

Then, move into the cloned directory:
````cd <repository-folder>````

### Setting Up Conda Environment:
Run the following command to create an environment with the necessary dependencies, including kpLogo, via conda:
````conda env create -f environment.yaml --solver=libmamba````
(On older versions of Conda where libmamba may not be installed, install it by running the following first:
````conda install conda-forge::conda-libmamba-solver````
then create the environment.)

Activate the conda environment with:
````conda activate motif_analysis_env````

### Setting Up UV Environment:
Simply run:
````uv sync````

## Initial Configuration
### Configuring Input & Output Directories/Files:
Please modify the ```config.yaml``` accordingly and confirm the correct paths are provided properly. You may also change the sample names and parameters in ````samples.csv```` if necessary.

### Downloading & Naming Samples:
To download the raw NGS data necessary for reproducing figures in our publication, please run the bash script fetch_raw.sh with the syntax:
````./fetch_raw.sh <path/tosamples.csv> <path/to/raw_dir>````

(Be sure to add execution rights before running: chmod +x ./fetch_raw.sh)
The download process will considerable time. If the process crashed mid-run or some samples failed to download, simply run it again. Files that are already downloaded will be skipped. Please check that a total of 22 files of the extension ````.fastq.gz```` are downloaded.

Alternately, you may also download the raw data via the [SRA Explorer](https://sra-explorer.info/#) and search for the accession number GSE313507. Include all 22 files within the series and click "Add to collection". Click on "22 saved datasets" and run the automatically generated download script of your choice.

Note that the script supports 1-2 replicates. To process more than 2 replicates, you may either edit the internal code, or simply repeat the process and specify a different set of samples to run (save the previous outputs first before re-running).

### Raw Sequence Processing --- Parameters for Cleavage and Dictionary Type Samples:
Please make sure filenames of RAW NGS sequences follow the pattern: {replicate}_{protein}-{type}.fastq.gz, e.g., 1_Lbu-NP.fastq.gz. Note that only read 1 (forward read) is required for our data. 

Parameters, including adapter sequences and length, for each sample can be changed in samples.csv. Simply use the csv file provided for reproducing the original figures. You may change the parameters based on your on experiment.

### Scatter Plot Generation:
If you decide to apply the pipeline on sequences of your choice, please modify the sample importing and exporting logic of ````generate_plots.py```` accordingly, particularly for the ````MOTIF_FIGS```` dictionary, and the sections Supp-Fig-6a and Fig-2f.

### pLogo Generation with kpLogo:


### Output Files:
A total of 3 subdirectories, namely ````main_fig````, ````supp_fig```` and ````motif_plogo````, and one xlsx sheet will be produced. If you choose to use our dataset from GSE313507, there will be 7 files in ```main_fig``` by default, each corresponding to parts of the main figures in the original publication; 6 files in ````supp_fig````, each corresponding to parts of the supplementary figures in the original publication; and 10 files in ````motif_plogo```` corresponding to the motifs of the 5 Cas13 protein tested in two replicates.

## Running the Script
To run the entire pipeline, execute the following:
````uv run snakemake --core 30````

## Cleanup Guide

### Removing Intermediate Output & Output Files:
Simply remove the directories you specified in ````config.yaml```` by running ````rm -rf <dir>````. Alternatively, you may use the automatic clean up rule in snakemake by running:
````uv run snakemake --cores 4 clean_all```` 

### Removing Environments
For the Conda environment, deactivate it first:
````conda deactivate````

Then, remove the conda environment:
````conda env remove -n motif_analysis_env````

and clear the package cache:
````conda clean --packages --tarballs````

For the UV environment, deactivate with:
````deactivate````

Then, remove the environment with:
````rm -rf .venv````