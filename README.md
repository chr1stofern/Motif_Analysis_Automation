# Title Here

## Overview
### Experiment
In our previous publication, we mapped the in vitro cleavage sites of various Cas13 proteins using random substrates. In total, 22 samples/sequencing files (11 for each of two replicate) were obtained for downstream analysis. The samples for each replicate are as follows:

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
### Conda
Note that this package is designed to be run on a Linux distribution as the kpLogo conda package is installed for linux-64 on the genomedk channel. You may also choose to install kpLogo with other methods, but this is not recommended. Alternatively, if you wish to obtain the logos manually using pLogo's [web server](https://plogo.uconn.edu/), you may refer to the "pLogo Generation with pLogo Web Page" section.

We use conda to manage the installation of kpLogo and other packages. Please refer to [Conda's official installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html) if you do not have Conda installed.
            
### UV
We also recommend using [UV](https://docs.astral.sh/uv/) for managing python-related enviroinments. It can be installed by simply running on the terminal:
````curl -LsSf https://astral.sh/uv/install.sh | sh```` 

## Installation Guide
Please clone the repository from GitHub to a desired location on your device or server:
````git clone <repository-url>````.

### Setting Up Conda Environment:
Please move into the cloned directory before proceeding::
````cd <repository-folder>````.

Run the following command to create an environment with the necessary dependencies as declared by environment.yaml, including kpLogo, via conda:
````conda env create -f environment.yaml --solver=libmamba````.
(On older versions of Conda where libmamba may not be installed, install it by running the following first:
````conda install conda-forge::conda-libmamba-solver````
then create the environment.)

Activate the conda environment with:
````conda activate motif_analysis_env````. Please run this everytime before running the script.

### Setting Up UV Environment:
Simply run:
````uv sync````.

## Initial Configuration
### Configuring Input & Output Directories/Files:
Please modify the ```config.yaml``` accordingly and confirm the correct paths are provided properly. 

### Setting Up Samples:
This package can be used in one of 3 ways: with example (toy) sequences, with raw sequences from which we produced the figures for our previous publication, and with your own data.

The example sequences are already included in the repository, and are obtained by running the command ````awk 'NR % 10 == 0' {ORIGINAL_SEQ.collapsed.tsv} > {EXAMPLE_SEQ.collapsed.tsv}```` on RfxCas13d and dictionary sequences.

To download the raw NGS data (deposited in the NCBI Gene Expression Omnibus, GEO, under the accession number GSE313507) necessary for reproducing figures in our publication, please add execution rights to the download script by running: ````chmod +x ./fetch_raw.sh````, then run the bash script fetch_raw.sh with the syntax:
````./fetch_raw.sh {path/tosamples.csv} {path/to/raw_dir}````

The download process will considerable time. If the process crashed mid-run or some samples failed to download, simply run it again. Files that are already downloaded will be skipped. Please check that a total of 22 files of the extension ````.fastq.gz```` are downloaded.

Alternately, you may also download the raw data via the [SRA Explorer](https://sra-explorer.info/#) and search for the accession number GSE313507. Include all 22 files within the series and click "Add to collection". Click on "22 saved datasets" and run the automatically generated download script of your choice. A total of 44 sequence files will be downloaded, since each sample produces two reads (read 1/forward read, and read 2/reverse read). Note that only read 1 (forward read) is required for our analysis. 

Note that the script supports 1-2 replicates. To process more than 2 replicates, you may either edit the internal code, or simply repeat the process and specify a different set of samples for the script to run (save the previous outputs first before re-running).

### Naming Samples
If you choose to use the toy sequences or the .sh download script, the downloaded files will already be in the correct naming format. If you choose to use your own data, please make sure your raw sequence files match the following naming format:
````{rep#}_{sample/protein name}-{type: CR (crRNA and protein added), or NP (no protein), or DT (dictionary)}.fastq.gz````. For example, for a crRNA added PspCas13b conducted in replicate 1, it shall be named 1_Psp-CR.fastq.gz ({sample/protein name} can be changed to any string of your liking, but must not include "-" or "_"); for a dictionary conducted in replicate 2, it shall be named 2_Dict-DT.fastq.gz (again, "Dict" can be changed to other strings except "-" and "_").

The full list of sequencing samples should also be documented in the sample sheet (by default ````samples.csv```` or ````examples.csv````, you can change where ````sample_sheet```` points to in ````config.yaml````). 
- Column 1 ````accession````: specifies the GEO assession number of each sample
- Column 2 ````replicate````: specifies the replicate number --- ````1````, or ````2````
- Column 3 ````protein````: specifies the protein or sample name (any string without "-" or "_")
- Column 4 ````type````: specifies the sample type --- ````CR```` (crRNA and protein added), or ````NP```` (no protein), or ````DT```` (dictionary)
- The remaining columns are described in the next section

If you are using our original data, simply point ````sample_sheet```` to ````samples.csv```` and need not change its content. Similarly, if you are using the example data, point ````sample_sheet```` to ````examples.csv```` and need not its content.

### Configuring Parameters for Raw Sequence Preprocessing:
Please make sure filenames of RAW NGS sequences follow the naming pattern as aforementioned. After that, please obtain the appropriate 5' and 3' adaptor sequences used in your experiment's sequencing setup, and paste them into column 5 and 6 of the sample sheet. 

Parameters, including adapter sequences and length, for each sample can be changed in samples.csv. Simply use the csv file provided for reproducing the original figures. You may change the parameters based on your on experiment.

### Scatter Plot Generation:
If you decide to apply the pipeline on sequences of your choice, please modify the sample importing and exporting logic of ````generate_plots.py```` accordingly, particularly for the ````MOTIF_FIGS```` dictionary, and the sections Supp-Fig-6a and Fig-2f.

### pLogo Generation with kpLogo Package:

### pLogo Generation with pLogo Web Page:


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