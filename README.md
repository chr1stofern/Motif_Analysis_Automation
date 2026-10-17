# Title Here

## Overview

## Prerequisites
Note that this package is designed to be run on a Linux distribution as the kplogo conda package is only available for linux-64 on the genomedk channel. You may also choose to install kplogo with other methods, but this is not recommended. Alternatively, if you wish to obtain the logos manually using pLogo's web server, you may skip "Setting Up kpLogo" in the "Installation Guide" and refer to the "How to Use" section.

We use conda to manage the installation of kpLogo. If you wish to opt for this route, you may install Conda by running:
conda install anaconda::conda
                          
## Installation Guide
Clone the repo from github

### Setting Up kpLogo:
Run the following command to install kplogo and its dependencies via conda:
conda env create -f environment.yaml

Then, activate the conda environment with:
conda activate kplogo_env                   

## How to Use


## Cleanup Guide
First, deactivate the conda environment:
conda deactivate.

Then, remove the conda environment:
conda env remove -n kplogo_env,

and clear the package cache:
conda clean --packages --tarballs.