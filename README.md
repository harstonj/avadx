# AVA,Dx

## Pipeline
AVA,Dx (Analysis of Variation for Association with Disease) a computational method for defining the functional role of DNA variation in complex diseases. AVA,Dx uses exonic variants from whole exome or genome sequencing data to extract a *disease* signal and predict *disease* status.

This repository hosts the implementation of the improved AVA,Dx pipeline. An inital version version of the AVA,Dx pipeline and its in depth documentation can found [here](https://bitbucket.org/ywang202/avadx-meta)

## Models
Currently pre-computed models are available for:
- [**Crohns Disease**](https://bitbucket.org/bromberglab/avad-cd) (ava-cd)

# Requirements
The AVA,Dx pipeline is implemented and python and uses pre-built images runnable via docker or singularity. Thus, the only dependencies to run the pipeline are:

* [Docker](https://www.docker.com/get-started) or [Singularity](https://sylabs.io/singularity)
* Python ≥ 3.8

## Inputs
For the individuals (samples) in the sequencing panel/cohort
* A single Variant Calling File [VCF] comprising all samples. Default assembly is GRCh37/hg19
* Sample Information file [comma separated, CSV file], which contains a list of sample ids and corresponding class (e.g. disease vs healthy) labels. *optional* group/family or fold labels for cross validation can also be provided

## Outputs
* list of *(normalized) gene scores* representing functional effects of variants. Normalization is done by protein length
* list of label(e.g. disease)-relevance gene scores (per fold)
* model performance evaluation, which informs the value of N top ranked selected genes
* list of genes responsible for top performance across folds (ranked by number of folds affected)

# Usage

## Overview
Confirm you have the required versions of Docker/Singularity and Python installed. The following steps are detailed in the *Run AVA,Dx* section below:

1. Install AVA,Dx python package
2. Init pipeline - retrieves pre-built images & required databases. A configuration template file named *avadx.ini* is created in the working directory)
3. Edit configuration file *avadx.ini* created in step (2)
4. Run the pipeline

## Run AVA,Dx

### 1. Install AVA,Dx python package
```bash
pip install https://bitbucket.org/bromberglab/avadx/get/master.zip
```

### 2. Init pipeline
```bash
avadx --init
```

### 3. Edit configuration file
Open the file *avadx.ini* in the working directory and update the paths to the two required input files (**vcf file** and **samples info file**) and other parameters as required. **Note:** You can change the configuration file at runtime by supplying the path to the config file as argument (see `--help`).

### 4. Run the pipeline
If no path to an user configuration file is supplied, a default config file *avadx.ini* is expected in the current working directory

```bash
avadx --wd <PATH_TO_OUTPUT_DIRECTORY> --uid <MY_PIPELINE_RUN_ID> <PATH_TO_CONFIG_FILE>
```

For an extensive list of parameters use `--help`.

## Update AVA,Dx databases/datasources

### Update VM images
```bash
avadx --update vm
```

### Update Databases
```bash
avadx --update data
```

# Notes

* The pipeline defaults to hg19. For hg18, we recommend lifting over to hg19 first. For hg38, we recommend changeing the reference databases to hg38.
* This pipeline is currently for regular VCF file input (modifications are needed in order to use gVCF files).
* Manual interpretation of quality outliers, ethnicity, is a manual step strongly ancouraged.
* When the input VCF contains variants with no SNAP score records available, the *varidb* database needs to be updated.

# Dependencies
Following frameworks, tools and libraries are used in the pre-built images:

* Python
* Java
* R
  * BiocManager (SNPRelate) 
  * Packages (optparse, data.table, tidyverse, seqinr, stringr, EthSEQ, SNPRelate, e1071, caret, ggfortify, R.utils, PRROC, xlsx, ranger)
* [samtools/tabix](https://github.com/samtools/tabix)
* [bcftools](https://samtools.github.io/bcftools/)
* [ANNOVAR](http://annovar.openbioinformatics.org)
* [SNAP](https://bitbucket.org/bromberglab/snap)
