# AVA,Dx pipeline

AVA,Dx (Analysis of Variation for Association with Disease) a computational method for defining the functional role of DNA variation in complex diseases. AVA,Dx uses exonic variants from whole exome or genome sequencing data to extract a *disease* signal and predict *disease* status.

This repository hosts the implementation of the AVA,Dx pipeline.

A inital version version of the AVA,Dx pipeline and its documentation can found [here](https://bitbucket.org/ywang202/avadx-meta).

## Requirements
* [Docker](https://www.docker.com/get-started) or [Singularity](https://sylabs.io/singularity)
* Python ≥ 3.8

### Required inputs
* Variant Calling File [VCF file]: Default assembly is GRCh37/hg19
* Sample Information file [CSV file]: Contains a list of sample ids to process and their class labels and *optional*: group labels, fold labels for cross validation

### Outputs
* a *gene score* table
* a **normed** *gene score* table, normalizing every *gene score* by protein length
* list with selected genes
* model performance evaluation

## First Steps

### Overview

1. Confirm you have the required versions of Docker/Singularity and Python installed
2. Install the AVA,Dx pipeline python package
3. Run the pipeline VM update to retrieve the latest pre-build images
4. Run the pipeline init command to retrieve and pre-process the required databases
5. Copy the pipeline.ini configuration file template and change the paths to the two required input files. Default parameters can be adusted here as well.
6. Run the pipeline by providing at minimum the path to the user configuration file.

### Install AVA,Dx python package

### Update VM images
```bash
avadx-meta --update vm
```
### Update Databases
```bash
avadx-meta --update data
```
### Create pipeline.ini config from template
Copy the provided configuration file template (`pipeline.ini`) and update the paths to the two input files (vcf file and samples info file); change any default parameters as required.
### Run the pipeline
```bash
avadx-meta --wd <PATH_TO_OUTPUT_DIRECTORY> --uid <MY_PIPELINE_RUN_ID> <PATH_TO_CONFIG_FILE>
```
## Notes

* The pipeline defaults to hg19. For hg18, we recommend lifting over to hg19 first. For hg38, we recommend changeing the reference databases to hg38.
* This pipeline is currently for regular VCF file input (modifications are needed in order to use gVCF files).
* Manual interpretation of quality outliers, ethnicity, is a manual step strongly ancouraged.
* When the input VCF contains variants with no SNAP score records available, the *varidb* database needs to be updated.

## Dependencies
Following tools and libraries are used in the workflow:

* Python
* Java
* R
  * BiocManager (SNPRelate) 
  * Packages (optparse, data.table, tidyverse, seqinr, stringr, EthSEQ, SNPRelate, e1071, caret, ggfortify, R.utils, PRROC, xlsx, ranger)
* [samtools/tabix](https://github.com/samtools/tabix)
* [bcftools](https://samtools.github.io/bcftools/)
* [ANNOVAR](http://annovar.openbioinformatics.org)
* SNAP
