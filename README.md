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
4. Edit sample info file
5. Run the pipeline

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

### 4. Edit sample info file
Create the samples info file specified in the *avadx.ini* configuration file as *samples = <PATH_TO_FILE>*. Minimum information required is **sample id** and **class label**, **group label** and manual **fold** assignment for cross-validation are oprional.

**Format:** File has to be in CSV (comma separated) format and include a header with column names.

| sampleid |class | *[group]* | *[fold]* |
| -------- | ---- | --------- | -------- |
| sample_1 |  0   |    *1*    |    *1*   |
| sample_2 |  0   |    *1*    |    *2*   |
| sample_3 |  0   |    *2*    |    *4*   |
| sample_4 |  1   |    *2*    |    *1*   |
| sample_5 |  1   |    *3*    |    *2*   |
| sample_6 |  1   |    *3*    |    *3*   |

### 5. Run the pipeline
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

### Change scoring function

AVA,Dx provides a plugin based, flexible way to alter variant scoring and gene scoring aggregation functions.
In the configuration (.ini) file simply specify the python to a python file containing a single scoring function.
Templates for both scoring functions can be found in this repository under [python/scoring_functions](https://bitbucket.org/bromberglab/avadx/src/master/python/scoring_functions/).
By default, AVA,Dx only scores *SNP* variants. TO also include *INDELS* for scoring update the config file accordingly and specify an approproate scoring in the varian scoring function.
Refer to the Annovar variant types table below for more details for scorable variant types.

### Annovar variation types available in the variant scoring function (key = "type")
| Annotation                       | Precedence | Explanation                                                                                                                                                                                                                                                                                         | Sequence Ontology                  |
|----------------------------------|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------|
| frameshift insertion             | 1          | an insertion of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                                    | frameshift_elongation (SO:0001909) |
| frameshift deletion              | 2          | a deletion of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                                      | frameshift_truncation (SO:0001910) |
| frameshift block substitution    | 3          | a block substitution of one or more nucleotides that cause frameshift changes in protein coding sequence                                                                                                                                                                                            | frameshift_variant (SO:0001589)    |
| stopgain                         | 4          | a nonsynonymous SNV; frameshift insertion/deletion; nonframeshift insertion/deletion or block substitution that lead to the immediate creation of stop codon at the variant site. For frameshift mutations; the creation of stop codon downstream of the variant will not be counted as "stopgain"! | stop_gained (SO:0001587)           |
| stoploss                         | 5          | a nonsynonymous SNV; frameshift insertion/deletion; nonframeshift insertion/deletion or block substitution that lead to the immediate elimination of stop codon at the variant site                                                                                                                 | stop_lost (SO:0001578)             |
| nonframeshift insertion          | 6          | an insertion of 3 or multiples of 3 nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                     | inframe_insertion (SO:0001821)     |
| nonframeshift deletion           | 7          | a deletion of 3 or mutliples of 3 nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                       | inframe_deletion (SO:0001822)      |
| nonframeshift block substitution | 8          | a block substitution of one or more nucleotides that do not cause frameshift changes in protein coding sequence                                                                                                                                                                                     | inframe_variant (SO:0001650)       |
| nonsynonymous SNV                | 9          | a single nucleotide change that cause an amino acid change                                                                                                                                                                                                                                          | missense_variant (SO:0001583)      |
| synonymous SNV                   | 10         | a single nucleotide change that does not cause an amino acid change                                                                                                                                                                                                                                 | synonymous_variant (SO:0001819)    |
| unknown                          | 11         | unknown function (due to various errors in the gene structure definition in the database file)                                                                                                                                                                                                      | sequence_variant (SO:0001060)      |

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
