# AVA,Dx Pipeline

AVA,Dx (Analysis of Variation for Association with Disease) is a computational method for defining the functional role of DNA variation in complex diseases. AVA,Dx uses exonic variants from whole exome or genome sequencing data to extract a *disease* signal and predict *disease* status.

A precursor of the AVA,Dx approach can be found [here](https://bitbucket.org/ywang202/avadx-meta)

## Inputs
* A single Variant Calling File [VCF] comprising all individuals (samples) in the sequencing panel/cohort. Default assembly is GRCh37/hg19
* Sample Information file [comma separated, CSV file], which contains a list of sample ids and corresponding class (e.g. disease vs healthy) labels. *optional* group/family or fold labels for cross validation can also be provided

## Outputs
* table of *gene scores* per sample for subset of most relevant genes
* AVA,Dx model for **D**iseases **x** and performance evaluation
* pathway enrichment annalysis for subset of most relevant genes

## Available Models
Currently pre-computed models are available for:
- [**Crohns Disease**](https://bitbucket.org/bromberglab/avad-cd) (ava-cd)

<br/>

# Installation
## Requirements
The AVA,Dx pipeline is implemented and python and uses pre-built images runnable via docker or singularity. Thus, the only dependencies to run the pipeline are:

* [Docker](https://www.docker.com/get-started) or [Singularity](https://sylabs.io/singularity)
* Python ≥ 3.8
* 10GB free storage

Confirm you have the required versions of Docker/Singularity and Python installed.

For macOS environments we recommend to use [Homebrew](https://brew.sh/) and [pyenv-virtualenv](https://medium.com/python-every-day/python-development-on-macos-with-pyenv-virtualenv-ec583b92934c) to install the required Python dependencies. Docker can be installed from [here](https://docs.docker.com/docker-for-mac/install/).


## Install python package
```bash
pip install https://bitbucket.org/bromberglab/avadx/get/master.zip
```
## Init pipeline
Init AVA,Dx pipeline, download required images and databases and preprocess datasets
```bash
avadx --init
```
A folder `data` and a config file (`avadx.ini`) will be created in the current working directory with the following content:

```
data
├── annovar
│   └── humandb
│       ├── annovar_downdb_gnomad_exome.log
│       ├── annovar_downdb_gnomad_genome.log
│       ├── annovar_downdb_refGene.log
│       ├── <assembly_version>_gnomad_exome.txt.gz
│       ├── <assembly_version>_gnomad_exome.txt.idx
│       ├── <assembly_version>_gnomad_genome.txt.gz
│       ├── <assembly_version>_gnomad_genome.txt.idx
│       ├── <assembly_version>_refGene.txt
│       ├── <assembly_version>_refGeneMrna.fa
│       └── <assembly_version>_refGeneVersion.txt
├── avadx
│   ├── CPDB_pathways_genesymbol.tab
│   ├── Transcript-ProtLength.csv
│   ├── Transcript-ProtLength_cleaned.csv
│   ├── <assembly_version>_gnomad_exome_allAFabove0.txt.gz
│   ├── <assembly_version>_gnomad_exome_allAFabove0.txt.gz.tbi
│   ├── <assembly_version>_gnomad_genome_allAFabove0.txt.gz
│   ├── <assembly_version>_gnomad_genome_allAFabove0.txt.gz.tbi
│   ├── prot_seqs.fa
│   ├── refseq_mapping.csv
│   ├── varidb.db
│   ├── varidb.log
│   └── varidb.md5
├── ethseq
│   └── models
│       └── Exonic.All.Model.gds
├── refseq
│   ├── <assembly_version>_feature_table.txt
│   └── <assembly_version>_protein.faa
└── opt
```
A detailed list of the data resources can also be found [here](https://bitbucket.org/bromberglab/avadx/src/master/data/README.md).

Any user specific data can be added to the `opt` folder and is made available to customized scoring functions and/or models at `/mnt/data/opt`.

**Note:**
- Depending on your download speed and compute resources initialization will run for **>30 minutes**.
- Initialization requires at least **10GB** free storage.

<br />

# Usage
The following steps are detailed below:

1. Edit configuration file
2. Create sample info file
3. Run the pipeline

<br />

## 1 - Edit configuration file
Open the file *avadx.ini* in the current directory (automatically created when the pipeline was initialized) and update the paths to the two required input files (**vcf file** and **samples info file**) and other parameters as required.

**Note:** You can change the configuration file at runtime by supplying the path to the config file as argument (see `--help`).

<br />

## 2 - Edit sample info file
Create the **samples info file** specified in the *avadx.ini* configuration file as `samples = <PATH_TO_FILE>`. For each *sampleid* of interest supply a *class* label and optionally a *group* label or manual *fold* assignment for cross-validation. If neiter *group* nor *fold* are specified, folds are assignet automatically based on the cross-validationn options defined in the *avadx.ini* configuration file.

**Format:** CSV (comma separated values) format with header idetifying column names (see example below).

**Example:** 
| sampleid |class | *[group]* | *[fold]* |
| -------- | ---- | --------- | -------- |
| sample_1 |  0   |    *1*    |    *1*   |
| sample_2 |  0   |    *1*    |    *2*   |
| sample_3 |  0   |    *2*    |    *4*   |
| sample_4 |  1   |    *2*    |    *1*   |
| sample_5 |  1   |    *3*    |    *2*   |
| sample_6 |  1   |    *3*    |    *3*   |

<br />

## 3 - Run the pipeline
*Different ways to run the AVA,Dx pipeline*

-- Run the pipeline using a config file named *avadx.ini* in the current working directory:

```bash
avadx
```

-- Run the pipeline using a specifig config file:
```bash
avadx <PATH_TO_CONFIG_FILE>
```

-- Write output to a specific directory:
```bash
avadx <PATH_TO_CONFIG_FILE> --wd <PATH_TO_OUTPUT_DIRECTORY>
```

-- Save run in a different subfolder (default is configfile name without extension):
```bash
avadx <PATH_TO_CONFIG_FILE> --wd <PATH_TO_OUTPUT_DIRECTORY> --uid <MY_PIPELINE_RUN_ID>
```

For an extensive list of arguments and parameters use `--help`.

<br />

# Update

Update pyton package:
```bash
pip install https://bitbucket.org/bromberglab/avadx/get/master.zip
```

Update VM images:
```bash
avadx --update vm
```

Update Databases:
```bash
avadx --update data
```

<br />

# Modify
## Change scoring functions

AVA,Dx provides a plugin based, flexible way to alter variant scoring and gene scoring aggregation functions.
In the configuration (.ini) file simply specify the path to a python file containing a single scoring function for variants (`variantscore.fn = <PATH>`) and/or genes (`genescore.fn = <PATH>`).
Templates for both scoring functions can be found in this repository under [python/scoring_functions](https://bitbucket.org/bromberglab/avadx/src/master/python/scoring_functions/).
By default, AVA,Dx only scores *SNP* variants. To also include *INDELS* for scoring update the config file accordingly and specify an approproate scoring in the varian scoring function.
Refer to the Annovar variant types table below for more details for scorable variant types.
Additional data resources required can be stored in the `opt` folder and are made available at `/mnt/data/opt` (see *Init pipeline*).

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

<br />

## Change feature selection and model implementations
Both feature selection and classifier/model implementations can be altered using the same plugin system as for the scoring functions.
In the configuration (.ini) file simply specify the path to a python file containing alternative feature selection (`fselection.class = <PATH>`) or classifier/model (`model.class = <PATH>`) implementations.
Templates for feature selection implementations can be found in this repository under [python/feature_selections](https://bitbucket.org/bromberglab/avadx/src/master/python/feature_selections/).
Templates for model implementations can be found in this repository under [python/models](https://bitbucket.org/bromberglab/avadx/src/master/python/models/).
Additional data resources required can be stored in the `opt` folder and are made available at `/mnt/data/opt` (see *Init pipeline*).

<br />

# Notes

* The pipeline defaults to hg19. For hg18, we recommend lifting over to hg19 first. For hg38, we recommend changeing the reference databases to hg38.
* This pipeline is currently for regular VCF file input (modifications are needed in order to use gVCF files).
* Manual interpretation of quality outliers, ethnicity, is a manual step strongly ancouraged.
* When the input VCF contains variants with no SNAP score records available, the *varidb* database needs to be updated.

<br />

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
