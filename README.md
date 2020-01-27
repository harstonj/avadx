# AVA,Dx workflow

---
* **Input**: VCF file, class labels, (cross-validation) data split schemes, *external gene set to use as features*.
* **Manual Processes Needed**: determination of all arbitrary thresholds along all steps; outlier identification (e.g. ethnicity check); SNAP score calculation (parallelization on amarel); gene score calculation (parallelization on amarel); feature selection and model selection in model training
* **Output**: selected genes, model performance.

---
## Step One: VCF file variant QC

1.Extract individuals of interest (diseased and healthy individuals, ideally from the same cohort, meaning sequenced and called in the same batch).
```
bcftools view -S sampleID.txt source.vcf.gz -Oz -o source_s-selected.vcf.gz
```

2.Remove variant sites which did not pass the VQSR standard.
```
bcftools filter -i 'FILTER="PASS"' source_s-selected.vcf.gz -Oz -o source_s-selected_v-PASS.vcf.gz
```

3.Split SNV and InDel calls to separated files because they use different QC thresholds. Current AVA,Dx workflow works mainly with SNPs.
```
bcftools view --types snps source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_snps.vcf.gz

bcftools view --types indels source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_indels.vcf.gz
```

4.Remove variant sites by site-wise quality. Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150. These thresholds are arbitrarily and empirically determined.
```
bcftools view -i 'QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6' source_s-selected_v-PASS_snps.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz
```

5.Check individual call quality. Good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4. These thresholds are arbitrarily and empirically determined. Bad individual GTs are converted into missing "./.". Remove variant sites with a low call rate. Low call rate is arbitrarily determined as a call rate < 80%, *i.e.* missing rate >= 20%.
```
python filterVCF_by_ABAD.py \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

---
## Step Two: Individual checks

Identify outliers in terms of **quality** and **ethnicity**.

1.Check quality outliers by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing.
```
bcftools stats -v -s - source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz > source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt

Rscript source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt output.pdf
```
The above result needs manual interpretation and selection of outliers. Alternatively, one can use below code to automatically select outliers using a distance-based statistical method.
```
Rscript
```

2.Check ethnicity with AIM (ancestry informative markers) or LD-pruned SNPs.

* *Method 1*: AIPS from [Byun *et al*](https://morgan1.dartmouth.edu/~f000q4v/html/aips.html).
```
# Convert VCF file into plink format.

# Merge user file and the reference file.
plink --bfile euro952samples --bmerge input.bed input.bim input.fam --recodeA --out outputA

# Run AIPS-PCA.R and AIPS-AI.R.
```

* *Method 2*: PCA with [SNPRelate package](http://corearray.sourceforge.net/tutorials/SNPRelate/) in R.
```
Rscript ethnicity_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

* *Method 3*: [EthSEQ](https://cran.r-project.org/web/packages/EthSEQ/index.html) R package.
```
# If the number of individuals exceeds certain number, "memory exhausted" error may occur. Manually divide input VCF into chunks of individuals and run EthSEQ separately for each chunk:
bcftools query -l source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz > sample_list.txt
csplit sample_list.txt 500  # outputs from xx00 to xx0n

# Clean VCF format for EthSEQ input (do the same thing for every chunk):
bcftools view -S xx00 source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz | bcftools annotate --remove 'ID,INFO,FORMAT' | bcftools view --no-header -Oz -o source_xx00_EthSEQinput.vcf.gz
# If no separation of individuals:
bcftools annotate --remove 'ID,INFO,FORMAT' source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz | bcftools view --no-header -Oz -o source_EthSEQinput.vcf.gz

# Run EthSEQ:
Rscript ethnicity_EthSEQ.R source_EthSEQinput.vcf.gz /path/to/output/folder
```
Results of ethnicity predictions are in `/path/to/output/folder/Report.txt` and the corresponding sample IDs are in `sample_list.txt`.

```
Rscript ethnicity_EthSEQ_summary.R /path/to/output/folder/Report.txt sample_list.txt /path/to/output/folder
```
Above returns two files: `sampleID_closest_EUR.txt` and `sampleID_inside_EUR.txt`. Customized script should be used for special requirements. For example, in tourette yale-1 dataset, we keep only trios when all people in the family are EUR.

* *Method 4*: Calculate probabilities of individuals being a [known ethnicity](https://frog.med.yale.edu/FrogKB/FrogServlet) by forensic marker [frequency production](https://frog.med.yale.edu/FrogKB/formula.jsp).
```
# Extract only the 55 markers from KiddLab.
bcftools view -R 55markers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_Kidd-markers.vcf.gz

# Calculate the probability using a production method.
Rscript forensic_method.R source_Kidd-markers.vcf.gz
```

**Note that:**

  * *Method 1* uses AIMs to infer ethnicity using reference labels (952 ancestry known samples).
  * *Method 2* takes all SNPs and do PCA on LD-pruned SNPs to infer ethnicity using reference labels (user defined).
  * *Method 3* uses pre-calculated reference SS2 model and gives prediction of five 1000Genomes population code.
  * *Method 4* uses AIMs to infer ethnicities (known ethnicities).

  * Technically, results should be very **consistent across all method**. But results need manual interpretation.


3.Check relatedness within datasets. A BID > 0.05 is considered to be related.
```
Rscript relatedness_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

  Additionally, check HW equilibrium within diseased / healthy cohort. Check sex (if sex labels are provided).
```
Rscript
```

4.Lastly, remove individual outliers from dataset.
```
bcftools -S ^outliers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.vcf.gz
```

---
## Step Three: Variant annotation and SNAP score calculation
ANNOVAR is used to do variant annotation:
```
# Convert VCF file into ANNOVAR input format:
convert2annovar.pl -format vcf4old source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.vcf.gz -outfile source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.avinput

# Annotate using hg19 human reference:
annotate_variation.pl -buildver hg19 source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.avinput humandb/
```
In the above commands, `vcf4old` argument is used to ensure that all records in the VCF file is output to one resulting file, which contains every possible variant in this dataset. The output file should include `source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function`, which is used as input for the next steps.


AVA,Dx has pre-calculated SNAP scores in the *db* folder (`Mutations.mutOut`). Additionally, `mRNA_identifiers.txt` stores the transcript identifiers (currently 46,327 mRNA transcripts), and `prot_seqs.txt` stores the protein sequences with "NM_XX_prot_NP_YY" in the header lines for mapping uses.

`Mutations.mutOut` has three columns: mRNA accession, amino acid mutation, pre-calculated SNAP score.

For a new dataset, we check if the pre-calculated SNAP score files already have recorded all the variants. If not, we generate a "missing SNAP" list and make SNAP input files for those variants:
```
Rscript check_missing_SNAP.R source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function /path/to/db/folder /path/to/output/folder
```
The above command first check if the SNAP score file `Mutations.mutOut` contains all the variants in the dataset. If not, the script prints out the number of "missing SNAP" mutations and starts to generate SNAP input for those variants. To do this, first step is to check if the `prot_seqs.txt` file contains all protein sequences of the missing SNAP mutations. If not, the script outputs a file with mRNA accession numbers (`TranscriptAccess_missing_prot_seq.txt`) in the output folder. User needs to manually retrieve the corresponding protein sequences at [NCBI Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez). Upload the `TranscriptAccess_missing_prot_seq.txt` file with "Choose File" and *Retrieve* in the *Nucleotide database*. Then click *Send to* at the resulting page, choose *Coding Sequences* and *FASTA Protein* to get the protein sequence. Click *Create File* to download. 

SNAP needs two input files: **amino acid mutation list** and the **protein fasta file** for each protein. Code below extracts all mutations from `.exonic_variant_function` file and write them into SNAP input files.
```
Rscript generate_SNAP_input.R source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function /path/to/output/folder
```
The above step generates one mutation file `geneA.mut` and one fasta file `geneA.fasta` for *geneA* to the output folder. There will be approximately 20,000 to 45,000 files generated depending on the genes covered by the original data.

SNAP is available on [amarel server](https://oarc.rutgers.edu/amarel/) and can be run using code below (submit.sh SBATCH submission shell script):
```
#!/bin/bash

#SBATCH --partition=bromberg_1,main
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --array=0-999
#SBATCH --requeue

module load singularity/.2.4-PR1106
inArray=(geneA geneB geneC)

#sbatch --partition=${4} --array=${1}-${2} submit_partial_test.sh $3 IL9R_HUMAN_1,IL9R_HUMAN_2
singularity exec /home/yw410/bromberglab_predictprotein_yanran-2017-12-06-fa6f97ee098c.img snapfun -i /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_input_part1/$input.fasta -m /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_input_part1/$input.mutation -o /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_output_part1/$input.out --print-collection --tolerate-sift-failure --tolerate-psic-failure

```

SNAP outputs text files of predicted variant scores. It needs to be converted to tab-separated format using below code:
```
python snap_output2tab.py snapOutput snapOutput.tab
```

After all SNAP output has been converted to tab-separated format, merge them into one mutation score table:
```
cat *.tab > mutScore.txt
```

---
## Step Four: Gene score calculation


---
## Step Five: Feature selection and training (cross-validation)
