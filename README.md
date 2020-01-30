# AVA,Dx workflow

---
## Inputs and ouputs

* **Input**:

 * VCF file
 * class labels
 * individual list (sample IDs of interest)
 * (cross-validation) data split schemes (optional)
 * external gene set to use as features (optional)

* **Output**:

 * a *gene score* table
 * selected genes (with default FS method: Kolmogorov–Smirnov test or DKM method from CORElearn R package)
 * model performance (with default settings: SVM from e1071 R package)
 
* **Check before running all steps**:

 * The current workflow works with hg19.
 * This pipeline is currently for regular VCF file input. If gVCF is provided, some of the scripts need to be re-written.
 * All thresholds along all steps are arbitrarily and empirically determined (details below).
 * Outlier identification (e.g. ethnicity check) may need human interpretation.
 * For a new input VCF, pre-calculated SNAP scores may not contain all mutations and SNAP needs to be re-run for the "missing SNAP" mutations.
 * Four different *gene score* calculation methods are available currently (details below).
 * Feature selection (FS) and model selection in model training, including FS method choosing, model choosing, model tuning, etc. need human interpretation.
---
## Prerequisite
* R and packages (data.table, tydiverse, seqinr, stringr, EthSEQ, SNPRelate, e1071)
* python
* [bcftools](https://samtools.github.io/bcftools/)
* [ANNOVAR](http://annovar.openbioinformatics.org)
* [PLINK](https://www.cog-genomics.org/plink2/)

---
## Step One: VCF file variant QC

* Extract individuals of interest (diseased and healthy individuals, ideally from the same cohort, meaning sequenced and called in the same batch).
```
bcftools view -S sampleID.txt source.vcf.gz -Oz -o source_s-selected.vcf.gz
```

* Remove variant sites which did not pass the VQSR standard.
```
bcftools filter -i 'FILTER="PASS"' source_s-selected.vcf.gz -Oz -o source_s-selected_v-PASS.vcf.gz
```

* Split SNV and InDel calls to separated files because they use different QC thresholds. Current AVA,Dx workflow works mainly with SNPs.
```
bcftools view --types snps source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_snps.vcf.gz
bcftools view --types indels source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_indels.vcf.gz
```

* Remove variant sites by site-wise quality. Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150. These thresholds are arbitrarily and empirically determined.
```
bcftools view -i 'QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6' source_s-selected_v-PASS_snps.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz
```

* Check individual call quality. Good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4. These thresholds are arbitrarily and empirically determined. Bad individual GTs are converted into missing "./.". Remove variant sites with a low call rate. Low call rate is arbitrarily determined as a call rate < 80%, *i.e.* missing rate >= 20%.
```
python filterVCF_by_ABAD.py \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

---
## Step Two: Individual checks

Identify outliers in terms of **quality** and **ethnicity**.

* Check quality outliers by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing.
```
bcftools stats -v -s - source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz > source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt
Rscript source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt output.pdf
```

* The above result needs manual interpretation and selection of outliers. Alternatively, one can use below code to automatically select outliers using a distance-based statistical method.
```
Rscript
```

* Check ethnicity with AIM (ancestry informative markers) or LD-pruned SNPs.

 * *Method 1*: AIPS from [Byun *et al*](https://morgan1.dartmouth.edu/~f000q4v/html/aips.html).
```
# Convert VCF file into plink format.
plink xx
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

 * Results of ethnicity predictions are in `/path/to/output/folder/Report.txt` and the corresponding sample IDs are in `sample_list.txt`.
```
Rscript ethnicity_EthSEQ_summary.R /path/to/output/folder/Report.txt sample_list.txt /path/to/output/folder
```
 * Above returns two files: `sampleID_closest_EUR.txt` and `sampleID_inside_EUR.txt`. Customized script should be used for special requirements.

  * *Method 4*: Calculate probabilities of individuals being a [known ethnicity](https://frog.med.yale.edu/FrogKB/FrogServlet) by forensic marker [frequency production](https://frog.med.yale.edu/FrogKB/formula.jsp).
```
# Extract only the 55 markers from KiddLab.
bcftools view -R 55markers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_Kidd-markers.vcf.gz
# Calculate the probability using a production method.
Rscript forensic_method.R source_Kidd-markers.vcf.gz
```

* **Note that:**
   * *Method 1* uses AIMs to infer ethnicity using reference labels (952 ancestry known samples).
   * *Method 2* takes all SNPs and do PCA on LD-pruned SNPs to infer ethnicity using reference labels (user defined).
   * *Method 3* uses pre-calculated reference SS2 model and gives prediction of five 1000Genomes population code.
   * *Method 4* uses AIMs to infer ethnicities (known ethnicities).
   * Technically, results should be very **consistent across all method**. But human interpretation may be needed for specific cases.


* Then, check relatedness within datasets. A BID > 0.05 is considered to be related.
```
Rscript relatedness_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

* Additionally, check HW equilibrium within diseased / healthy cohort. Also, check sex if sex labels are provided.
```
Rscript
```

* Lastly, remove individual outliers from dataset.
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
In the above commands, `vcf4old` argument is used to ensure that all records in the VCF file is output to one resulting file, which contains every possible variant in this dataset. The output file should include `source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function`, which is used as input for steps below.

AVA,Dx has pre-calculated SNAP scores stored in the *db* folder (`Mutations.mutOut`). Additionally in the same folder, `mRNA_identifiers.txt` stores the transcript identifiers (currently 46,327 mRNA transcripts), and `prot_seqs.txt` stores the protein sequences with "NM_XX_prot_NP_YY" in the header lines for mapping uses.

`Mutations.mutOut` has three columns: mRNA accession, amino acid mutation, pre-calculated SNAP score.

For a new dataset, we first check if the pre-calculated SNAP score file already has recorded all the variants. If not, we generate a "missing SNAP" list and make SNAP input files for those variants:
```
Rscript check_missing_SNAP.R source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function /path/to/db/folder /path/to/output/folder
```
The above command first check if the SNAP score file `Mutations.mutOut` contains all the variants in the dataset. If not, the script prints out the number of "missing SNAP" mutations and generates SNAP input for those variants. To do this, the script first checks if `prot_seqs.txt` file contains all protein sequences and then generates SNAP input files for new variants.

If `prot_seqs.txt` is not complete, the script outputs a file with mRNA accession numbers (`TranscriptAccess_missing_prot_seq.txt`) in the output folder. User needs to retrieve the corresponding protein sequences at [NCBI Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez). Upload the `TranscriptAccess_missing_prot_seq.txt` file with "Choose File" and *Retrieve* in the *Nucleotide database*. Then click *Send to* at the resulting page, choose *Coding Sequences* and *FASTA Protein* to get the protein sequence. Click *Create File* to download. Then, append the protein sequences to the original `prot_seqs.txt` file by:
```
cat prot_seqs.txt sequence.txt > prot_seqs_new.txt
rm prot_seqs.txt
mv prot_seqs_new.txt prot_seqs.txt
```

Then update the `Transcript-ProtLength.csv` file in *db* folder:
```
Rscript update_Transcript-ProtLength.R /path/to/db/folder
# Check if the output Transcript-ProtLength_update.csv is correct
rm Transcript-ProtLength.csv
mv Transcript-ProtLength_update.csv Transcript-ProtLength.csv
```

After the protein sequences are updated by above steps, run the `check_missing_SNAP.R` again to generate mutation files for SNAP input:
```
Rscript check_missing_SNAP.R source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.exonic_variant_function /path/to/db/folder /path/to/output/folder
```
All SNAP input files (*amino acid mutation list* `geneA.mutation` and the *protein fasta file* `geneA.fasta`) will be output to `/path/to/output/folder`, This process might take some time if there are many missing SNAPs.


After all SNAP input files are ready, SNAP is available on [amarel server](https://oarc.rutgers.edu/amarel/) and can be run using code below (submit.sh SBATCH submission shell script):
```
#!/bin/bash
#SBATCH --partition=bromberg_1,main
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --array=0-999
#SBATCH --requeue
module load singularity/.2.4-PR1106
input=(geneA geneB geneC ...)
singularity exec /home/yw410/bromberglab_predictprotein_yanran-2017-12-06-fa6f97ee098c.img snapfun -i /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_input/$input.fasta -m /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_input/$input.mutation -o /home/yw410/singularity_in/SNAPinput-dbGaP/SNAP_output/$input.out --print-collection --tolerate-sift-failure --tolerate-psic-failure
```

SNAP outputs plain text files `geneA.out` of predicted variant scores. It needs to be converted to tab-separated format using below code for SNAP output of all genes:
```
#!/bin/bash
for f in /path/to/snap/output/*.out
do
	python snap-scores-mutOut.py $f
done
# Outputs Mutations.mutOut file
```

`Mutations.mutOut` has three columns: gene ID (transcript accession), amino acid mutation, SNAP score. After all new SNAP output has been converted to tab-separated format, merge them with the original SNAP scores stored in the *db* folder:
```
cat /path/to/db/folder/Mutations.mutOut Mutations.mutOut > /path/to/db/folder/Mutations_new.mutOut
cd /path/to/db/folder/
rm Mutations.mutOut
mv Mutations_new.mutOut Mutations.mutOut
```
Now the *db* folder should have:

 * updated SNAP score file `Mutations.mutOut`
 * updated transcript - protein - protein length file `Transcript-ProtLength.csv`

---
## Step Four: Gene score calculation

Before *gene score* calculation, cleaned VCF file should be converted to individual ANNOVAR annotations by:
```
ulimit -n 2048  # set number of open files limit to max
convert2annovar.pl -format vcf4 source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.vcf.gz -outfile /path/to/output/folder/sample -allsample
```
The `vcf4` and `-allsample` arguments specify that one `sample*.avinput` output file will be generated for every individual in the VCF. Note that `convert2annovar.pl` script opens new files simultaneously for all samples, and when the sample number exceeds the maximum number (2048 set by `ulimit -n 2048`), the script returns error because it cannot open a larger number of files. User needs to split the VCF file into chunks of samples and run `convert2annovar.pl` separately for each chunk.

Next, do annotation with hg19 assembly to all `sample*.avinput` files (preferably on amarel using job arrays):
```
#!/bin/bash
#SBATCH --partition=bromberg_1,main
#SBATCH --time=24:00:00
#SBATCH --mem=12288
#SBATCH --array=0-371
#SBATCH --requeue
inArray=(sample1.avinput sample2.avinput sample3.avinput ...)
input=${inArray[$SLURM_ARRAY_TASK_ID]}
annotate_variation.pl -build hg19 $input /humandb
```

*gene score* is defined as a **sum** or **production** of the all variant scores within the gene coding region.

*EQUATIONS*

All gene score calculation functions are stored at `gene_score_schemes.R` script. There are some arguments required for the gene score calculation: `-scheme` asks user to choose either to sum or multiply variant scores into *gene score*; `-heti` asks for the coefficient used for the heterozygous genotypes and `HIPred` means using the happloinsufficienty predictions from [HIPred](https://www.ncbi.nlm.nih.gov/pubmed/28137713); `-input` asks for the ".exonic_variant_function" file of that individual; `-snap_score` asks for the path to the tab-separated file that stores all SNAP score predictions; `-db` asks for the path to the *db* folder; lastly, `-output` asks user to specify the outputting path.
```
# Run this for every person:
Rscript gene_score_schemes.R -scheme sum -heti 0.25 -input sample1.exonic_variant_function -snap_score Mutations.mutOut -db /path/to/db/folder -output /path/to/output
# Or:
Rscript gene_score_schemes.R -scheme prod -heti HIPred -input sample1.exonic_variant_function -snap_score Mutations.mutOut -db /path/to/db/folder -output /path/to/output
```
Assuming there are N individuals in the dataset, N resulting files will be generated with calculated *gene_scores* for all available genes. Below script will read-in all files and merge them into a file table where a row is an individual and a column is a gene (protein). To avoid *redundancy*, we are currently only considering the "canonical" protein of each gene, *i.e.* the longest protein sequence.
```
Rscript merge_individual_scores.R /path/to/individual/score/folder /path/to/output
```

---
## Step Five: Feature selection (FS) and training (cross-validation)

There are three types of FS types: filter, wrapper, and embedded.

Current AVA,Dx does not use the wrapper method because wrapper method usually takes long time to search feature subset and there are always multiple optimal feature sets, which are hard to interpret from a biological perspective.

Current AVA,Dx uses both filter and embedded methods. Filter method includes,
