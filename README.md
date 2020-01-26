# AVA,Dx workflow

---
* **Input**: VCF file, class labels, (cross-validation) data split schemes, *external gene set to use as features*.

* **Output**: selected genes, model performance.

---
## Step One: VCF file variant QC

Extract individuals of interest (diseased and healthy individuals, ideally from the same cohort, meaning sequenced and called in the same batch).
```
bcftools view -S sampleID.txt source.vcf.gz -Oz -o source_s-selected.vcf.gz
```

Remove variant sites which did not pass the VQSR standard.
```
bcftools filter -i 'FILTER="PASS"' source_s-selected.vcf.gz -Oz -o source_s-selected_v-PASS.vcf.gz
```

Split SNV and InDel calls to separated files because they use different QC thresholds. Current AVA,Dx workflow works mainly with SNPs.
```
bcftools view --types snps source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_snps.vcf.gz

bcftools view --types indels source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_indels.vcf.gz
```

Remove variant sites by site-wise quality. Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150. The thresholds are arbitrarily and empirically determined.
```
bcftools view -i 'QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6' source_s-selected_v-PASS_snps.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz
```

Check individual call quality. Good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4. The thresholds are arbitrarily and empirically determined. Bad individual GTs are converted into missing "./.". Remove variant sites with a low call rate. Low call rate is arbitrarily determined as a call rate < 80%, i.e. missing rate >= 20%.
```
python filterVCF_by_ABAD.py \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

---
## Step Two: Individual checks

Identify outliers in terms of **quality** and **ethnicity**.

Check quality outliers by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing.
```
bcftools stats stats -v -s - source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz > source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt

Rscript source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt output.pdf
```


Check ethnicity with AIM (ancestry informative markers).

* Method 1: AIPS from [Byun *et al*](https://morgan1.dartmouth.edu/~f000q4v/html/aips.html).
```
# Convert VCF file into plink format.

# Merge user file and the reference file.
plink --bfile euro952samples --bmerge input.bed input.bim input.fam --recodeA --out outputA

# Run AIPS-PCA.R and AIPS-AI.R.
```

* Method 2: PCA with SNPRelate package in R.
```
Rscript ethnicity_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

* Method 3: Calculate probabilities of individuals being a known ethnicity by forensic marker frequency production.
```
# Extract only the 55 markers from KiddLab.
bcftools view -R 55markers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_Kidd-markers.vcf.gz

# Calculate the probability using a production method.
Rscript forensic_method.R source_Kidd-markers.vcf.gz
```

**Note that:**

  * Method 1 uses AIMs to infer ethnicity using reference labels (952 ancestry known samples).
  * Method 2 takes all SNPs and do PCA on LD-pruned SNPs to infer ethnicity using reference labels (user defined). 
  * Method 3 uses AIMs to infer ethnicities (known ethnicities).
  
  * Technically, results should be very consistant across all method. But results need manual interpretation by human at current stage.


Check relatedness within datasets.
```
Rscript relatedness_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```


Lastly, remove individual outliers from dataset.
```
bcftools -S ^outliers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.vcf.gz
```
---
## Step Three: Variant annotation and SNAP score calculation
```
convert2annovar.pl -format vcf4old 
```

---
## Step Four: Gene score calculation


---
## Step Five: Feature selection and training (cross-validation)




