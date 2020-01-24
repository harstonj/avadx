# AVA,Dx workflow

---
* Input: VCF file, class labels, split schemes, (external gene set to use as features)

* Output: selected genes, model performance

---
## Step One: VCF file clean-up

Extract individuals of interest (diseased and healthy individuals, ideally from the same cohort, meaning sequenced and called in the same batch).
```
bcftools view -S sampleID.txt source.vcf.gz -Oz -o source_s-selected.vcf.gz
```

Remove variant sites which did not pass the VQSR standard.
```
bcftools filter -i 'FILTER="PASS"' source_s-selected.vcf.gz -Oz -o source_s-selected_v-PASS.vcf.gz
```

Split SNV and InDel calls to separated files because they use different QC thresholds.
```
bcftools view --types snps source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_snps.vcf.gz

bcftools view --types indels source_s-selected_v-PASS.vcf.gz -Oz -o source_s-selected_v-PASS_indels.vcf.gz
```

Remove variant sites by site-wise quality. Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150. The thresholds are arbitrarily and empirically determined.
```
bcftools view -i 'QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6' source_s-selected_v-PASS_snps.vcf.gz -Oz -o source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz
```

Check individual call quality. Good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4. The thresholds are arbitrarily and empirically determined. The bad individual GTs should be converted into missing "./.". Remove variant sites with a low call rate. Low call rate is arbitrarily determined as a call rate < 80% or missing rate >= 20%.
```
python filterVCF_by_ABAD.py \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

---
## Step Two: Individual checks

Identify outliers in terms of calling quality and ethnicity.

Check calling quality by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing.
```
bcftools stats stats -v -s - source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz > source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt

Rscript source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt output.pdf
```


Check ethnicity with AIM (ancestry informative markers).
Method 1: AIPS from [Byun *et al*](https://morgan1.dartmouth.edu/~f000q4v/html/aips.html).
```
# Convert VCF file into plink format.

# Merge user file and the reference file.
plink --bfile euro952samples --bmerge input.bed input.bim input.fam --recodeA --out outputA

# Run AIPS-PCA.R and AIPS-AI.R.
```

Method 2: SNPRelate package in R.
```
Rscript ethnicity_SNPRelate.R \
  source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz
```

Method 3: Calculate probabilities individuals being a known ethnicity by forensic marker frequency production.
```
# Extract only the 55 markers from KiddLab.
bcftools view -R 55markers.txt source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz -Oz -o source_Kidd-markers.vcf.gz

# Calculate the probability using a production method.
Rscript forensic_method.R source_Kidd-markers.vcf.gz
```

Note that:
  * Method 1 uses AIMs to infer ethnicity using reference labels (952 ancestry known samples).
  * Method 2 takes all SNPs and do PCA on LD-pruned SNPs to infer ethnicity using reference labels (user defined). 
  * Method 3 uses AIMs to infer ethnicities (known ethnicities).
  
  * Technically, results should be very consistant across all method.

---

**Edit**

When you're done, you can delete the content in this README and update the file with details for others getting started with your repository.

*We*

---

## Edit a file

You’ll start by editing this README file to learn how to edit a file in Bitbucket.

1. Click **Source** on the left side.
2. Click the README.md link from the list of files.
3. Click the **Edit** button.
4. Delete the following text: *Delete this line to make a change to the README from Bitbucket.*
5. After making your change, click **Commit** and then **Commit** again in the dialog. The commit page will open and you’ll see the change you just made.
6. Go back to the **Source** page.

---

## Create a file

Next, you’ll add a new file to this repository.

1. Click the **New file** button at the top of the **Source** page.
2. Give the file a filename of **contributors.txt**.
3. Enter your name in the empty file space.
4. Click **Commit** and then **Commit** again in the dialog.
5. Go back to the **Source** page.

Before you move on, go ahead and explore the repository. You've already seen the **Source** page, but check out the **Commits**, **Branches**, and **Settings** pages.

---

## Clone a repository

Use these steps to clone from SourceTree, our client for using the repository command-line free. Cloning allows you to work on your files locally. If you don't yet have SourceTree, [download and install first](https://www.sourcetreeapp.com/). If you prefer to clone from the command line, see [Clone a repository](https://confluence.atlassian.com/x/4whODQ).

1. You’ll see the clone button under the **Source** heading. Click that button.
2. Now click **Check out in SourceTree**. You may need to create a SourceTree account or log in.
3. When you see the **Clone New** dialog in SourceTree, update the destination path and name if you’d like to and then click **Clone**.
4. Open the directory you just created to see your repository’s files.

Now that you're more familiar with your Bitbucket repository, go ahead and add a new file locally. You can [push your change back to Bitbucket with SourceTree](https://confluence.atlassian.com/x/iqyBMg), or you can [add, commit,](https://confluence.atlassian.com/x/8QhODQ) and [push from the command line](https://confluence.atlassian.com/x/NQ0zDQ).