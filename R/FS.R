#library(CORElearn)
library(Boruta)
library(data.table)
library(dplyr)
library(ggplot2)
library(caret)

# args <- commandArgs(trailingOnly=TRUE)
args <- c("/Users/WangYanran/Desktop/tmp/GeneScoreTable_normed.txt", "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Phenotype/cv scheme.txt", "/Users/WangYanran/Desktop/tmp/FS_result")

gs_fp <- args[1]
cvsch_fp <- args[2]
out_fp <- args[3]

# Read-in files:
df <- fread(gs_fp, data.table=F)
cvsch <- fread(cvsch_fp, data.table=F)
cvsch$sample_id <- paste0("sample.", cvsch$sample_ID)

# Extract individuals from cv_scheme file only:
df <- df[, c("Gene", "Transcript", cvsch$sample_id)]

# Count NA gene score cells between disease and healthy:
count_na <- df %>%
  select(-one_of(c("Gene", "Transcript"))) %>%
  summarise_all(funs(sum(is.na(.))))

count_na_df <- data.frame("Sample_ID" = colnames(count_na),
                          "Count_NA" = unlist(count_na[1,]))
count_na_df$status <- cvsch$status[match(count_na_df$Sample_ID, cvsch$sample_id)]

setwd(out_fp)
# Compare:
# test_ks <- ks.test(count_na_df$Count_NA[count_na_df$status==unique(count_na_df$status)[1]],
#                    count_na_df$Count_NA[count_na_df$status==unique(count_na_df$status)[2]])
# 
# ggplot(count_na_df, aes(x=Count_NA, color=status)) + 
#   geom_histogram(fill="white", alpha=0.5, position="identity", binwidth=50) +
#   facet_grid(status ~ .) +
#   theme_classic(base_size=20)



# FS by attrEval:
df_fs <- as.data.frame(t(df[, 3:ncol(df)]))
df_fs[is.na(df_fs)] <- 0
# Remove nzv:
nzv <- nearZeroVar(df_fs)
df.X <- df_fs[, -nzv]
df.y <- cvsch$status[match(rownames(df_fs), cvsch$sample_id)]

df.nzv <- df_fs[, nzv]

df_input <- cbind(df.X, df.y)
colnames(df_input)[ncol(df_input)] <- "status"

df.nzv_input <- cbind(df.nzv, df.y)
colnames(df.nzv_input)[ncol(df.nzv_input)] <- "status"

# Clean-up variables:
# rm(count_na, count_na_df, df)

# fs1_ReliefFexpRank <- attrEval(status ~ ., data=df_fs[rownames(df_fs) %in% cvsch$sample_id[cvsch$fold_num==1],], estimator="ReliefFexpRank")
# setwd(out_fp)
# saveRDS(fs1_ReliefFexpRank, "fs1_ReliefFexpRank.rds")

setwd(out_fp)

# boruta_fs12 <- Boruta(status ~ ., 
#                       data=df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 2)],], 
#                       doTrace=2,
#                       maxRuns=150)
# saveRDS(boruta_fs12, "boruta_fs12.rds")
# 
# boruta_fs13 <- Boruta(status ~ ., 
#                       data=df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 3)],], 
#                       doTrace=2,
#                       maxRuns=150)
# saveRDS(boruta_fs12, "boruta_fs13.rds")
# 
# boruta_fs23 <- Boruta(status ~ ., 
#                       data=df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2, 3)],], 
#                       doTrace=2,
#                       maxRuns=150)
# saveRDS(boruta_fs12, "boruta_fs23.rds")

# FS with anovaScores:
# anova_FS <- function(df){
#   p_values <- c()
#   for(i in 1:(ncol(df)-1)){
#     pval <- anovaScores(df[,i], df$status)
#     p_values <- c(p_values, pval)
#   }
#   return(p_values)
# }
# 
# anova_fs12 <- anova_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 2)],])
# anova_fs13 <- anova_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 3)],])
# anova_fs23 <- anova_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2, 3)],])
# #  Summary:
# # sum(anova_fs12<0.05)
# # sum(anova_fs13<0.05)
# # sum(anova_fs23<0.05)
# anova_FS_result <- data.frame("Gene"=df$Gene[-nzv],
#                               "Transcript"=df$Transcript[-nzv],
#                               "Pvalue_folds12"=anova_fs12,
#                               "Pvalue_folds13"=anova_fs13,
#                               "Pvalue_folds23"=anova_fs23)
# Pcut <- 0.3
# anova_FS_result_sig <- anova_FS_result %>%
#   filter(Pvalue_folds12<Pcut & Pvalue_folds13<Pcut & Pvalue_folds23<Pcut)
# write.table(anova_FS_result_sig$Gene, "/Users/WangYanran/Desktop/tmp/FS_result/anova_list.txt", quote=F, col.names=F, row.names=F)
# 
# tmp <- anova_FS_result %>% arrange(Pvalue_folds12)



# FS with KStest:
ks_FS <- function(df){
  p_values <- c()
  for(i in 1:(ncol(df)-1)){
    ks.res <- ks.test(df[df$status=="proband", i], 
                      df[df$status=="parent", i])
    pval <- ks.res$p.value
    p_values <- c(p_values, pval)
  }
  return(p_values)
}


ks_fs12 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 2)],])
ks_fs13 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 3)],])
ks_fs23 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2, 3)],])
ks_fs1 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1)],])
ks_fs2 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2)],])
ks_fs3 <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(3)],])

ks_fs123 <- ks_FS(df_input)

#  Summary:
sum(ks_fs12<0.05)
sum(ks_fs13<0.05)
sum(ks_fs23<0.05)
sum(ks_fs1<0.05)
sum(ks_fs2<0.05)
sum(ks_fs3<0.05)
ks_FS_result <- data.frame("Gene"=df$Gene[-nzv],
                           "Transcript"=df$Transcript[-nzv],
                           "Pvalue_folds12"=ks_fs12,
                           "Pvalue_folds13"=ks_fs13,
                           "Pvalue_folds23"=ks_fs23,
                           "Pvalue_folds1"=ks_fs1,
                           "Pvalue_folds2"=ks_fs2,
                           "Pvalue_folds3"=ks_fs3,
                           "Pvalue_folds123"=ks_fs123)
ks_FS_result <- ks_FS_result %>% mutate_if(is.factor, as.character)

Pcut <- 0.2
ks_FS_result_sig <- ks_FS_result %>%
  filter(Pvalue_folds1<Pcut & Pvalue_folds2<Pcut & Pvalue_folds3<Pcut)
ks_FS_result_sig$Gene
#write.table(ks_FS_result_sig$Gene, "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/ks_list.txt", quote=F, col.names=F, row.names=F)

ks_genes_folds12 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds12)][1:350]
ks_genes_folds13 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds13)][1:350]
ks_genes_folds23 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds23)][1:350]
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result")
write.table(ks_genes_folds12, "ks_genes_folds12.txt", quote=F, row.names=F, col.names=F)
write.table(ks_genes_folds13, "ks_genes_folds13.txt", quote=F, row.names=F, col.names=F)
write.table(ks_genes_folds23, "ks_genes_folds23.txt", quote=F, row.names=F, col.names=F)
write.table(ks_FS_result$Gene, "ks_background.txt", quote=F, row.names=F, col.names=F)
Reduce(intersect, list(ks_genes_folds12, ks_genes_folds13, ks_genes_folds23)) # ASCC3

ks_genes_folds1 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds1)][1:100]
ks_genes_folds2 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds2)][1:100]
ks_genes_folds3 <- ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds3)][1:100]
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result")
write.table(ks_genes_folds1, "ks_genes_folds1.txt", quote=F, row.names=F, col.names=F)
write.table(ks_genes_folds2, "ks_genes_folds2.txt", quote=F, row.names=F, col.names=F)
write.table(ks_genes_folds3, "ks_genes_folds3.txt", quote=F, row.names=F, col.names=F)
write.table(ks_FS_result$Gene, "ks_background.txt", quote=F, row.names=F, col.names=F)
Reduce(intersect, list(ks_genes_folds1, ks_genes_folds2, ks_genes_folds3)) # OR5R1, ADAMTS16


# Check ASCC3 gene:
match("ASCC3", ks_FS_result$Gene)
df_ASCC3 <- df_input[, c(519, 12261)]
ggplot(df_ASCC3, aes(x=status, y=V529)) + geom_boxplot()

ggplot(df_ASCC3, aes(x=V529, color=status)) +
  geom_histogram(fill="white", alpha=1, position="identity", binwidth=0.01) +
  facet_grid(status ~ .) + xlab("ASCC3 gene score") +
  theme_classic(base_size=20)

# Check TSPO gene: (selected from entire)
match("TSPO", ks_FS_result$Gene)
df_TSPO <- df_input[, c(6779, 12261)]
ggplot(df_TSPO, aes(x=status, y=V7005)) + geom_boxplot()
ggplot(df_TSPO, aes(x=V7005, color=status)) +
  geom_histogram(fill="white", alpha=1, position="identity", binwidth=0.01) +
  facet_grid(status ~ .) + xlab("TSPO gene score") +
  theme_classic(base_size=20)

# Check HS3ST6 gene:
match("HS3ST6", ks_FS_result$Gene)
df_HS3ST6 <- df_input[, c(2880, 12261)]
ggplot(df_HS3ST6, aes(x=status, y=V2975)) + geom_boxplot()

# Check WASHC2C gene:
match("WASHC2C", ks_FS_result$Gene)
df_WASHC2C <- df_input[, c(9052, 12261)]
ggplot(df_WASHC2C, aes(x=status, y=V9408)) + geom_boxplot()

# Check OR5R1 gene:
match("OR5R1", ks_FS_result$Gene)
df_OR5R1 <- df_input[, c(4543, 12261)]
ggplot(df_OR5R1, aes(x=status, y=V4686)) + geom_boxplot()

# Check ADAMTS16 gene:
match("ADAMTS16", ks_FS_result$Gene)
df_ADAMTS16 <- df_input[, c(146, 12261)]
ggplot(df_ADAMTS16, aes(x=status, y=V149)) + geom_boxplot()

############################
# train and test:
library(e1071)
dat12 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 2)],]
dat13 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1, 3)],]
dat23 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2, 3)],]

dat1 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1)],]
dat2 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2)],]
dat3 <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(3)],]


gn <- 8
gs <- c("ASCC3", "HS3ST6", "WASHC2C")
gs <- c("OR5R1", "ADAMTS16")

gs <- c("SLITRK1", "HDC", "DRD3", "DRD2", "SCL6A3", "MET", "COL27A1", "POU1F1")
gs %in% df$Gene
gs %in% ks_FS_result$Gene
gs <- c("HDC", "DRD3", "DRD2", "MET", "COL27A1")

# train:
model12 <- svm(status ~ .,
               dat12[, c(paste0("V", 
                                #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds12)][1:gn], df$Gene)),
                                match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat12$status))

model13 <- svm(status ~ .,
               dat13[, c(paste0("V", 
                                #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds13)][1:gn], df$Gene)),
                                match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat13$status))

model23 <- svm(status ~ .,
               dat23[, c(paste0("V", 
                                #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds23)][1:gn], df$Gene)),
                                match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat23$status))

model1 <- svm(status ~ .,
               dat1[, c(paste0("V", 
                                #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds1)][1:gn], df$Gene)),
                                match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat1$status))

model2 <- svm(status ~ .,
               dat2[, c(paste0("V", 
                               #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds2)][1:gn], df$Gene)),
                               match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat2$status))

model3 <- svm(status ~ .,
               dat3[, c(paste0("V", 
                                #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds3)][1:gn], df$Gene)),
                                match(gs, df$Gene)),
                         "status")],
               class.weights = 100 / table(dat3$status))

# test:
pred3 <- predict(model12, 
                 dat3[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds12)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))])
pred2 <- predict(model13, 
                 dat2[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds13)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))])
pred1 <- predict(model23, 
                 dat1[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds23)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))])

pred23 <- predict(model1, 
                 dat23[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds1)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))],
                 decision.values=T)

pred13 <- predict(model2, 
                 dat13[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds2)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))],
                 decision.values=T)

pred12 <- predict(model3, 
                 dat12[, c(paste0("V", 
                                 #match(ks_FS_result$Gene[order(ks_FS_result$Pvalue_folds3)][1:gn], df$Gene)))])
                                 match(gs, df$Gene)))],
                 decision.values=T)

pred23.num <- attr(pred23, "decision.values")
pred12.num <- attr(pred12, "decision.values")
pred13.num <- attr(pred13, "decision.values")

# eval1 <- confusionMatrix(pred1, dat1$status, positive="proband")
# eval2 <- confusionMatrix(pred2, dat2$status, positive="proband")
# eval3 <- confusionMatrix(pred3, dat3$status, positive="proband")
eval12 <- confusionMatrix(pred12, dat12$status, positive="proband")
eval23 <- confusionMatrix(pred23, dat23$status, positive="proband")
eval13 <- confusionMatrix(pred13, dat13$status, positive="proband")

# eval1$byClass["Balanced Accuracy"]
# eval2$byClass["Balanced Accuracy"]
# eval3$byClass["Balanced Accuracy"]
eval13$byClass["Balanced Accuracy"]
eval23$byClass["Balanced Accuracy"]
eval13$byClass["Balanced Accuracy"]

library(PRROC)
roc.curve(pred13.num[which(dat13$status=="proband")],
          pred13.num[which(dat13$status=="parent")])
roc.curve(pred12.num[which(dat12$status=="proband")],
          pred12.num[which(dat12$status=="parent")])
roc.curve(pred23.num[which(dat23$status=="proband")],
          pred23.num[which(dat23$status=="parent")])




boruta1 <- readRDS("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs1.rds")
boruta2 <- readRDS("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs2.rds")
boruta3 <- readRDS("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs3.rds")

boruta1
boruta2
boruta3

write.table(df$Gene[c(1782     ,2702     ,4843     ,5932    ,10137    ,10283    ,11552 )],
            "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs1.txt",
            quote=F, col.names=F, row.names=F)
write.table(df$Gene[c(1945     ,3216     ,6496     ,7663     ,8177     ,9873    ,11377)],
            "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs2.txt",
            quote=F, col.names=F, row.names=F)
write.table(df$Gene[c(116     ,2435     ,2928     ,5734     ,6852     ,9317 )],
            "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/boruta_fs3.txt",
            quote=F, col.names=F, row.names=F)
write.table(df$Gene[-nzv],
            "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper_TS/Filtered_Data/gene_scores_sum/FS_result/gene_rm-nzv_background.txt",
            quote=F, col.names=F, row.names=F)





### Do KS-test for the nzv genes:
ks_nzv_fs1 <- ks_FS(df.nzv_input[rownames(df.nzv_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(1)],])
ks_nzv_fs2 <- ks_FS(df.nzv_input[rownames(df.nzv_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(2)],])
ks_nzv_fs3 <- ks_FS(df.nzv_input[rownames(df.nzv_input) %in% cvsch$sample_id[cvsch$fold_num%in%c(3)],])

ks_nzv_FS_result <- data.frame("Gene"=df$Gene[nzv],
                           "Transcript"=df$Transcript[nzv],
                           "Pvalue_folds1"=ks_nzv_fs1,
                           "Pvalue_folds2"=ks_nzv_fs2,
                           "Pvalue_folds3"=ks_nzv_fs3)
ks_nzv_FS_result <- ks_nzv_FS_result %>% mutate_if(is.factor, as.character)



