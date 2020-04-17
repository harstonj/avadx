suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(caret))
suppressMessages(library(PRROC))
suppressMessages(library(ranger))
suppressMessages(library(xlsx))

# This script takes in GeneScoreTable_normed.txt, cv-scheme1.txt, and an output folder path
# It outputs 10-fold cross-validation performance with randomforest the enriched gene lists with KS test as FS
# Default gene number tested: seq(5, 300, 5)

# Note that some fields need to be changed in the script:
# ~ line 25: "parent" "proband" only for TS
# ~ line 46: "parent" "proband" only for TS

args <- commandArgs(trailingOnly=TRUE)

gs_fp <- args[1]    # path to GeneScoreTable_normed.txt
cvsch_fp <- args[2] # path to cv-scheme1.txt
out_fp <- args[3]   # path to output folder
fp_map_gene_col <- args[4]
method <- args[5]  # DKM, DKMc, Relieff

setwd(out_fp)

if(file.exists(paste0(gsub(".txt","",gs_fp), ".NAto0.nzv85-15.txt"))){
  df_input <- read.csv(paste0(gsub(".txt","",gs_fp), ".NAto0.nzv85-15.txt"), stringsAsFactors=F)
  colname_map <- read.table(fp_map_gene_col, stringsAsFactors=F, header=T)
  
  cvsch <- read.table(cvsch_fp, stringsAsFactors=F)
  cvsch$sample_id <- paste0("sample.", cvsch$SampleID)
  
  cvsch$status <- plyr::mapvalues(cvsch$Phenotype, c(1, 2), c("parent", "proband"))
  
  df <- fread(gs_fp, data.table=F)
  df <- df[, c("Gene", "Transcript", cvsch$sample_id)]
}else{
  # Read-in files:
  df <- fread(gs_fp, data.table=F)
  cvsch <- read.table(cvsch_fp, stringsAsFactors=F)
  cvsch$sample_id <- paste0("sample.", cvsch$SampleID)
  
  cvsch$status <- plyr::mapvalues(cvsch$Phenotype, c(1, 2), c("parent", "proband"))
  
  # Extract individuals from cv_scheme file only:
  df <- df[, c("Gene", "Transcript", cvsch$sample_id)]
  
  colname_map <- df[, c("Gene", "Transcript")]
  colname_map$col_name <- paste("V", c(1:nrow(df)), sep="")
  write.table(colname_map,
              "map_colname_gene_transcript.txt", quote=F, col.names=T, sep="\t", row.names=F)
  
  df_fs <- as.data.frame(t(df[, 3:ncol(df)]))
  df_fs[is.na(df_fs)] <- 0
  # Remove nzv:
  nzv <- nearZeroVar(df_fs, 85/15)
  df.X <- df_fs[, -nzv]
  df.y <- cvsch$status[match(rownames(df_fs), cvsch$sample_id)]
  
  df_input <- cbind(df.X, df.y)
  colnames(df_input)[ncol(df_input)] <- "status"
  
  write.table(df_input, paste0(gsub(".txt","",gs_fp), ".NAto0.nzv85-15.txt"), sep=",", quote=F, col.names=T, row.names=T)
  
}



k = 10
DKM <- list()
for(i in 1:10){
  DKM[[i]] <- readRDS(paste0("DKM/fs_DKM_fold",i,".rds"))
}

DKMc <- list()
for(i in 1:10){
  DKMc[[i]] <- readRDS(paste0("DKMc/fs_DKMcost_fold",i,".rds"))
}

Relieff <- list()
for(i in 1:10){
  Relieff[[i]] <- readRDS(paste0("ReliefFexpRank/fs_ReliefFexpRank_fold", i, ".rds"))
}


#  Summary:
left_gene.cols <- names(DKM[[1]])
if(method=="DKM"){
  FS_result <- data.frame("Gene"=colname_map$Gene[match(left_gene.cols, colname_map$col_name)],
                          "Transcript"=colname_map$Transcript[match(left_gene.cols, colname_map$col_name)],
                          "DKM_folds1out"=DKM[[1]],
                          "DKM_folds2out"=DKM[[2]],
                          "DKM_folds3out"=DKM[[3]],
                          "DKM_folds4out"=DKM[[4]],
                          "DKM_folds5out"=DKM[[5]],
                          "DKM_folds6out"=DKM[[6]],
                          "DKM_folds7out"=DKM[[7]],
                          "DKM_folds8out"=DKM[[8]],
                          "DKM_folds9out"=DKM[[9]],
                          "DKM_folds10out"=DKM[[10]],
                          stringsAsFactors=F)
}else if(method=="DKMc" | method=="DKMcost"){
  FS_result <- data.frame("Gene"=colname_map$Gene[match(left_gene.cols, colname_map$col_name)],
                          "Transcript"=colname_map$Transcript[match(left_gene.cols, colname_map$col_name)],
                          "DKMc_folds1out"=DKMc[[1]],
                          "DKMc_folds2out"=DKMc[[2]],
                          "DKMc_folds3out"=DKMc[[3]],
                          "DKMc_folds4out"=DKMc[[4]],
                          "DKMc_folds5out"=DKMc[[5]],
                          "DKMc_folds6out"=DKMc[[6]],
                          "DKMc_folds7out"=DKMc[[7]],
                          "DKMc_folds8out"=DKMc[[8]],
                          "DKMc_folds9out"=DKMc[[9]],
                          "DKMc_folds10out"=DKMc[[10]],
                          stringsAsFactors=F)
  
}else if(method=="ReliefF" | method=="Relieff" | method=="relieff"){
  FS_result <- data.frame("Gene"=colname_map$Gene[match(left_gene.cols, colname_map$col_name)],
                          "Transcript"=colname_map$Transcript[match(left_gene.cols, colname_map$col_name)],
                          "Relieff_folds1out"=Relieff[[1]],
                          "Relieff_folds2out"=Relieff[[2]],
                          "Relieff_folds3out"=Relieff[[3]],
                          "Relieff_folds4out"=Relieff[[4]],
                          "Relieff_folds5out"=Relieff[[5]],
                          "Relieff_folds6out"=Relieff[[6]],
                          "Relieff_folds7out"=Relieff[[7]],
                          "Relieff_folds8out"=Relieff[[8]],
                          "Relieff_folds9out"=Relieff[[9]],
                          "Relieff_folds10out"=Relieff[[10]],
                          stringsAsFactors=F)
}

cv_performance_results <- list()
for(gn in seq(5, 300, 5)){
  # gn <- 10
  print(gn)
  pred_res_rf <- list()
  for(i in 1:k){
    # print(i)
    train_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],]
    test_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% i],]
    
    w <- 1/table(train_dat$status)
    w <- w/sum(w)
    
    weights <- rep(0, nrow(train_dat))
    weights[train_dat$status=="proband"] <- w["proband"]
    weights[train_dat$status=="parent"] <- w["parent"]
    
    model_rf <- ranger(status ~ .,
                       train_dat[, c(paste0("V", 
                                            match(FS_result$Gene[order(-FS_result[,(2+i)])][1:gn], df$Gene)),
                                     "status")],
                       probability=T,
                       case.weights=weights
    )
    
    pred_rf <- predict(model_rf, 
                       test_dat[, c(paste0("V", 
                                           match(FS_result$Gene[order(-FS_result[,(2+i)])][1:gn], df$Gene)))])
    
    pred_num <- pred_rf$predictions[,"proband"]
    pred_res_rf[[i]] <- data.frame("status"=test_dat$status,
                                   "prediction"=pred_num,
                                   stringsAsFactors=F)
    rm(train_dat, test_dat, model_rf, pred_rf, pred_num, weights, w)
  }
  
  # RF performance evaluation:
  pred_res_rf.df <- bind_rows(pred_res_rf)
  roc_rf <- roc.curve(pred_res_rf.df$prediction[pred_res_rf.df$status=="proband"],
                      pred_res_rf.df$prediction[pred_res_rf.df$status=="parent"])
  # print(paste0(gn, " genes: ", roc_rf$auc))
  cv_performance_results[[paste0("GeneNumber.", gn)]] <- roc_rf$auc
}


cv_performance_results.df <- as.data.frame(do.call("rbind", cv_performance_results))
colnames(cv_performance_results.df)[1] <- "randomForest_AUC"
cv_performance_results.df$GeneNumber <- seq(5, 300, 5)

if(method=="DKM"){
  write.xlsx(cv_performance_results.df, "10F-CV-rf-performance.xlsx", sheetName="DKM", append=T)
}else if(method=="DKMc" | method=="DKMcost"){
  write.xlsx(cv_performance_results.df, "10F-CV-rf-performance.xlsx", sheetName="DKMc", append=T)
}else if(method=="ReliefF" | method=="Relieff" | method=="relieff"){
  write.xlsx(cv_performance_results.df, "10F-CV-rf-performance.xlsx", sheetName="ReliefF", append=T)
}

print("Done.")
