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

setwd(out_fp)

if(file.exists(paste0(gsub(".txt","",gs_fp), ".NAto0.nzv85-15.txt"))){
  df_input <- read.csv(paste0(gsub(".txt","",gs_fp), ".NAto0.nzv85-15.txt"), stringsAsFactors=F)
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
if(file.exists("10F-CV-KS-selectedGenes.xlsx")){
  ks_FS_result <- read.xlsx("10F-CV-KS-selectedGenes.xlsx", 1, colIndex=c(2:13))
  ks_FS_result$Gene <- as.character(ks_FS_result$Gene)
  ks_FS_result$Transcript <- as.character(ks_FS_result$Transcript)
}else{
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

  ks_fs_res <- list()
  for(i in 1:k){
    print(i)
    ks_fs <- ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],])
    ks_fs_res[[i]] <- ks_fs
    rm(ks_fs)
  }
  
  #  Summary:
  # unlist(lapply(ks_fs_res, function(x) sum(x < 0.05)))
  
  ks_FS_result <- data.frame("Gene"=df$Gene[-nzv],
                             "Transcript"=df$Transcript[-nzv],
                             "Pvalue_folds1out"=ks_fs_res[[1]],
                             "Pvalue_folds2out"=ks_fs_res[[2]],
                             "Pvalue_folds3out"=ks_fs_res[[3]],
                             "Pvalue_folds4out"=ks_fs_res[[4]],
                             "Pvalue_folds5out"=ks_fs_res[[5]],
                             "Pvalue_folds6out"=ks_fs_res[[6]],
                             "Pvalue_folds7out"=ks_fs_res[[7]],
                             "Pvalue_folds8out"=ks_fs_res[[8]],
                             "Pvalue_folds9out"=ks_fs_res[[9]],
                             "Pvalue_folds10out"=ks_fs_res[[10]],
                             stringsAsFactors=F)
  
  write.xlsx(ks_FS_result, "10F-CV-KS-selectedGenes.xlsx")
}


cv_performance_results <- list()
for(gn in seq(5, 300, 5)){
  print(paste0("Now running gene number: ", gn, " ..."))
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
                                            match(ks_FS_result$Gene[order(ks_FS_result[,(2+i)])][1:gn], df$Gene)),
                                     "status")],
                       probability=T,
                       case.weights=weights
    )
    
    pred_rf <- predict(model_rf, 
                       test_dat[, c(paste0("V", 
                                           match(ks_FS_result$Gene[order(ks_FS_result[,(2+i)])][1:gn], df$Gene)))])
    
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

write.xlsx(cv_performance_results.df, "10F-CV-rf-performance.xlsx", sheetName="KS", append=T)

print("Done.")
