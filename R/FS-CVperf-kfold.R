suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(caret))
suppressMessages(library(PRROC))
suppressMessages(library(xlsx))
library(parallel)

# This script takes in GeneScoreTable_normalized.txt, cv-scheme1.txt, and an output folder path
# It outputs 10-fold cross-validation performance with randomforest the enriched gene lists with KS test as FS
# Default gene number tested: seq(5, 300, 5)

# Note that some fields need to be changed in the script:
# ~ line 25: "parent" "proband" only for TS
# ~ line 46: "parent" "proband" only for TS

option_list = list(
  make_option(c("-f", "--input_file"), type="character", default=NULL,
              help="path to input file; e.g. GeneScoreTable_normalized.txt", metavar="character"),
  make_option(c("-m", "--fs_method"), type="character", default=NULL,
              help="feature selection method to use: 'ks' - Kolmogorov–Smirnov test between disease and control samples", metavar="character"),
  make_option(c("-M", "--ml_method"), type="character", default=NULL,
              help="machine learning method to use: 'rf' - random forest from package ranger; 'svm' - SVM from package e1071", metavar="character"),
  make_option(c("-v", "--variance_cutoff"), type="numeric", default=99,
              help="a cutoff percentage value to remove features/genes with low variance; default is 99, i.e. the feature/gene is removed if over 99% individuals have the same gene score"),
  make_option(c("-s", "--cv_scheme"), type="character", default=NULL,
              help="path to cross-validation scheme file; should contain columns 'SampleID', 'Phenotype', 'fold'. Note that Phenotype should be 0 and 1, where 0 means control and 1 means case/sick.", metavar="character"),
  make_option(c("-k", "--k_fold"), type="numeric", default=NULL,
              help="feature selection method to use: 'ks' - Kolmogorov–Smirnov test between disease and control samples", metavar="character"),
  make_option(c("-l", "--protlength_file"), type="character", default=NULL, 
              help="protein length file (in the db folder) where columns are: gene name, NM_ number, corresponding protein length from RefSeq; this file helps retain only the longest protein if there are duplicated gene/protein in the input data frame", metavar="character"),
  make_option(c("-t", "--step_of_top_genes"), type="numeric", default=5,
              help="step genes to try out; e.g. 5 (default) means increasing by 5 genes", metavar="character"),
  make_option(c("-n", "--number_of_top_genes"), type="numeric", default=100,
              help="number of top-ranked genes to try out; e.g. 100 (default) means try until the top 100 genes", metavar="character"),
  make_option(c("-K", "--ks_pval_cutoff"), type="numeric", default=0.2,
              help="only keep genes with p-val < 0.2 for model training", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="path to output folder", metavar="character"),
  make_option(c("-w", "--wd"), type="character", default=NULL,
              help="path to working directory", metavar="character")
)

# args <- commandArgs(trailingOnly=TRUE)
# gs_fp <- args[1]    # path to GeneScoreTable_normalized.txt
# cvsch_fp <- args[2] # path to cv-scheme1.txt
# out_fp <- args[3]   # path to output folder

opt = parse_args(OptionParser(option_list=option_list))

fs_method <- opt$fs_method
ml_method <- opt$ml_method
gs_fp <- opt$input_file
cvsch_fp <- opt$cv_scheme
prot_len_fp <- opt$protlength_file
out_fp <- opt$out
wd_fp <- opt$wd
k <- opt$k_fold
ks_pval_cutoff <- opt$ks_pval_cutoff
numCores <- detectCores()

dir.create(file.path(wd_fp), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_fp), recursive = TRUE, showWarnings = FALSE)

setwd(wd_fp)

# Read-in files & Extract individuals from cv_scheme file only:
cvsch <- fread(cvsch_fp, data.table=F)
setnames(cvsch, c('SampleID', 'fold', 'Phenotype'))
if (!grepl("sample.", cvsch$SampleID[1])) {
  cvsch$sample_id <- paste0("sample.", cvsch$SampleID) # In case the sample IDs starts with numbers
} else {
  cvsch$sample_id <- cvsch$SampleID
}

# get folds from cv-scheme.txt
if (is.null(k)) {
  k = length(unique(cvsch$fold))
}

if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx"))){
  file.remove(paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx"))
}
if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"))){
  file.remove(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"))
}
if(file.exists(paste0(k, "F-CV-", ml_method, "-performance.csv"))){
  file.remove(paste0(k, "F-CV-", ml_method, "-performance.csv"))
}
if(file.exists(paste0(gsub(".csv","",gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")))){
  file.remove(paste0(gsub(".csv","",gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")))
}

# Move features/genes with low variance:
# Default variance cutoff is 5-95%; i.e. if over 95% of individuals have the same value for the feature/gene, it is removed.
if(file.exists(paste0(gsub(".csv","",gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")))){
  df_input <- read.csv(paste0(gsub(".csv","",gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")), stringsAsFactors=F, check.names=F)
  df_input$status <- ifelse(df_input$status==1, "Positive", ifelse(df_input$status==0, "Negative", df_input$status))
  df_input$status <- as.factor(df_input$status)
}else{
  df <- fread(gs_fp, data.table=F)
  if(any(duplicated(df$Gene))){ # Remove the few duplicated genes based on protein lengths:
    prot_len <- fread(prot_len_fp, data.table=F)
    df_dups <- df[df$Gene %in% df$Gene[duplicated(df$Gene)], c(1,2)]
    df_dups$Prot_length <- prot_len$Prot_length[match(df_dups$Transcript, prot_len$Transcript)]
    df_dups_keep <- df_dups %>% group_by(Gene) %>% top_n(1, wt=Prot_length) %>% sample_n(1)
    transcript_to_remove <- setdiff(df_dups$Transcript, df_dups_keep$Transcript)
    
    df <- df[!df$Transcript%in%transcript_to_remove,]
    rownames(df) <- NULL
    rm(prot_len, df_dups, df_dups_keep, transcript_to_remove)
  }
  df <- df[, c("Gene", "Transcript", cvsch$sample_id)]
  colname_map <- df[, c("Gene", "Transcript")]
  # colname_map$col_name <- paste("V", c(1:nrow(df)), sep="")
  # write.table(colname_map, "map_colname_gene_transcript.txt", quote=F, col.names=T, sep="\t", row.names=F)
  df_fs <- as.data.frame(t(df[, 3:ncol(df)]))
  df_fs[is.na(df_fs)] <- 0
  colnames(df_fs) <- colname_map$Gene
  # Remove nzv:
  nzv <- nearZeroVar(df_fs, opt$variance_cutoff/(100-opt$variance_cutoff))
  df.X <- df_fs[, -nzv]
  df.y <- cvsch$Phenotype[match(rownames(df_fs), cvsch$sample_id)]
  df_input <- cbind(df.X, df.y)
  colnames(df_input)[ncol(df_input)] <- "status"
  df_input$status <- ifelse(df_input$status==1, "Positive", ifelse(df_input$status==0, "Negative", df_input$status))
  write.table(df_input, paste0(gsub(".csv", "", gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")), sep=",", quote=F, col.names=T, row.names=T)
  df_input$status <- as.factor(df_input$status)
}
file.rename(paste0(gsub(".csv", "", gs_fp), paste0(".NAto0.nzv", opt$variance_cutoff,"-", (100-opt$variance_cutoff), ".csv")), paste0(gsub(".csv", "", gs_fp), "_variation_filtered.csv"))

# Feature selection in a k-fold fashion:
if(fs_method %in% c("ks", "KS")){
  # Do KS test in a k-fold fashion:
  if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx"))){
    ks_FS_result <- read.xlsx(paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx"), 1)
    ks_FS_result$Gene <- as.character(ks_FS_result$Gene)
  }else if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"))){
    ks_FS_result <- read.csv(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"), stringsAsFactors=F, header=T)
  }else{
    # FS with KStest:
    ks_FS <- function(df){
      p_values <- c()
      for(i in 1:(ncol(df)-1)){ # (last col is status; remove)
        ks.res <- ks.test(df[df$status=="Positive", i], 
                          df[df$status=="Negative", i])
        pval <- ks.res$p.value
        p_values <- c(p_values, pval)
      }
      return(p_values)
    }
    
    ks_fs_parallel = function(i) {
      print(paste0("Performing k-fold KS FS: fold ", i, " ..."))
      return(ks_FS(df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],]))
    }
    ks_fs_res = mclapply(1:k, ks_fs_parallel, mc.cores = numCores)
    
    ks_FS_result <- as.data.frame(do.call("cbind", ks_fs_res))
    ks_FS_result <- cbind(colnames(df_input)[-ncol(df_input)], ks_FS_result)
    colnames(ks_FS_result) <- c("Gene", paste0("Pvalue_fold", 1:k, "out"))
    ks_FS_result$Gene <- as.character(ks_FS_result$Gene)
    
    write.csv(ks_FS_result, paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"), quote=F, row.names=F)
    print(paste0("KS FS results has been output to ", paste0(k, "F-CV-", fs_method, "-selectedGenes.csv")))
  }
}else if(fs_method %in% c("DKM", "DKMcost")){
  # Do CORElearn FS in a k-fold fashion:
  if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx"))){
    cl_fs_result <- read.xlsx(paste0(k, "F-CV-", fs_method,"-selectedGenes.xlsx"), 1)
    cl_fs_result$Gene <- as.character(cl_fs_result$Gene)
    cl_fs_result$Gene <- gsub("\\.", "-", cl_fs_result$Gene)
  }else if(file.exists(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"))){
    cl_fs_result <- read.csv(paste0(k, "F-CV-", fs_method,"-selectedGenes.csv"), stringsAsFactors=F, header=T)
    cl_fs_result$Gene <- gsub("\\.", "-", cl_fs_result$Gene)
  }else{
    suppressMessages(library(CORElearn))
    cl_fs_res <- list()
    for(i in 1:k){ # this may take a while ...
      print(paste0("Performing k-fold ", fs_method, " FS: fold ", i, " ..."))
      cl_fs <- attrEval(status ~ .,
                        df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],],
                        estimator=fs_method)
      cl_fs_res[[i]] <- cl_fs
      rm(cl_fs)
    }
    cl_fs_result <- do.call("cbind", cl_fs_res)
    cl_fs_result <- cbind(colnames(df_input)[-ncol(df_input)], cl_fs_result)
    colnames(cl_fs_result) <- c("Gene", paste0("Merit_fold", 1:k, "out"))
    cl_fs_result <- as.data.frame(cl_fs_result)
    # CORElearn automatically checkes column name; Needed to be converted back from "." to "-" for gene names:
    cl_fs_result$Gene <- gsub("\\.", "-", cl_fs_result$Gene)
    cl_fs_result$Gene <- as.character(cl_fs_result$Gene)
    
    #write.xlsx(cl_fs_result, paste0(k, "F-CV-", fs_method,"-selectedGenes.xlsx"), row.names=F)
    #print(paste0("CORElearn FS results has been output to ", paste0(k, "F-CV-", fs_method, "-selectedGenes.xlsx")))
    write.csv(cl_fs_result, paste0(k, "F-CV-", fs_method,"-selectedGenes.csv"), quote=F, row.names=F)
    print(paste0("CORElearn FS results has been output to ", paste0(k, "F-CV-", fs_method, "-selectedGenes.csv")))
  }
  
}

# Cross-validation with FS genes:
genes_considered = seq(opt$step_of_top_genes, opt$number_of_top_genes, opt$step_of_top_genes)
remaining = opt$number_of_top_genes - genes_considered[length(genes_considered)]
if (remaining > 0) {
  genes_considered[length(genes_considered) + 1] = opt$number_of_top_genes + remaining
}

genes_considered_parallel = function(gn) {
  if(gn>=2){
    print(paste0("Now running cross-validation with ", gn, " genes ..."))

    if(ml_method=="rf"){
      pred_res_rf <- list()
      suppressMessages(library(ranger))
      for(i in 1:k){
        
        train_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],]
        test_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% i],]
        
        w <- 1/table(train_dat$status)
        w <- w/sum(w)
        
        weights <- rep(0, nrow(train_dat))
        weights[train_dat$status=="Negative"] <- w["Negative"]
        weights[train_dat$status=="Positive"] <- w["Positive"]
        
        if(fs_method %in% c("ks")){

          genes_list = ks_FS_result[, paste0("Pvalue_fold", i, "out")]
          genes_ordered = order(genes_list)
          genes_subset = intersect(genes_ordered, which(genes_list < ks_pval_cutoff))
          genes_subset = genes_subset[1:min(gn, length(genes_subset))]
          gn_cutoff = min(gn, length(genes_subset))

          if (gn_cutoff == 0) {
            print(paste0("KS p-val cutoff (", ks_pval_cutoff, ") removed all genes in fold ", i, " ... aborting."))
            pred_res_rf[[i]] <- data.frame("status"=c(), "prediction"=c(), stringsAsFactors=F)
            rm(train_dat, test_dat, weights, w)
            next
          } else {
            removed_cnt = gn - length(genes_subset)
            print(paste0("KS p-val cutoff (", ks_pval_cutoff, ") removed ", removed_cnt, " (", format(round(100*removed_cnt/gn, 2), nsmall = 2), "%) of ", gn, " genes in fold ", i, " ..."))
          }
          
          model_rf <- ranger(dependent.variable.name="status",
                             data=train_dat[, c(ks_FS_result$Gene[genes_subset],
                                                "status")],
                             probability=T,
                             case.weights=weights
          )
          
          pred_rf <- predict(model_rf, 
                             test_dat[, c(ks_FS_result$Gene[genes_subset])])
          
        }else if(fs_method %in% c("DKM", "DKMcost")){
          model_rf <- ranger(dependent.variable.name="status",
                             data=train_dat[, c(cl_fs_result$Gene[order(cl_fs_result[, paste0("Merit_fold", i, "out")], decreasing=T)][1:gn],
                                                "status")],
                             probability=T,
                             case.weights=weights
          )
          
          pred_rf <- predict(model_rf, 
                             test_dat[, c(cl_fs_result$Gene[order(cl_fs_result[, paste0("Merit_fold", i, "out")], decreasing=T)][1:gn])])
          
        }
        
        pred_num <- pred_rf$predictions[, "Positive"]
        pred_res_rf[[i]] <- data.frame("status"=test_dat$status,
                                       "prediction"=pred_num,
                                       stringsAsFactors=F)
        rm(train_dat, test_dat, model_rf, pred_rf, pred_num, weights, w)
      }
      # Performance evaluation:
      pred_res_rf.df <- bind_rows(pred_res_rf)
      roc_rf <- roc.curve(pred_res_rf.df$prediction[pred_res_rf.df$status=="Positive"],
                          pred_res_rf.df$prediction[pred_res_rf.df$status=="Negative"])
      # print(paste0(gn, " genes: ", roc_rf$auc))
      return(roc_rf$auc)
      
    }else if(ml_method %in% c("svm", "SVM")){
      pred_res_svm <- list()
      suppressMessages(library(e1071))
      for(i in 1:k){

        train_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% setdiff(c(1:k), c(i))],]
        test_dat <- df_input[rownames(df_input) %in% cvsch$sample_id[cvsch$fold %in% i],]
        
        if(fs_method %in% c("ks")){

          genes_list = ks_FS_result[, paste0("Pvalue_fold", i, "out")]
          genes_ordered = order(genes_list)
          genes_subset = intersect(genes_ordered, which(genes_list < ks_pval_cutoff))
          genes_subset = genes_subset[1:min(gn, length(genes_subset))]
          gn_cutoff = min(gn, length(genes_subset))

          if (gn_cutoff == 0) {
            print(paste0("KS p-val cutoff (", ks_pval_cutoff, ") removed all genes in fold ", i, " ... aborting."))
            pred_res_svm[[i]] = data.frame("probabilities"=c(), "SampleID"=c(), stringsAsFactors=F)
            rm(train_dat, test_dat)
            next
          } else {
            removed_cnt = gn - length(genes_subset)
            print(paste0("KS p-val cutoff (", ks_pval_cutoff, ") removed ", removed_cnt, " (", format(round(100*removed_cnt/gn, 2), nsmall = 2), "%) of ", gn, " genes in fold ", i, " ..."))
          }

          model_svm <- svm(x=train_dat[, c(ks_FS_result$Gene[genes_subset])],
                           y=train_dat[, "status"],
                           probability=T,
                           class.weights = 100 / table(train_dat$status)
          )
          pred_svm <- predict(model_svm, 
                              test_dat[, c(ks_FS_result$Gene[genes_subset])],
                              decision.values=T, probability=T)
          
        }else if(fs_method %in% c("DKM", "DKMcost")){
          model_svm <- svm(x=train_dat[, c(cl_fs_result$Gene[order(cl_fs_result[, paste0("Merit_fold", i, "out")], decreasing=T)][1:gn])],
                           y=train_dat[, "status"],
                           probability=T,
                           class.weights = 100 / table(train_dat$status)
          )
          
          pred_svm <- predict(model_svm, 
                              test_dat[, c(cl_fs_result$Gene[order(cl_fs_result[, paste0("Merit_fold", i, "out")], decreasing=T)][1:gn])],
                              decision.values=T, probability=T)
          
        }
        #pred_num <- attr(pred_svm, "decision.values")
        pred_num <- attr(pred_svm, "probabilities")
        pred_res_svm[[i]] <- as.data.frame(pred_num)
        pred_res_svm[[i]]$SampleID <- rownames(pred_res_svm[[i]])
        rm(train_dat, test_dat, model_svm, pred_svm, pred_num)
      }
      # SVM performance evaluation:
      pred_res_svm.df <- bind_rows(pred_res_svm)
      pred_res_svm.df$status <- cvsch$Phenotype[match(pred_res_svm.df$SampleID, cvsch$sample_id)]
      roc_svm <- roc.curve(pred_res_svm.df$Positive[pred_res_svm.df$status==1],
                           pred_res_svm.df$Positive[pred_res_svm.df$status==0])
      # print(paste0(gn, " genes: ", roc_svm$auc))
      return(roc_svm$auc)
    }
  }
}

cv_performance_results = mclapply(genes_considered, genes_considered_parallel, mc.cores = numCores)
cv_performance_results.df <- as.data.frame(do.call("rbind", cv_performance_results))
colnames(cv_performance_results.df)[1] <- "AUC"

#write.xlsx(cv_performance_results.df, paste0(k, "F-CV-", ml_method, "-performance.xlsx"), sheetName=fs_method, append=T)

cv_performance_ordered = setDT(cv_performance_results.df, keep.rownames = TRUE)
cv_performance_ordered$rn = genes_considered
cv_performance_ordered = cv_performance_ordered[order(-cv_performance_ordered$AUC),]
setnames(cv_performance_ordered, "rn", "selected_genes")
write.table(cv_performance_ordered, paste0(k, "F-CV-", ml_method, "-performance.csv"), quote=F, row.names=F, col.names=T, sep=",")
for (ridx in 1:nrow(cv_performance_ordered)) {
  row = cv_performance_ordered[ridx,]
  gn = row[[1]]
  top_genes_by_auc = c()
  for(i in 1:k){
    if (fs_method %in% c("ks")) {
      fs_result = ks_FS_result$Gene[order(ks_FS_result[, paste0("Pvalue_fold", i, "out")])][1:gn]
    } else if (fs_method %in% c("DKM", "DKMcost")) {
      fs_result = cl_fs_result$Gene[order(cl_fs_result[, paste0("Merit_fold", i, "out")], decreasing=T)][1:gn]
    }
    top_genes_by_auc = c(top_genes_by_auc, fs_result)
  }
  top_genes_by_auc.df = as.data.frame(table(top_genes_by_auc))
  top_genes_by_auc.df = top_genes_by_auc.df[order(-top_genes_by_auc.df$Freq),]
  write.table(top_genes_by_auc[order(top_genes_by_auc)], paste0("AUC_rank.", ridx, "_", "top.", gn, "_", k, "F-CV-", fs_method,"-selectedGenes.csv"), quote=F, row.names=F, col.names=F)
  if (ridx == 1) {
    write.table(top_genes_by_auc.df$top_genes_by_auc, file.path(wd_fp, "AUC_rank.1-genes-list.csv"), quote=F, row.names=F, col.names=F)
    write.table(top_genes_by_auc.df, file.path(out_fp, "AUC_rank.1-genes.csv"), quote=F, row.names=F, col.names=F, sep = ",")
  }
}

file.copy(paste0(k, "F-CV-", ml_method, "-performance.csv"), file.path(out_fp, "performance.csv"), overwrite=T)
file.copy(paste0(k, "F-CV-", fs_method, "-selectedGenes.csv"), file.path(out_fp, "performance_folds_genes.csv"), overwrite=T)

print("Done.")
