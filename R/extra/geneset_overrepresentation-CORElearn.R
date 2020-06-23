suppressMessages(library(data.table))
suppressMessages(library(xlsx))

args <- commandArgs(trailingOnly=TRUE)

fp <- args[1]   # path to folder of CORElearn results
gn <- as.numeric(args[2])   # number of genes to select for analysis
fp_bgdf <- args[3]  # path to the background df e.g. GeneScoreTable_normed.NAto0.nzv85-15.txt
fp_map_gene <- args[4] # path to the map_colname_gene_transcript.txt file
method <- args[5]
fp_out <- args[6]

setwd(fp)

bg.df <- read.csv(fp_bgdf, stringsAsFactors=F, nrows=5)
bg.colnames <- colnames(bg.df)

map_col_gene <- read.table(fp_map_gene, stringsAsFactors=F, header=T)
background_genes <- map_col_gene$Gene[map_col_gene$col_name%in%bg.colnames]


if(method=="DKM"){
  DKM <- list()
  for(i in 1:10){
    DKM[[i]] <- readRDS(paste0("DKM/fs_DKM_fold",i,".rds"))
  }
  FS_result <- data.frame("Gene"=map_col_gene$Gene[match(names(DKM[[1]]), map_col_gene$col_name)],
                          "Transcript"=map_col_gene$Transcript[match(names(DKM[[1]]), map_col_gene$col_name)],
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
}else if(method=="DKMc"){
  DKMc <- list()
  for(i in 1:10){
    DKMc[[i]] <- readRDS(paste0("DKMc/fs_DKMcost_fold",i,".rds"))
  }
  FS_result <- data.frame("Gene"=map_col_gene$Gene[match(names(DKMc[[1]]), map_col_gene$col_name)],
                          "Transcript"=map_col_gene$Transcript[match(names(DKMc[[1]]), map_col_gene$col_name)],
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
}else if(method=="ReliefF"){
  Relieff <- list()
  for(i in 1:10){
    Relieff[[i]] <- readRDS(paste0("ReliefFexpRank/fs_ReliefFexpRank_fold", i, ".rds"))
  }
  FS_result <- data.frame("Gene"=map_col_gene$Gene[match(names(Relieff[[1]]), map_col_gene$col_name)],
                          "Transcript"=map_col_gene$Transcript[match(names(Relieff[[1]]), map_col_gene$col_name)],
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


k <- 10

genes_all <- list()
for(i in 1:k){
  genes <- FS_result$Gene[order(-FS_result[,(2+i)])][1:gn] # large to small (rank merit)
  genes_all[[i]] <- genes
}
genes_inter <- Reduce("intersect", genes_all)
genes_union <- Reduce("union", genes_all)


# Import CPDB:
cpdb.gene_symbol <- fread("/Users/WangYanran/Documents/Bitbucket/repos/avadx-meta/db/cpdb/CPDB_pathways_genesymbol.tab",
                          stringsAsFactors=F)

cpdb_genes <- strsplit(cpdb.gene_symbol$hgnc_symbol_ids, ",")
cpdb_genes_all <- unique(unlist(cpdb_genes))

cpdb_analysis <- function(input.fs.result, # a data frame containing gene and transcript names and their FS results; see below
                          i=1, # the column number for FS cross-validation result; 1 means fold 1 taken out
                          gene.number, # number of genes to use
                          descending=F, # whether to rank the FS result by descending order; defult is FALSE for KStest results where p-vals will be ranked ascendingly
                          fdr.cutoff=0.05){
  # input.fs.result is a data.frame with 1st col "Gene" and 2nd col "Transcript"
  # the FS result (e.g. p-val from K-S test or merit from CORElearn attrEval) starts with the 3rd col
  # For 10-fold cv, the FS results should be from Col 3 to Col 12; for 20-fold cv, the FS results should be from Col 3 to Col 22
  
  genes <- input.fs.result$Gene[order(input.fs.result[,(2+i)])][1:gene.number]
  if(descending==T){
    genes <- input.fs.result$Gene[order(-input.fs.result[,(2+i)])][1:gene.number]
  }
  genes_intersect.cpdb <- lapply(cpdb_genes, function(x) intersect(x, genes))
  genes_intersect.cpdb_length <- lapply(cpdb_genes, function(x) length(intersect(x, genes)))
  cpdb.gene_symbol$pathway[which(unlist(genes_intersect.cpdb_length)>2)]
  
  bg_cpdb_overlap <- intersect(background_genes, cpdb_genes_all)
  
  # phyper: 
  #  white balls: genes from pathway X
  #  black balls: genes not from pathway X
  #  draw: FS
  genes_intersect.cpdb_pval <- lapply(cpdb_genes, function(x) phyper(q=length(intersect(x, genes)), # intersect of pathway X genes and FS gene set
                                                                     m=length(intersect(x, background_genes)), # (# of white balls) # of genes from pathway X
                                                                     n=length(bg_cpdb_overlap)-length(intersect(x, background_genes)),  # (# of black balls) # of genes not from pathway X
                                                                     k=length(intersect(genes, cpdb_genes_all)), # (# of drawn) size of the FS gene set
                                                                     lower.tail=F))
  
  # adjust p-value:
  fdr_pval <- p.adjust(unlist(genes_intersect.cpdb_pval), method="fdr")
  #bonfer_pval <- p.adjust(unlist(genes_intersect.cpdb_pval), method="bonferroni")
  
  cpdb_result <- cpdb.gene_symbol
  cpdb_result$pathway.size <- unlist(lapply(cpdb_genes, function(x) length(x)))
  cpdb_result$pathway.size.in.background <- unlist(lapply(cpdb_genes, function(x) length(intersect(x, background_genes))))
  cpdb_result$fdr.pval <- fdr_pval
  #cpdb_result$bonferroni.pval <- bonfer_pval
  cpdb_result$num_overlap <- unlist(genes_intersect.cpdb_length)
  cpdb_result$fdr.pval <- ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$fdr.pval)
  #cpdb_result$bonferroni.pval <- ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$bonferroni.pval)
  cpdb_result$overlapping.genes <- unlist(lapply(cpdb_genes, function(x) paste(intersect(x, genes), collapse=",")))
  
  cpdb_result_significant <- cpdb_result[cpdb_result$fdr.pval<fdr.cutoff,]
  cpdb_result_significant$hgnc_symbol_ids <- NULL
  
  return(cpdb_result_significant)
}

# cpdb_analysis(FS_result, i=1, gene.number=gn, fdr.cutoff=.05)

cpdb_results <- list()
for(i in 1:k){
  cpdb_results[[i]] <- cpdb_analysis(FS_result, i=i, gene.number=gn, fdr.cutoff=.05, descending=T)
}
pathway_overlap_cv <- as.data.frame(table(unlist(lapply(cpdb_results, function(x) x[, "pathway"]))))

for(i in 1:k){
  write.xlsx(cpdb_results[[i]], fp_out, sheetName=paste0("Fold", i, "_out"), append=T)
}

write.xlsx(pathway_overlap_cv, fp_out, sheetName="Overlapped", append=T)



