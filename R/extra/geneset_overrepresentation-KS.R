suppressMessages(library(data.table))
suppressMessages(library(xlsx))

args <- commandArgs(trailingOnly=TRUE)

fp <- args[1]   # path to 10F-CV-KS-selectedGenes.xlsx
gn <- as.numeric(args[2])   # number of genes to select for analysis
fp_bgdf <- args[3]  # path to the background df e.g. GeneScoreTable_normed.NAto0.nzv85-15.txt
fp_map_gene <- args[4] # path to the map_colname_gene_transcript.txt file
fp_out <- args[5]

FS_result <- read.xlsx(fp, 1, colIndex=c(2:13))
FS_result$Gene <- as.character(FS_result$Gene)
FS_result$Transcript <- as.character(FS_result$Transcript)

bg.df <- read.csv(fp_bgdf, stringsAsFactors=F, nrows=5)
bg.colnames <- colnames(bg.df)

map_col_gene <- read.table(fp_map_gene, stringsAsFactors=F, header=T)
background_genes <- map_col_gene$Gene[map_col_gene$col_name%in%bg.colnames]


k <- 10

genes_all <- list()
for(i in 1:k){
  genes <- FS_result$Gene[order(FS_result[,(2+i)])][1:gn] # small to large (rank p-val)
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
  cpdb_results[[i]] <- cpdb_analysis(FS_result, i=i, gene.number=gn, fdr.cutoff=.05)
}
pathway_overlap_cv <- as.data.frame(table(unlist(lapply(cpdb_results, function(x) x[, "pathway"]))))

for(i in 1:k){
  write.xlsx(cpdb_results[[i]], fp_out, sheetName=paste0("Fold", i, "_out"), append=T)
}

write.xlsx(pathway_overlap_cv, fp_out, sheetName="Overlapped", append=T)



