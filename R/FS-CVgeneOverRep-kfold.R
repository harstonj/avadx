suppressMessages(library(xlsx))
suppressMessages(library(optparse))
suppressMessages(library(data.table))

option_list = list(
  make_option(c("-f", "--input_file"), type="character", default=NULL,
              help="path to input file; e.g. 10F-CV-ks-selectedGenes.xlsx", metavar="character"),
  make_option(c("-b", "--input_file_genescore"), type="character", default=NULL,
              help="path to the cleaned gene score table; e.g. GeneScoreTable_normed.NAto0.nzv85-15.txt", metavar="character"),
  make_option(c("-n", "--number_of_top_genes"), type="numeric", default=100,
              help="number of top-ranked genes for over-representation analysis; e.g. 100 (default) means using the top 100 genes as over-representation analysis input", metavar="character"),
  make_option(c("-d", "--cpdb_database"), type="character", default=NULL,
              help="path to the CPDB_pathways_genesymbol.tab file", metavar="character"),
  make_option(c("-a", "--ascending"), type="logical", default=T,
              help="T: list gene ranking ascendingly (from small to large; for ks); F: list descendingly (from large to small; for DKM)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="path to output folder", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))

print(opt$input_file)
print(opt$number_of_top_genes)
print(opt$cpdb_database)
print(opt$out)

setwd(opt$out)

# Read in CPDB:
cpdb <- fread(opt$cpdb_database, data.table=F)
cpdb_genes <- strsplit(cpdb$hgnc_symbol_ids, ",")
cpdb_genes_all <- unique(unlist(cpdb_genes))

# Read in FS result:
fs <- read.xlsx(opt$input_file, 1, stringsAsFactors=F)

# Read in gene score table for all background genes:
df <- read.csv(opt$input_file_genescore, stringsAsFactors=F, check.names=F)


cpdb_analysis <- function(input.fs.result, # a data frame containing gene and transcript names and their FS results; see below
                          i=1, # the column number for FS cross-validation result; 1 means fold 1 taken out
                          gene.number, # number of genes to use
                          descending=F, # whether to rank the FS result by descending order; defult is FALSE for KStest results where p-vals will be ranked ascendingly
                          fdr.cutoff=0.05){
  # input.fs.result is a data.frame with 1st col "Gene" and 2nd col "Transcript"
  # the FS result (e.g. p-val from K-S test or merit from CORElearn attrEval) starts with the 3rd col
  # For 10-fold cv, the FS results should be from Col 3 to Col 12; for 20-fold cv, the FS results should be from Col 3 to Col 22
  
  genes <- input.fs.result$Gene[order(input.fs.result[,(1+i)])][1:gene.number]
  if(descending==T){
    genes <- input.fs.result$Gene[order(-input.fs.result[,(1+i)])][1:gene.number]
  }
  genes_intersect.cpdb <- lapply(cpdb_genes, function(x) intersect(x, genes))
  genes_intersect.cpdb_length <- lapply(cpdb_genes, function(x) length(intersect(x, genes)))
  cpdb$pathway[which(unlist(genes_intersect.cpdb_length)>2)]
  
  background_genes <- colnames(df)[-ncol(df)]
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
  
  cpdb_result <- cpdb
  cpdb_result$pathway.size <- unlist(lapply(cpdb_genes, function(x) length(x)))
  cpdb_result$pathway.size.in.background <- unlist(lapply(cpdb_genes, function(x) length(intersect(x, background_genes))))
  cpdb_result$fdr.pval <- fdr_pval
  #cpdb_result$bonferroni.pval <- bonfer_pval
  cpdb_result$num_overlap <- unlist(genes_intersect.cpdb_length)
  cpdb_result$fdr.pval <- ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$fdr.pval)
  #cpdb_result$bonferroni.pval <- ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$bonferroni.pval)
  cpdb_result$overlapping.genes <- unlist(lapply(cpdb_genes, function(x) paste(intersect(x, genes), collapse=",")))
  
  cpdb_result_significant <- cpdb_result[which(cpdb_result$fdr.pval<fdr.cutoff),]
  #cpdb_result_significant$hgnc_symbol_ids <- NULL
  
  return(cpdb_result_significant)
}

cpdb_overrep <- list()
for(i in 1:(ncol(fs)-1)){
  cpdb_overrep[[i]] <- cpdb_analysis(fs, i=i, gene.number=opt$number_of_top_genes, descending=!(opt$ascending))
}

cpdb_overrep_pathway <- as.data.frame(table(unlist(lapply(cpdb_overrep, function(x) x[, "pathway"]))))

colnames(cpdb_overrep_pathway) <- c("PathwayName", "Freq") 

# Output k-fold details:
for(i in 1:(ncol(fs)-1)){
  write.xlsx(cpdb_overrep[[i]], 
             paste0("PahtwayOverRepresentation_GeneN-", opt$number_of_top_genes, "_", gsub("-selectedGenes.xlsx", "",basename(opt$input_file)), ".xlsx"), 
             sheetName=paste0("Fold", i, "out"), append=T)
}

# Output the summary:
write.xlsx(cpdb_overrep_pathway,
           paste0("PahtwayOverRepresentation_GeneN-", opt$number_of_top_genes, "_", gsub("-selectedGenes.xlsx", "",basename(opt$input_file)), "_summary.xlsx")
           )




