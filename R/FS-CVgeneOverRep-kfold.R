suppressMessages(library(xlsx))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
library(parallel)

option_list = list(
  make_option(c("-f", "--input_file"), type="character", default=NULL,
              help="path to input file; e.g. selectedGenes.csv", metavar="character"),
  make_option(c("-b", "--input_file_genescore"), type="character", default=NULL,
              help="path to the cleaned gene score table; e.g. GeneScoreTable_normalized.NAto0.nzv85-15.txt", metavar="character"),
  make_option(c("-t", "--top_genes"), type="character", default=NULL,
              help="path to the top AUC rank file", metavar="character"),
  make_option(c("-n", "--number_of_top_genes"), type="numeric", default=100,
              help="number of top-ranked genes for over-representation analysis; e.g. 100 (default) means using the top 100 genes as over-representation analysis input", metavar="character"),
  make_option(c("-d", "--cpdb_database"), type="character", default=NULL,
              help="path to the CPDB_pathways_genesymbol.tab file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="path to output folder", metavar="character")
)
opt = parse_args(OptionParser(option_list=option_list))
numCores = detectCores()

# change to output directory
setwd(opt$out)

# Read in CPDB:
cpdb = fread(opt$cpdb_database, data.table=F)
cpdb_genes = strsplit(cpdb$hgnc_symbol_ids, ",")
cpdb_genes_all = unique(unlist(cpdb_genes))

# Read in FS result:
fs = fread(opt$input_file, stringsAsFactors=F)

# Read AUC top genes list
top_genes = fread(opt$top_genes, header=F, stringsAsFactors=F)

# check topgenes selection
if (opt$number_of_top_genes == 0) {
  opt$number_of_top_genes = nrow(top_genes)
}

# Read in gene score table for all background genes:
df = fread(opt$input_file_genescore, stringsAsFactors=F, check.names=F)

# Background genes
background_genes = colnames(df)[-1]
bg_cpdb_overlap = intersect(background_genes, cpdb_genes_all)

# update cpdb
cpdb$pathway.size = unlist(lapply(cpdb_genes, function(x) length(x)))
cpdb$pathway.size.in.background = unlist(lapply(cpdb_genes, function(x) length(intersect(x, background_genes))))

fdr.cutoff = 0.05
cpdb_analysis_parallel = function(i) {
 
  if (i > 0) {
    print(paste0("Performing enrichment analysis for fold ", i, " ..."))
    fold = c(as.character(1))
    genes = fs[[fold]][1:opt$number_of_top_genes]
  } else {
    print("Performing enrichment analysis for rank 1 AUC ...")
    genes = top_genes[[1]]
  }
  genes = genes[!is.na(genes)]

  genes_intersect.cpdb = lapply(cpdb_genes, function(x) intersect(x, genes))
  genes_intersect.cpdb_length = lapply(cpdb_genes, function(x) length(intersect(x, genes)))
  # cpdb$pathway[which(unlist(genes_intersect.cpdb_length)>2)]
  
  # phyper: 
  #  white balls: genes from pathway X
  #  black balls: genes not from pathway X
  #  draw: FS
  genes_intersect.cpdb_pval = lapply(cpdb_genes, function(x) phyper(q=length(intersect(x, genes)), # intersect of pathway X genes and FS gene set
                                                                     m=length(intersect(x, background_genes)), # (# of white balls) # of genes from pathway X
                                                                     n=length(bg_cpdb_overlap)-length(intersect(x, background_genes)),  # (# of black balls) # of genes not from pathway X
                                                                     k=length(intersect(genes, cpdb_genes_all)), # (# of drawn) size of the FS gene set
                                                                     lower.tail=F))
  
  # adjust p-value:
  fdr_pval = p.adjust(unlist(genes_intersect.cpdb_pval), method="fdr")
  #bonfer_pval = p.adjust(unlist(genes_intersect.cpdb_pval), method="bonferroni")
  
  cpdb_result = copy(cpdb)
  cpdb_result$fdr.pval = fdr_pval
  #cpdb_result$bonferroni.pval = bonfer_pval
  cpdb_result$num_overlap = unlist(genes_intersect.cpdb_length)
  cpdb_result$fdr.pval = ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$fdr.pval)
  #cpdb_result$bonferroni.pval = ifelse(cpdb_result$num_overlap%in%c(0, 1), NA, cpdb_result$bonferroni.pval)
  cpdb_result$overlapping.genes = unlist(lapply(cpdb_genes, function(x) paste(intersect(x, genes), collapse=",")))
  
  cpdb_result_significant = cpdb_result[which(cpdb_result$fdr.pval<fdr.cutoff),]
  #cpdb_result_significant$hgnc_symbol_ids = NULL
  
  return(cpdb_result_significant)
}

cpdb_overrep = mclapply(colnames(fs)[-1], cpdb_analysis_parallel, mc.cores=numCores)

cpdb_overrep_pathway = as.data.frame(table(unlist(lapply(cpdb_overrep, function(x) x[, "pathway"]))))

if (nrow(cpdb_overrep_pathway) > 0) {
  colnames(cpdb_overrep_pathway) = c("PathwayName", "Freq")
} else {
  cpdb_overrep_pathway = data.frame(PathwayName=NA, Freq=NA)
}

dir.create(file.path('PathwayOverRepresentation_by_folds'), recursive = TRUE, showWarnings = FALSE)
# Output k-fold details:
for(i in 1:(ncol(fs)-1)){
  cpdb_overrep_df_i = data.frame(cpdb_overrep[[i]])
  if (nrow(cpdb_overrep_df_i) == 0) {
    cpdb_overrep_i = cpdb_overrep_df_i
    cpdb_overrep_i[nrow(cpdb_overrep_i)+1,] = NA
  } else {
    cpdb_overrep_i = cpdb_overrep[[i]]
  }
  cpdb_overrep_dt_i = data.table(cpdb_overrep_i)
  cpdb_overrep_dt_i[, overlapping.genes:=gsub(overlapping.genes, pattern=",", replacement=";")]
  cpdb_overrep_dt_i[, hgnc_symbol_ids:=gsub(hgnc_symbol_ids, pattern=",", replacement=";")]
  write.table(cpdb_overrep_dt_i, 
             file.path("PathwayOverRepresentation_by_folds", paste0("fold_", i, ".csv")),
             quote=F, row.names=F, col.names=T, sep=",")
}

if(file.exists(paste0("PathwayOverRepresentation_by_folds_summary.csv"))){
  file.remove(paste0("PathwayOverRepresentation_by_folds_summary.csv"))
}
# Output the summary:
write.table(cpdb_overrep_pathway,
           paste0("PathwayOverRepresentation_by_folds_summary.csv"),
           quote=F, row.names=F, col.names=T, sep=",")

# OUTPUT AUC rank.1
auc_rank1 = data.table(cpdb_analysis_parallel(0))
auc_rank1[, overlapping.genes:=gsub(overlapping.genes, pattern=",", replacement=";")]
auc_rank1[, hgnc_symbol_ids:=gsub(hgnc_symbol_ids, pattern=",", replacement=";")]
write.table(auc_rank1, "PathwayOverRepresentation_AUC_rank.1-genes.csv", quote=F, row.names=F, col.names=T, sep=",")
