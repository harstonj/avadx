#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

exonic.file <- args[1]
wd <- args[2]

exonic <- data.table::fread(exonic.file, header=F, data.table=F)

variants <- strsplit(exonic$V3, ",")
variants <- lapply(variants, function(x) strsplit(x, ":"))
variants_unlist <- unlist(unlist(variants))

all_mRNA <- variants_unlist[which(startsWith(variants_unlist, "exon"))-1]
all_mRNA <- unique(all_mRNA)


# For query upload at https://www.ncbi.nlm.nih.gov/sites/batchentrez
setwd(wd)

# Esearch:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=NM_152486,NM_015658,NM_198317

# Elink:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=protein&id=1519314899

# Efetch:
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=1519314899








write.table(all_mRNA[1:2000], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part1.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[2001:4500], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part2.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[4501:7500], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part3.txt", quote=F, row.names=F, col.names=F)
#Part3: The following records can't be retrieved:
#       Id=NM_016170
write.table(all_mRNA[7501:10500], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part4.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[10501:13000], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part5.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[13001:16000], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part6.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[16001:17000], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part7.txt", quote=F, row.names=F, col.names=F)
write.table(all_mRNA[17001:19000], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part8.txt", quote=F, row.names=F, col.names=F)

length(19001:46327)/2500
starts = seq(19001, 46327, 2500)
ends = seq(21500, 46327, 2500)

for(i in 1:10){
  write.table(all_mRNA[starts[i]:ends[i]], paste0("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part", i+8,".txt"), quote=F, row.names=F, col.names=F)
}
write.table(all_mRNA[44001:46327], "/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/All_mRNA_part19.txt", quote=F, row.names=F, col.names=F)




# Read in fasta files:
library(seqinr)
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/NCBI query/")
f1 <- read.fasta("prot_coding_sequence_part1.txt", as.string=T, seqtype="AA", set.attributes=F)
f1[[1]][1] #sequence
f1[[2]][1]
attr(f1[[1]])
tmp <- unlist(lapply(f1, function(y) y[[1]]))

dat <- list()
for(i in 1:19){
  dat[[i]] <- read.fasta(paste0("prot_coding_sequence_part",i,".txt"), as.string=T, seqtype="AA", set.attributes=F)
}

all_prot_names <- unlist(lapply(dat, names))
all_prot_sequences <- unlist(unlist(dat))
rm(dat)

variants[[1]]
variant_lengths <- lapply(variants, function(x) unlist(lapply(x, length)))
unique(unlist(variant_lengths))
for(i in 1:743201){
  if(1 %in% variant_lengths[[i]]){
    print(i)
  }
}

variants[[741781]]

#### Start writing SNAP input ####

# Start writing SNAP input:
exonic <- exonic[exonic$V3!="UNKNOWN",]
variants <- strsplit(exonic$V3, ",")
variants <- lapply(variants, function(x) strsplit(x, ":"))
variant_lengths <- lapply(variants, function(x) unlist(lapply(x, length)))
unique(unlist(variant_lengths)) # all is 5

# Decide which transcript to use for each gene:
# e.g.:
variants[962:963]
#TNFRSF18: NM_148901, NM_004195, NM_148902; all c.A67G and p.S23G
prot_idx <- grep(paste0("NM_148902", "\\."), all_prot_names)
prot <- all_prot_sequences[[prot_idx]]
prot 

variants_mRNA <- lapply(variants, function(x) lapply(x, function(y) y[2]))
variants_aaMut <- lapply(variants, function(x) lapply(x, function(y) y[5]))
variants_gene <- lapply(variants, function(x) lapply(x, function(y) y[1]))


variants_mRNA[962:963]

####### Find longest transcript for every gene #######
# Store a list of RefSeq transcripts for gene names:
gene_mRNA_list <- list()
for(i in 1:735146){
  g <- unlist(variants_gene[[i]])
  trans <- variants_mRNA_lov[[i]]
  dfi <- data.frame("Gene"=g, "Transcript"=trans)
  gene_mRNA_list[[i]] <- dfi
}
library(data.table)
gene_mRNA <- rbindlist(gene_mRNA_list, use.names=T)
gene_mRNA <- unique(gene_mRNA)
# Add length of the transcript to the df:
#gene_mRNA$Prot_length <- NA
for(i in 5873:nrow(gene_mRNA)){
  o <- gene_mRNA$Transcript[i]
  idx <- grep(paste0(o, "\\."), all_prot_names)
  #print(idx)
  if(length(idx) > 0){
    prot <- all_prot_sequences[[(idx)]]
    #print(prot)
    l <- nchar(prot)
    gene_mRNA$Prot_length[i] <- l
  }
}

gene_mRNA.cl <- gene_mRNA[!is.na(gene_mRNA$Prot_length),]
fwrite(gene_mRNA.cl, "../table-Transcript-ProtLength.csv")




####### Below: write SNAP input files #######
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part1") #1:100000
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part2") #100001:200000
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part3") #200001:300000
setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part4") #300001:400000
for(i in 300001:400000){
  print(i)
  mRNA <- variants_mRNA[i]
  mRNA <- unlist(mRNA)
  aaMut <- variants_aaMut[i]
  aaMut <- unlist(lapply(aaMut, function(x) gsub("p.", "", x)))
  
  #prot_idxs <- grep(paste(paste0(mRNA, "\\."), collapse="|"), all_prot_names)
  prot_idxs <- c()
  for(o in mRNA){
    idx <- grep(paste0(o, "\\."), all_prot_names)
    prot_idxs <- c(prot_idxs, idx)
  }
  
  if(length(prot_idxs)>0){
    # Finding the longest transcript:
    prot <- ""
    transcript <- ""
    mutation <- ""
    fasta_title <- ""
    for(j in 1:length(prot_idxs)){
      prot_new <- all_prot_sequences[[(prot_idxs[j])]]
      if(nchar(prot_new) > nchar(prot)){
        #print(j)
        prot <- prot_new
        transcript <- mRNA[j]
        mutation <- aaMut[j]
        fasta_title <- names(all_prot_sequences)[prot_idxs[j]]
      }
    }
    ref <- substr(mutation, 1, 1)
    alt <- substr(mutation, nchar(mutation), nchar(mutation))
    pos <- as.numeric(substr(mutation, 2, nchar(mutation)-1))
    
    # check if the mutation is correct:
    if((substr(prot, pos, pos)==ref) & (alt!="X")){
      write.table(mutation, paste0(transcript, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
    }
    if(!file.exists(paste0(transcript, ".fasta"))){
      write.fasta(prot, fasta_title, file.out=paste0(transcript, ".fasta"))
    }
  }
}

setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part5") #400001:500000
for(i in 400001:500000){
  print(i)
  mRNA <- variants_mRNA[i]
  mRNA <- unlist(mRNA)
  aaMut <- variants_aaMut[i]
  aaMut <- unlist(lapply(aaMut, function(x) gsub("p.", "", x)))
  
  #prot_idxs <- grep(paste(paste0(mRNA, "\\."), collapse="|"), all_prot_names)
  prot_idxs <- c()
  for(o in mRNA){
    idx <- grep(paste0(o, "\\."), all_prot_names)
    prot_idxs <- c(prot_idxs, idx)
  }
  
  if(length(prot_idxs)>0){
    # Finding the longest transcript:
    prot <- ""
    transcript <- ""
    mutation <- ""
    fasta_title <- ""
    for(j in 1:length(prot_idxs)){
      prot_new <- all_prot_sequences[[(prot_idxs[j])]]
      if(nchar(prot_new) > nchar(prot)){
        #print(j)
        prot <- prot_new
        transcript <- mRNA[j]
        mutation <- aaMut[j]
        fasta_title <- names(all_prot_sequences)[prot_idxs[j]]
      }
    }
    ref <- substr(mutation, 1, 1)
    alt <- substr(mutation, nchar(mutation), nchar(mutation))
    pos <- as.numeric(substr(mutation, 2, nchar(mutation)-1))
    
    # check if the mutation is correct:
    if((substr(prot, pos, pos)==ref) & (alt!="X")){
      write.table(mutation, paste0(transcript, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
    }
    if(!file.exists(paste0(transcript, ".fasta"))){
      write.fasta(prot, fasta_title, file.out=paste0(transcript, ".fasta"))
    }
  }
}

setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part6") #500001:600000
for(i in 500001:600000){
  print(i)
  mRNA <- variants_mRNA[i]
  mRNA <- unlist(mRNA)
  aaMut <- variants_aaMut[i]
  aaMut <- unlist(lapply(aaMut, function(x) gsub("p.", "", x)))
  
  #prot_idxs <- grep(paste(paste0(mRNA, "\\."), collapse="|"), all_prot_names)
  prot_idxs <- c()
  for(o in mRNA){
    idx <- grep(paste0(o, "\\."), all_prot_names)
    prot_idxs <- c(prot_idxs, idx)
  }
  
  if(length(prot_idxs)>0){
    # Finding the longest transcript:
    prot <- ""
    transcript <- ""
    mutation <- ""
    fasta_title <- ""
    for(j in 1:length(prot_idxs)){
      prot_new <- all_prot_sequences[[(prot_idxs[j])]]
      if(nchar(prot_new) > nchar(prot)){
        #print(j)
        prot <- prot_new
        transcript <- mRNA[j]
        mutation <- aaMut[j]
        fasta_title <- names(all_prot_sequences)[prot_idxs[j]]
      }
    }
    ref <- substr(mutation, 1, 1)
    alt <- substr(mutation, nchar(mutation), nchar(mutation))
    pos <- as.numeric(substr(mutation, 2, nchar(mutation)-1))
    
    # check if the mutation is correct:
    if((substr(prot, pos, pos)==ref) & (alt!="X")){
      write.table(mutation, paste0(transcript, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
    }
    if(!file.exists(paste0(transcript, ".fasta"))){
      write.fasta(prot, fasta_title, file.out=paste0(transcript, ".fasta"))
    }
  }
}

setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part7") #600001:700000
for(i in 600001:700000){
  print(i)
  mRNA <- variants_mRNA[i]
  mRNA <- unlist(mRNA)
  aaMut <- variants_aaMut[i]
  aaMut <- unlist(lapply(aaMut, function(x) gsub("p.", "", x)))
  
  #prot_idxs <- grep(paste(paste0(mRNA, "\\."), collapse="|"), all_prot_names)
  prot_idxs <- c()
  for(o in mRNA){
    idx <- grep(paste0(o, "\\."), all_prot_names)
    prot_idxs <- c(prot_idxs, idx)
  }
  
  if(length(prot_idxs)>0){
    # Finding the longest transcript:
    prot <- ""
    transcript <- ""
    mutation <- ""
    fasta_title <- ""
    for(j in 1:length(prot_idxs)){
      prot_new <- all_prot_sequences[[(prot_idxs[j])]]
      if(nchar(prot_new) > nchar(prot)){
        #print(j)
        prot <- prot_new
        transcript <- mRNA[j]
        mutation <- aaMut[j]
        fasta_title <- names(all_prot_sequences)[prot_idxs[j]]
      }
    }
    ref <- substr(mutation, 1, 1)
    alt <- substr(mutation, nchar(mutation), nchar(mutation))
    pos <- as.numeric(substr(mutation, 2, nchar(mutation)-1))
    
    # check if the mutation is correct:
    if((substr(prot, pos, pos)==ref) & (alt!="X")){
      write.table(mutation, paste0(transcript, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
    }
    if(!file.exists(paste0(transcript, ".fasta"))){
      write.fasta(prot, fasta_title, file.out=paste0(transcript, ".fasta"))
    }
  }
}

setwd("/Users/WangYanran/Documents/BrombergLab/AVA_Method_paper/Filtered_Data/dbGaP/SNAP_scores/SNAP_input_part8") #700001:length(variants)
for(i in 700001:length(variants)){
  print(i)
  mRNA <- variants_mRNA[i]
  mRNA <- unlist(mRNA)
  aaMut <- variants_aaMut[i]
  aaMut <- unlist(lapply(aaMut, function(x) gsub("p.", "", x)))
  
  #prot_idxs <- grep(paste(paste0(mRNA, "\\."), collapse="|"), all_prot_names)
  prot_idxs <- c()
  for(o in mRNA){
    idx <- grep(paste0(o, "\\."), all_prot_names)
    prot_idxs <- c(prot_idxs, idx)
  }
  
  if(length(prot_idxs)>0){
    # Finding the longest transcript:
    prot <- ""
    transcript <- ""
    mutation <- ""
    fasta_title <- ""
    for(j in 1:length(prot_idxs)){
      prot_new <- all_prot_sequences[[(prot_idxs[j])]]
      if(nchar(prot_new) > nchar(prot)){
        #print(j)
        prot <- prot_new
        transcript <- mRNA[j]
        mutation <- aaMut[j]
        fasta_title <- names(all_prot_sequences)[prot_idxs[j]]
      }
    }
    ref <- substr(mutation, 1, 1)
    alt <- substr(mutation, nchar(mutation), nchar(mutation))
    pos <- as.numeric(substr(mutation, 2, nchar(mutation)-1))
    
    # check if the mutation is correct:
    if((substr(prot, pos, pos)==ref) & (alt!="X")){
      write.table(mutation, paste0(transcript, ".mutation"), quote=F, row.names=F, col.names=F, append=T)
    }
    if(!file.exists(paste0(transcript, ".fasta"))){
      write.fasta(prot, fasta_title, file.out=paste0(transcript, ".fasta"))
    }
  }
}

######
# Check for error after SNAP run:
# NM_004759.fasta has title: >lcl|NM_032960.4_prot_NP_116584.2_1
grep("NM_004759\\.", all_prot_names)
all_prot_sequences[4288]

#tmp <- exonic[grep("NM_004759:", exonic$V3),]
tmp <- exonic[grep("MAPKAPK2:", exonic$V3),]

# Reasons:
# NM_004759 and NM_032960 are from the same gene
# NM_032960 is longer than NM_004759
# So, files NM_004759.fasta and NM_032960.fasta both contain the exact same sequence for NM_032960.4_prot_NP_116584.2_1
#




