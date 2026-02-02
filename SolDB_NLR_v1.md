---
title: "SolDB_NLR_v1"
author: "AmirAli Toghani"
date: "2025-07-21"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r}
library(tidyverse)
library(readxl)
library(Biostrings)
library(ggtree)
library(pheatmap)
library(reshape2)
library(svglite)
library(wordcloud2)
```

## 1. import the species/genome metadata:
```{r}
metadata <- read_xlsx("/path/to/Data_SX1.xlsx")

species_meta <- metadata[,c(2,3,4,5,6,14,20)]
```


## 2. Import sequences and metadata:
```{r}
setwd("/path/to/nlrtracker")

nlrtracker_files <- list.files(pattern = "_Domains.tsv", recursive = T)

nlrtracker_list <- list()

# Loop through each file
for (i in seq_along(nlrtracker_files)) {
  # Read the Excel file
  nlrtracker_list[[i]] <- read_delim(nlrtracker_files[i])
  
  # Extract the specific part of the file name and add it as a new column
  nlrtracker_list[[i]] <- nlrtracker_list[[i]] %>%
    mutate('Proteome File Name' = sub("_Domains\\.tsv$", ".gz", sub("^nlrtracker_", "", basename(nlrtracker_files[i]))))

  # Print a message for each file (optional)
  cat("Read file:", nlrtracker_files[i], "and added formatted filename column\n")
}


# merge all and add the metadata
NLRtracker <- do.call(rbind,nlrtracker_list) %>% inner_join(species_meta, by = "Proteome File Name")


NLRtracker <- NLRtracker %>% 
  mutate(ID = paste(`Assembly Accession`, Species, seqname, Simple, sep = "_")) %>%
  # Replace spaces in the "Species" column with underscores for the "ID" column
  mutate(ID = gsub(" ", "_", ID))

write_csv(NLRtracker, "/path/to/NLRtracker.csv")

```

## 2. Prepare the functions and prerequisites:
```{r}

domains <- c("CNL","CNLO","CN","OCNL","CONL",
             "NL","NLO","ONL",
             "BCNL","BCNLO","BNL","BBNL","BCN","BCCNL","BBCNL","BNLO","BOCNL","BBCNLO",
             "RNL",
             "TN","TNL","TNLO","TNLJ")

```


## 4. Extract the NLRs and NBARCs:
```{r}
NLR <- filter(NLRtracker, NLRtracker$type == "CHAIN" & NLRtracker$Status == "NLR")
NLR_seqnames <- NLR[,c(1,17)]

# only keep the ones with desired domain architechture
NLR_filtered <- NLR[NLR$Simple %in% domains,]

NLR_filtered_seq <- AAStringSet(NLR_filtered$sequence)
NLR_filtered_seq@ranges@NAMES <- NLR_filtered$ID

writeXStringSet(NLR_filtered_seq, "/path/to/NLR_filtered_seq.fasta")


# remove duplicate sequences using CD-HIT
NLR_filtered_deduplicated_seq <- readAAStringSet("/path/to/NLR_filtered_seq_clu.fasta")

NLR_filtered_deduplicated <- NLR_filtered[NLR_filtered$ID %in% NLR_filtered_deduplicated_seq@ranges@NAMES,]

# extract the NBARC domains
NBARC <- NLRtracker[NLRtracker$description == "NBARC",]
write_csv(NBARC, "/path/to/NBARC.csv")

NBARC_filtered <- NBARC[NBARC$ID %in% NLR_filtered$ID,]
NBARC_filtered_seq <- AAStringSet(NBARC_filtered$sequence)
NBARC_filtered_seq@ranges@NAMES <- NBARC_filtered$ID

NBARC_filtered_deduplicated <- NBARC[NBARC$ID %in% NLR_filtered_deduplicated$ID,]

# filter out NLRs with truncated NBARC domain
NBARC_filtered_deduplicated_len <- filter(NBARC_filtered_deduplicated, 
                                                        NBARC_filtered_deduplicated$end -
                                                        NBARC_filtered_deduplicated$start > 250 &
                                                          NBARC_filtered_deduplicated$end -
                                                          NBARC_filtered_deduplicated$start < 400) # NBARCs



NLR_filtered_deduplicated_len <- NLR_filtered_deduplicated[NLR_filtered_deduplicated$ID %in% NBARC_filtered_deduplicated_len$ID,] # NLRs


# convert the final data to biostring objects
NLR_filtered_deduplicated_len_seq <- AAStringSet(NLR_filtered_deduplicated_len$sequence)
NLR_filtered_deduplicated_len_seq@ranges@NAMES <- NLR_filtered_deduplicated_len$ID

NBARC_filtered_deduplicated_len_seq <- AAStringSet(NBARC_filtered_deduplicated_len$sequence)
NBARC_filtered_deduplicated_len_seq@ranges@NAMES <- NBARC_filtered_deduplicated_len$ID

```

## 5. Import RefPlantNLR and add it to the main dataset and export for phylogenetics analysis:
```{r}
RefPlantNLR <- readAAStringSet("/path/to/RefPlantNLR.fasta")

RefPlantNLR_NBARC <- readAAStringSet("/path/to/RefPlantNLR_NBARC.fasta")

NLR_filtered_deduplicated_len_seq_ref <- c(NLR_filtered_deduplicated_len_seq, RefPlantNLR)
writeXStringSet(NLR_filtered_deduplicated_len_seq_ref, "/path/to/NLR_ref.fasta")

NBARC_filtered_deduplicated_len_seq_ref <- c(NBARC_filtered_deduplicated_len_seq, RefPlantNLR_NBARC)
writeXStringSet(NBARC_filtered_deduplicated_len_seq_ref, "/path/to/NBARC_ref.fasta")



NBARC_filtered_seq_ref <- c(NBARC_filtered_seq, RefPlantNLR_NBARC)
writeXStringSet(NBARC_filtered_seq_ref, "/path/to/NBARC_filtered_seq_ref.fasta")

```