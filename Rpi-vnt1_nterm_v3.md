---
title: "Rpi-vnt1_nterm_v3"
author: "AmirAli Toghani"
date: "2025-07-21"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r}
library(tidyverse)
library(Biostrings)
library(readxl)
library(svglite)
```

## 0. Extract the CC domains:
```{r}
CC <- NLRtracker %>% filter(type == "REGION" & description == "CC")
```


## 1. Import the HMMsearch output for MADA:
```{r}
# SolDB
# import the hmm output for MADA in SolDB database
MADA_hmm_SolDB_raw <- read.table("/path/to/MADA/MADA_tbl.out", skip = 3) # table output
MADA_hmm_SolDB_dom_raw <- read.table("/path/to/MADA_domtbl.out", skip = 3) # domain table output

# clean up the domain table
MADA_hmm_SolDB_dom <- MADA_hmm_SolDB_dom_raw[,c(1,3,7,8,18,19)] %>% setNames(c("seqname","length","evalue","score","from","to"))
MADA_hmm_SolDB_dom$domainL <- (MADA_hmm_SolDB_dom$to - MADA_hmm_SolDB_dom$from + 1) # add a domain length column

MADA_SolDB_meta <- NLR_filtered_deduplicated_len[NLR_filtered_deduplicated_len$ID %in% MADA_hmm_SolDB_dom$seqname,]

MADA_SolDB_ext <- MADA_hmm_SolDB_dom[MADA_hmm_SolDB_dom$from >= 2,] # anything with MADA starting from the second residue
MADA_SolDB_ext_seq <- NLR_filtered_deduplicated_len_seq_ref[MADA_SolDB_ext$seqname]

MADA_SolDB_ext_startM_seq <- subseq(MADA_SolDB_ext_seq, start = MADA_SolDB_ext$from, end = MADA_SolDB_ext$from)
MADA_SolDB_ext_startM_df <- as.data.frame(MADA_SolDB_ext_startM_seq) %>% setNames("MADA starting AA")
MADA_SolDB_ext_startM_df$ID <- MADA_SolDB_ext_startM_seq@ranges@NAMES


```

```{r}
# N-terminus through HMM start (INCLUSIVE)
MADA_SolDB_leader <- subseq(MADA_SolDB_ext_seq, start = 1, end = MADA_SolDB_ext$from)

# counts & flags
first_AA          <- subseq(MADA_SolDB_ext_seq, start = 1, end = 1)
first_is_M        <- as.vector(first_AA == "M")

m_counts_mat      <- letterFrequency(MADA_SolDB_leader, letters = "M", as.prob = FALSE)
m_count           <- as.integer(m_counts_mat[, 1])

domstart_AA       <- subseq(MADA_SolDB_ext_seq,
                            start = MADA_SolDB_ext$from,
                            end   = MADA_SolDB_ext$from) %>% as.vector()
M_at_domstart     <- domstart_AA == "M"

# annotate the domain table with sequence-derived features
MADA_SolDB_ext_annot <- MADA_SolDB_ext %>%
  mutate(
    first_is_M     = first_is_M,
    m_count        = m_count,              # number of 'M' in positions 1..from (inclusive)
    domstart_AA    = domstart_AA,
    M_at_domstart  = M_at_domstart
  )

# keep: starts with M and has a second M somewhere in 1..from (inclusive)
MADA_SolDB_keep <- MADA_SolDB_ext_annot %>%
  filter(first_is_M, m_count >= 2)

MADA_SolDB_keep_meta <- NLR_filtered_deduplicated_len[NLR_filtered_deduplicated_len$ID %in% MADA_SolDB_keep$seqname,]

# keep only NLRs with CNL architecture
MADA_SolDB_keep_meta_filtered <- MADA_SolDB_keep_meta[MADA_SolDB_keep_meta$Simple == "CNL",]

MADA_SolDB_keep_filtered <- MADA_SolDB_keep[match(MADA_SolDB_keep_meta_filtered$ID,MADA_SolDB_keep$seqname),]

MADA_SolDB_keep_meta_filtered_CC <- CC[match(MADA_SolDB_keep_meta_filtered$ID,CC$ID),]

# make sure the MADA hit is within the CC domain boundaries
MADA_SolDB_keep_filtered <- MADA_SolDB_keep_filtered[MADA_SolDB_keep_filtered$to < MADA_SolDB_keep_meta_filtered_CC$end,]

MADA_SolDB_keep_meta_filtered <- MADA_SolDB_keep_meta_filtered[match(MADA_SolDB_keep_filtered$seqname, MADA_SolDB_keep_meta_filtered$ID),]
write_csv(MADA_SolDB_keep_meta_filtered, "/path/to/SolDB_MADA_ext_meta.csv")

plot(table(MADA_SolDB_keep_filtered$from), ylab = "Frequency", xlab = "MADA Motif Starting Position")


write_csv(MADA_SolDB_keep_filtered, "/path/to/SolDB_MADA_ext.csv")
write(MADA_SolDB_keep_filtered$seqname, "/path/to/SolDB_MADA_ext_IDs.txt")

MADA_SolDB_keep_filtered <- MADA_SolDB_keep_filtered %>% mutate(m_count_annot = case_when(
  m_count == 2 ~ "#ef3d3a",
  m_count > 2 ~ "#ffd111"))

MADA_SolDB_keep_filtered_annot <- MADA_SolDB_keep_filtered[,c(1,12)] %>% setNames(c("ID","color"))

## sequences with normal MADA hits
MADA_SolDB_normal <- MADA_SolDB_meta[!MADA_SolDB_meta$ID %in% MADA_SolDB_keep_filtered$seqname,]

MADA_SolDB_normal$annot <- "#3a3a39"
MADA_SolDB_normal_annot <- MADA_SolDB_normal[,c(17,18)] %>% setNames(c("ID","color"))

annot <- rbind(MADA_SolDB_normal_annot, MADA_SolDB_keep_filtered_annot)

write_csv(annot, "~/Desktop/external_projects/Rpi-vnt1_nterm/analyses/v2/MADA_ext/SolDB_MADA_annot.csv")


# sequences that passed the filter
MADA_SolDB_keep_seq <- NLR_filtered_deduplicated_len_seq_ref[MADA_SolDB_keep_filtered$seqname]
writeXStringSet(MADA_SolDB_keep_seq, "/path/to/SolDB_MADA_ext.fasta")

MADA_SolDB_keep_NBARC_seq <- NBARC_filtered_deduplicated_len_seq_ref[NBARC_filtered_deduplicated_len_seq_ref@ranges@NAMES %in% MADA_SolDB_keep_filtered$seqname]
writeXStringSet(MADA_SolDB_keep_NBARC_seq, "/path/to/SolDB_MADA_ext_NBARC.fasta")

```


## Import the clustering data from CD-HIT and extract the redundant hits:
```{r}
parse_cdhit_clusters <- function(clstr_file) {
  
  # Read all lines from the file
  lines <- readLines(clstr_file)
  
  # Initialize variables
  current_cluster <- NA
  representative_seq <- NA
  results <- list()
  result_index <- 1
  
  # Process each line
  for (line in lines) {
    line <- trimws(line)
    
    # Skip empty lines
    if (line == "") next
    
    # Check if this is a cluster header
    if (startsWith(line, ">Cluster")) {
      # Extract cluster number
      current_cluster <- as.integer(str_extract(line, "\\d+"))
      representative_seq <- NA
      next
    }
    
    # Parse sequence line
    # Extract sequence ID (between > and ...)
    seq_match <- str_match(line, ">([^.\\s]+(?:\\.[^.\\s]+)*)")
    
    if (!is.na(seq_match[1])) {
      seq_id <- seq_match[2]
      
      # Check if this is a representative sequence (ends with *)
      if (str_detect(line, "\\*$")) {
        # This is the representative sequence
        representative_seq <- seq_id
      } 
      # Check if this is a removed sequence (contains "at X%")
      else if (str_detect(line, "at\\s+[\\d.]+%")) {
        # Extract similarity percentage
        sim_match <- str_match(line, "at\\s+([\\d.]+)%")
        similarity <- as.numeric(sim_match[2])
        
        # Add to results
        if (!is.na(representative_seq)) {
          results[[result_index]] <- data.frame(
            representative_seq = representative_seq,
            removed_seq = seq_id,
            cluster_number = current_cluster,
            similarity_percentage = similarity,
            stringsAsFactors = FALSE
          )
          result_index <- result_index + 1
        }
      }
    }
  }
  
  # Combine all results into a single data frame
  if (length(results) > 0) {
    result_df <- bind_rows(results)
  } else {
    # Return empty data frame with correct columns
    result_df <- data.frame(
      representative_seq = character(),
      removed_seq = character(),
      cluster_number = integer(),
      similarity_percentage = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  return(result_df)
}

#' Remove suffix from gene IDs
#'
#' @param gene_ids Vector of gene IDs
#' @return Vector of gene IDs with suffixes removed
strip_suffix <- function(gene_ids) {
  # Remove suffix like _CNL, _TNL, _NL, etc.
  str_replace(gene_ids, "_[A-Z]+$", "")
}

#' Main function to run the script
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript parse_cdhit.R <input.clstr> <output.csv> [--strip-suffix]\n")
    cat("\nOptions:\n")
    cat("  --strip-suffix    Remove domain suffixes from gene IDs\n")
    cat("\nExample:\n")
    cat("  Rscript parse_cdhit.R clusters.clstr removed_genes.csv\n")
    cat("  Rscript parse_cdhit.R clusters.clstr removed_genes.csv --strip-suffix\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  output_file <- args[2]
  strip_suffixes <- "--strip-suffix" %in% args
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    stop(paste("Error: Input file", input_file, "does not exist"))
  }
  
  cat("Reading cluster file:", input_file, "\n")
  
  # Parse the cluster file
  df <- parse_cdhit_clusters(input_file)
  
  # Strip suffixes if requested
  if (strip_suffixes) {
    cat("Stripping domain suffixes from gene IDs...\n")
    df$representative_seq <- strip_suffix(df$representative_seq)
    df$removed_seq <- strip_suffix(df$removed_seq)
  }
  
  # Print summary
  cat("\n=== Summary ===\n")
  cat("Total removed sequences:", nrow(df), "\n")
  if (nrow(df) > 0) {
    cat("Unique clusters:", length(unique(df$cluster_number)), "\n")
    cat("Unique representative sequences:", length(unique(df$representative_seq)), "\n")
    cat("Similarity range:", min(df$similarity_percentage), "-", 
        max(df$similarity_percentage), "%\n")
  }
  
  # Save to CSV
  cat("\nSaving results to:", output_file, "\n")
  write.csv(df, output_file, row.names = FALSE)
  
  cat("Done!\n")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
```


```{r}

NLR_filtered_CDHIT <- parse_cdhit_clusters("/path/to/NLR_filtered_seq_clu.fasta.clstr")

MADA_SolDB_keep_filtered_redundant <- NLR_filtered_CDHIT[NLR_filtered_CDHIT$representative_seq %in% MADA_SolDB_keep_meta_filtered$ID,]

MADA_SolDB_keep_meta_filtered_redundant <- NLR_filtered[NLR_filtered$ID %in% MADA_SolDB_keep_filtered_redundant$removed_seq,]

MADA_SolDB_keep_filtered_all_meta <- rbind(MADA_SolDB_keep_meta_filtered, MADA_SolDB_keep_meta_filtered_redundant)
write_csv(MADA_SolDB_keep_filtered_all_meta, "/path/to/MADA_ext/SolDB_MADA_ext_meta_redundant.csv")

```

