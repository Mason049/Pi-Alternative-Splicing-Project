# Pi-Alternative-Splicing-Project

A bioinformatics pipeline for identifying and analyzing Plant NLR genes. Includes scripts for genomic coordinate matching, sequence extraction, and downstream analysis of N-terminal MADA motifs and alternative splicing signatures.

---

# SolDB NLR & MADA Motif Analysis Pipeline

This repository contains a comprehensive bioinformatics workflow for identifying and analyzing Plant Nucleotide-binding Leucine-rich Repeat (NLR) genes. The pipeline spans from genomic coordinate extraction to deep motif analysis of the N-terminal MADA motif.



## Project Structure

### Python Scripts (Data Processing)
* **`extract_gene_coords.py`**: Extracts genomic coordinates by matching metadata entries to GFF files. It supports recursive directory searching and various matching modes like `strict-id` or `id+name`.
* **`extract_gene_sequences.py`**: Fetches genomic sequences from FASTA files based on coordinates. It is strand-aware and allows for the inclusion of upstream and downstream flanking regions.

### R Analysis (Downstream Pipeline)
* **`SolDB_NLR_v1.Rmd`**: Handles the initial integration of `NLRtracker` data with species metadata. It performs classification (CNL, TNL, etc.), deduplication using CD-HIT results, and NBARC domain extraction.
* **`Rpi-vnt1_nterm_v3.Rmd`**: A specialized script for analyzing the MADA motif. It filters sequences based on HMMsearch outputs, identifies specific N-terminal features (e.g., secondary 'M' residues), and manages redundant hits.

---

## Getting Started

### 1. Prerequisites
**Python Environment**:
* Python 3.x
* `pandas`

**R Environment**:
* `tidyverse`, `Biostrings`, `readxl`, `ggtree`, `svglite`

### 2. Recommended Usage Order
1. **Coordinate Mapping**: Use `extract_gene_coords.py` with your metadata and GFF root folder to generate a coordinate table.
2. **Sequence Extraction**: Use the coordinate table with `extract_gene_sequences.py` to pull FASTA sequences from your genome files.
3. **NLR Filtering**: Run `SolDB_NLR_v1.Rmd` to consolidate `NLRtracker` outputs and filter for specific NLR architectures.
4. **Motif Analysis**: Run `Rpi-vnt1_nterm_v3.Rmd` to perform high-resolution analysis of the MADA motif and N-terminal leaders.

---

## Datasets
The sequence coordinates and extracted gene sequences used in this pipeline are expected to be organized in a local `data/` folder. For large-scale genomic datasets (e.g., GFF and FASTA files), please refer to your respective institutional genomic databases or public repositories like **Sol Genomics Network**.

## Author
This pipeline is developed by **AmirAli Toghani**. For questions, comments, or technical support, please contact the author or submit an issue on this GitHub repository.

## Copyright
Copyright (C) 2025 **AmirAli Toghani**. Distributed under the MIT License.
