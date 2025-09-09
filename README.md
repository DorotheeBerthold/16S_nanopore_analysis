# 16S_nanopore_analysis

## Overview

**16S_nanopore_analysis** provides a streamlined and reproducible workflow for analyzing 16S rRNA gene sequencing data generated with Oxford Nanopore’s MinION. It processes raw reads and annotates them using EMU for taxonomic profiling. This toolkit is ideal for researchers working on microbial community analysis, metagenomics, and microbiome studies using portable long-read sequencing.

## Features

- **End-to-end workflow**: From raw FASTQ reads to annotated taxonomic profiles using EMU.
- **Multi-language integration**: Combines Shell scripts for data handling with R for downstream analysis and visualization.
- **Reproducible and modular**: Structured code organization and clear documentation support repeatable analysis.

## Repository Structure

## Installation & Setup

### Dependencies

This workflow relies on the following tools and packages:

- **Shell environment** (Linux/macOS) with:
  - `bash`
  - `awk`, `sed`, `grep` (standard GNU command-line tools)
  - `cutadapt` -- create own venv with cutadapt
  -  `emu` -- create own venv with emu [githublink] (https://github.com/treangenlab/emu)
- **R** (version ≥ 4.0.0) with the following packages:
  - `tidyverse`
  - `data.table`
  - `ggplot2` *(and any others used in your R scripts)*

### Setting Up

1. Clone the repo:

   ```bash
   git clone https://github.com/DorotheeBerthold/16S_nanopore_analysis.git
   cd 16S_nanopore_analysis
