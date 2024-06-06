#!/bin/bash

# Directory containing the FASTQ files - change to where the filtered reads are
input_dir="./data/DB044/240529_run1/filtered_reads"

# Parameters for the emu abundance command
threads=10
db="emu_database"
min_abundance=0.00001

# Loop over each FASTQ file in the directory
for file in "$input_dir"/*_filtered.fastq; do
    # Run the emu abundance command on the current file
    emu abundance "$file" --keep-counts --threads $threads --db $db --min-abundance $min_abundance
done

