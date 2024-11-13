#!/bin/bash

# Directory containing the FASTQ files
input_dir="./emu_database/data/DB059/filtered_reads"

# Directory to store the results
results_dir="./emu_database/data/DB059/results"
mkdir -p "$results_dir"  # Create the results directory if it doesn't exist

# Parameters for the emu abundance command
threads=106
db="emu_database"
min_abundance=0.00001

# Loop over each FASTQ file in the directory
for file in "$input_dir"/*_filtered.fastq; do
    # Extract the filename without the path
    filename=$(basename "$file")
    
    # Define the output path for the results
    output_file="$results_dir/${filename%.fastq}_abundance_results.txt"
    
    # Run the emu abundance command on the current file, directing output to the results directory
    emu abundance "$file" --keep-counts --threads $threads --db $db --min-abundance $min_abundance > "$output_file"
done

emu combine-outputs $results_dir species
