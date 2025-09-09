#!/bin/bash

# Directory containing the FASTQ files
input_dir="./emu_database/data/DB087/filtered_reads"

# Parameters for the emu abundance command
threads=106
db="emu_database"
min_abundance=0.00001

for i in $(seq 1 96); do
    if [ $i -lt 10 ]; then
        barcode="barcode0$i"
    else
        barcode="barcode$i"
    fi

    file="$input_dir/${barcode}_filtered.fastq"

    if [[ -f "$file" ]]; then
        emu abundance "$file" --keep-counts --threads "$threads" --db "$db" --min-abundance "$min_abundance"
    fi
done
