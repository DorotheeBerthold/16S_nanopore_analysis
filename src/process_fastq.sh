# Directory for concatenated files - where the fastq files should go
output_dir="concatenated_files_run"

# Create the concatenated_files_run1 directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each barcode folder, properly formatted with leading zeros
for i in $(seq -w 1 22); do
    barcode="barcode${i}"

    # Check if the barcode directory exists
    if [ -d "$barcode" ]; then
        echo "Processing $barcode"

        # Change to the barcode directory
        cd "$barcode"

        # Check if there are any .fastq.gz files to process
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            # Gunzip all .fastq.gz files in the current barcode folder
            gunzip *.fastq.gz
        else
            echo "No .fastq.gz files found in $barcode"
        fi

        # Check if there are any .fastq files to concatenate
        if ls *.fastq 1> /dev/null 2>&1; then
            # Concatenate all .fastq files in the current barcode folder into one .fastq file named after the folder
            cat *.fastq > ../"$output_dir"/${barcode}.fastq
        else
            echo "No .fastq files found in $barcode after decompression"
        fi

        # Change back to the parent directory
        cd ..
    else
        echo "Directory $barcode does not exist."
    fi
done

