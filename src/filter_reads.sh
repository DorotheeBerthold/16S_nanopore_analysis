# Directory containing the concatenated fastq files
input_dir="concatenated_files_run"

# Directory for filtered reads
output_dir="filtered_reads"

# Create the filtered_reads directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each .fastq file in the input directory
for fastq_file in "$input_dir"/*.fastq; do
    # Extract the base name of the file (e.g., barcode01 from barcode01.fastq)
    base_name=$(basename "$fastq_file" .fastq)

    # Define the output file name
    output_file="$output_dir/${base_name}_filtered.fastq"

    # Run cutadapt with the specified parameters
    cutadapt -m 1300 -M 1700 -o "$output_file" "$fastq_file"

    # Print a message indicating the completion of the file processing
    echo "Filtered reads saved to $output_file"
done

