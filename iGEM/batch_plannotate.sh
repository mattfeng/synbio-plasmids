#!/bin/bash

# Check for correct number of arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <output_directory>"
    exit 1
fi

# Get the output directory from the command line argument
OUTPUT_DIR="$1"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if Plannotate is installed
if ! command -v plannotate &> /dev/null; then
    echo "Error: Plannotate is not installed or not in PATH."
    exit 1
fi

# Process each FASTA file in the input directory
for fasta_file in ./*.fasta; do
    if [[ -f "$fasta_file" ]]; then
        echo "Processing $fasta_file..."

        # Run Plannotate and save results in the output directory
        plannotate batch -i "$fasta_file" --output "$OUTPUT_DIR" > /dev/null

        if [[ $? -eq 0 ]]; then
            echo "Finished processing $fasta_file"
        else
            echo "Error processing $fasta_file"
        fi
    fi
done

echo "All files processed. Results are in $OUTPUT_DIR."