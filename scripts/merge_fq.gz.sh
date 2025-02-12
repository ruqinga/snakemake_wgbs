#!/bin/bash

# Loop through all fastq.gz files
find ./ -name "*rep[0-9]*.fastq.gz" | while read file; do
    # Define the merged file by removing the _rep[0-9]* part of the filename
    merged_file="$(echo "$file" | sed 's/_rep[0-9]*//' | sed 's/.gz//')"

    # Check if the merged file already exists
    if [ ! -f "$merged_file" ]; then
        # If it doesn't exist, start by extracting the first file
        zcat "$file" > "$merged_file"
        echo "Merging $file into $merged_file"
    else
        # If it exists, append the contents of the current file
        zcat "$file" >> "$merged_file"
        echo "Merging $file into $merged_file"
    fi
done

gzip *.fastq

# After all files are merged
echo "Merging completed!"
