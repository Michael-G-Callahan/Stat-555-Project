#!/bin/bash

# Directory containing your BigBed files
# Replace "/path/to/your/directory" with the actual path to your BigBed files
DIRECTORY="."

# Loop through each BigBed file in the directory
for input_file in "$DIRECTORY"/*; do
    # Check if the file is a BigBed file
    if [[ $input_file == *.bigBed ]]; then
        # Construct the output .bed file name
        output_file="${input_file}.bed"
        # Execute the conversion command
        echo "Converting $input_file to $output_file"
        ./bigBedToBed "$input_file" "$output_file"
    fi
done

echo "Conversion complete."
