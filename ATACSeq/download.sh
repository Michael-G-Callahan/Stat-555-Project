#!/bin/bash

# Path to the file with the bigBed file IDs
ID_FILE="bigBedIDs.txt"

# Base URL for ENCODE downloads
BASE_URL="https://www.encodeproject.org/files"

# Read file line by line
while IFS= read -r ID; do
    # Construct the download URL
    DOWNLOAD_URL="${BASE_URL}/${ID}/@@download/${ID}.bigBed"
    
    # Download the file using curl
    echo "Downloading ${ID}.bigBed..."
    wget "$DOWNLOAD_URL"
done < "$ID_FILE"

echo "Download complete."
