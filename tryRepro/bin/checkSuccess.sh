#!/bin/bash

ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/tryRepro/bin"
PATTERN="qdd-*"

# Create an array of directories to process
declare -a dirs=()

# Loop to append directories matching the pattern and having repli11 to repli50
for folder in $ROOT_DIR/$PATTERN; do
    for i in {11..50}; do
        if [ -d "$folder/repli$i" ]; then
            dirs+=("$folder/repli$i/")
        fi
    done
done

# Check for the presence of consumption.txt in each directory
for dir in "${dirs[@]}"; do
    if [ -f "$dir/consumption.txt" ]; then
        echo "Directory: $dir - Found: consumption.txt"
    else
        echo "Directory: $dir - Not found: consumption.txt"
    fi
done

