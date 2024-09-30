#!/bin/bash

ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/tryRepro/bin"
PATTERN="qdd-*"
SOURCE_FILE="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/tryRepro/bin/qdd-gfortran-linux-FFTW-noOMP-noStatic-noDebug/repli2/"

# Create an array of directories to process
declare -a dirs=()
echo "Step 1: Identifying directories"

# Loop to append directories matching the pattern and having repli11 to repli50
for folder in "$ROOT_DIR"/$PATTERN; do
    echo "Checking folder: $folder"
    for i in {11..50}; do
        target_dir="$folder/repli$i"
        if [ ! -d "$target_dir" ]; then
            echo "Directory does not exist: $target_dir. Creating directory."
            mkdir -p "$target_dir"
        fi
        echo "Found or created directory: $target_dir"
        dirs+=("$target_dir/")
    done
done

echo "Step 2: Copying files"
# Loop through the determined directories
for dir in "${dirs[@]}"; do
    echo "Processing directory: $dir"
    cp "${SOURCE_FILE}for005"* "$dir"
    if [ $? -eq 0 ]; then
        echo "Files copied successfully to $dir"
    else
        echo "Failed to copy files to $dir"
    fi
done

echo "Script completed"

