#!/bin/bash

ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/QDDpackage_withexamples/bin"
PATTERN="qdd-*"
SOURCE_FILE="/home/benjamin/experience-QDD-ServeurAMD/QDD-MarcExperiment/jpwzm9knmb-1/QDDpackage_withexamples/QDDpackage_withexamples/bin/qdd-gfortran-linux-FFTW-noOMP-noStatic-noDebug/repli2/"

# Create an array of directories to process
declare -a dirs=()

# Loop to append directories matching the pattern and having repli1 to repli10
for folder in $ROOT_DIR/$PATTERN; do
    for i in {11..99}; do
        if [ -d "$folder/repli$i" ]; then
            dirs+=("$folder/repli$i/")
        fi
    done
done

# Loop through the determined directories
for dir in "${dirs[@]}"; do
    echo "Processing directory: $dir"
    
    # Copy the files beginning with "for005" from the source directory to the current directory in the loop
    cp "${SOURCE_FILE}for005"* "$dir"
done

