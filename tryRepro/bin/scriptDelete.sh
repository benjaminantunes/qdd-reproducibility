#!/bin/bash

ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/QDDpackage_withexamples/bin"
PATTERN="qdd-*"


# Create an array of directories to process
declare -a dirs=()

# Loop to append directories matching the pattern and having repli1 to repli10
for folder in $ROOT_DIR/$PATTERN; do
    for i in {1..10}; do
        if [ -d "$folder/repli$i" ]; then
            dirs+=("$folder/repli$i/")
        fi
    done
done

# Loop through the determined directories
for dir in "${dirs[@]}"; do
    echo "Processing directory: $dir"
    
    # Remove the files from the specific directory
    rm "${dir}"* 2>/dev/null
    #rm "${dir}"*.txt 2>/dev/null
    #rm "${dir}"*.dat 2>/dev/null
    #rm "${dir}"qdd 2>/dev/null
    

done

