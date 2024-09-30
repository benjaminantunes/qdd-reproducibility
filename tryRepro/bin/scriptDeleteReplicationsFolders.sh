#!/bin/bash

# Define the root directory and pattern to find subdirectories
ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/tryRepro/bin"
PATTERN="qdd-*"

# Loop through all subdirectories matching the pattern
for folder in $ROOT_DIR/$PATTERN; do
    echo "Processing subdirectory: $folder"

    # Loop through repli folders from repli51 to repli99
    for i in {51..99}; do
        REPLI_DIR="$folder/repli$i"
        if [ -d "$REPLI_DIR" ]; then
            echo "Deleting folder: $REPLI_DIR"
            rm -rf "$REPLI_DIR"
        else
            echo "Folder $REPLI_DIR does not exist."
        fi
    done
done

echo "Deletion process completed."

