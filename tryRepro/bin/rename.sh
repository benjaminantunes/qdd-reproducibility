#!/bin/bash

# Function to generate a unique temporary name
generate_temp_name() {
    local baseName=$1
    local count=0
    local tempName="${baseName}-temp-${count}"
    while [ -d "$tempName" ]; do
        ((count++))
        tempName="${baseName}-temp-${count}"
    done
    echo "$tempName"
}

# Iterate over directories in the current directory
for dir in */ ; do
    dir=${dir%/} # Remove the trailing slash for easier handling
    if [[ "$dir" = *-Debug ]]; then
        # Generate a temporary unique name
        tempName=$(generate_temp_name "$dir")
        # Rename to the temporary unique name
        mv "$dir" "$tempName"
    elif [[ "$dir" = *-noDebug ]]; then
        # Generate a temporary unique name
        tempName=$(generate_temp_name "$dir")
        # Rename to the temporary unique name
        mv "$dir" "$tempName"
    fi
done

# Rename from the temporary unique name to the final name
for dir in *-temp-*; do
    dir=${dir%/} # Remove the trailing slash for easier handling
    if [[ "$dir" = *-Debug-temp-* ]]; then
        # Replace the temporary name part with -noDebug
        finalName="${dir%-Debug-temp-*}-noDebug"
        mv "$dir" "$finalName"
    elif [[ "$dir" = *-noDebug-temp-* ]]; then
        # Replace the temporary name part with -Debug
        finalName="${dir%-noDebug-temp-*}-Debug"
        mv "$dir" "$finalName"
    fi
done

