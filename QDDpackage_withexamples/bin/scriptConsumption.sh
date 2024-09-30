#!/bin/bash

ROOT_DIR="/home/benjamin/experience-QDD-ServeurIntel/QDD-MarcExperimentVERIF/jpwzm9knmb-1/QDDpackage_withexamples/QDDpackage_withexamples/bin"
PATTERN="qdd-*"

# Set the environment variable for OpenMP
export OMP_NUM_THREADS=16
source /opt/intel/oneapi/setvars.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/mkl/latest/lib/intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin
# Create an array of directories to process
declare -a dirs=()

# Loop to append directories matching the pattern and having repli1 to repli10
for folder in $ROOT_DIR/$PATTERN; do
    for i in {31..50}; do
        if [ -d "$folder/repli$i" ]; then
            dirs+=("$folder/repli$i/")
        fi
    done
done

# Loop through the determined directories
for dir in "${dirs[@]}"; do
    echo "Processing directory: $dir"
    
    # Create a temporary directory
    TMP_DIR=$(mktemp -d)
    echo "Created temporary directory: $TMP_DIR"

    # Copy the contents of the target directory to the temporary directory
    echo "Copying contents of $dir to $TMP_DIR"
    cp -r "${dir}." "$TMP_DIR/"

    # Change to the temporary directory
    cd "$TMP_DIR"

    # Run ./qdd in the background for the PowerJoular measurement
    echo "Running ./qdd in the background"
    ./qdd > /dev/null 2>&1 &
    PROGRAM_PID=$!

    # Define where the PowerJoular output should be saved
    POWER_FILE="consumption.txt"

    # Run PowerJoular in the background
    echo "Running powerjoular in the background"
    powerjoular -p $PROGRAM_PID >> $POWER_FILE &
    POWERJOULAR_PID=$!

    # Sleep for 60 seconds
    echo "Sleeping for 60 seconds"
    sleep 60

    # Kill PowerJoular and the ./qdd process
    echo "Killing powerjoular and ./qdd processes"
    kill -INT $POWERJOULAR_PID
    kill -INT $PROGRAM_PID

    # Copy consumption.txt back to the original directory
    echo "Copying $POWER_FILE back to $dir"
    cp "$POWER_FILE" "$dir"

    # Return to the parent directory
    cd -

    # Delete the temporary directory
    echo "Deleting temporary directory: $TMP_DIR"
    rm -rf "$TMP_DIR"
done

