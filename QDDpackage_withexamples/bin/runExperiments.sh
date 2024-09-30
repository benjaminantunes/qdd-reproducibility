#!/bin/bash
source /opt/intel/oneapi/setvars.sh


# Set the environment variable for OpenMP
export OMP_NUM_THREADS=16


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/mkl/latest/lib/intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/lib/intel64_lin

# Save the initial directory
initial_dir=$(pwd)

# Loop through all directories in the current folder
for dir in */ ; do
    echo "Entering main directory: $dir"
    
    # Change to the directory
    pushd "$dir" || { echo "Error changing directory to $dir"; continue; }

    # Loop through the subdirectories "repli1" and "repli2"
    for i in {31..50}; do
        sub_dir="repli$i"  # Construct the directory name dynamically
        echo "Entering subdirectory: $sub_dir"

        # Change to the subdirectory
        pushd "$sub_dir" || { echo "Error changing directory to $sub_dir"; continue; }


	

        # Run the command and time it
        { time ./qdd > /dev/null 2>&1; } 2> time.txt


        # Pop back to the main directory
        popd || { echo "Error changing back to the main directory"; exit 1; }

        echo "Completed work in subdirectory: $sub_dir"
    done

    # Pop back to the initial directory
    popd || { echo "Error changing back to the original directory"; exit 1; }

    echo "Completed work in main directory: $dir"
done

