#!/bin/bash

# Define arrays of options
compilers=("gfortran" "ifort" "ifx")
debug_options=("YES" "NO")
omp_options=("YES" "NO")
dynomp_options=("YES" "NO") # Only relevant if OMP is YES
fft_types=("FFTW" "MKL")
link_static="NO" # Assuming you always want NO for static linking

# Base directory for the executable and folders
base_dir="../../bin"

# Loop through each combination of options
for compiler in "${compilers[@]}"; do
    for debug in "${debug_options[@]}"; do
        for omp in "${omp_options[@]}"; do
            for fft_type in "${fft_types[@]}"; do
                # DYNOMP only makes sense if OMP is YES
                if [ "$omp" == "YES" ]; then
                    dynomp_options=("YES" "NO")
                else
                    dynomp_options=("NO")
                fi
                for dynomp in "${dynomp_options[@]}"; do
                    # Modify Makefile with current options
                    sed -i "s/^COMPILER = .*/COMPILER = $compiler/" Makefile
                    sed -i "s/^DEBUG = .*/DEBUG = $debug/" Makefile
                    sed -i "s/^OMP = .*/OMP = $omp/" Makefile
                    sed -i "s/^DYNOMP = .*/DYNOMP = $dynomp/" Makefile
                    sed -i "s/^FFT_TYPE = .*/FFT_TYPE = $fft_type/" Makefile
                    sed -i "s/^LINK_STATIC = .*/LINK_STATIC = $link_static/" Makefile

                    # Clean and make
                    make clean
                    make

                    # Construct folder name based on options
                    folder_name="qdd-${compiler}-linux-${fft_type}-"
                    if [ "$omp" == "YES" ]; then
                        folder_name+="OMP-"
                    else
                        folder_name+="noOMP-"
                    fi
                    if [ "$dynomp" == "YES" ]; then
                        folder_name+="DYN-"
                    else
                        folder_name+=""
                    fi
                    
                    folder_name+="noStatic-"
                    
                    if [ "$debug" == "NO" ]; then
                        folder_name+="noDebug"
                    else
                        folder_name+="Debug"
                    fi

                    # Ensure the target base directory exists and navigate to it
                    mkdir -p "$base_dir"
                    cd "$base_dir" || exit

                    # Create the main folder if it doesn't exist
                    mkdir -p "$folder_name"

                    # Loop through each repli folder
                    for i in {11..99}; do
                        repli_folder="${folder_name}/repli$i"
                        mkdir -p "$repli_folder" # Ensure the repli folder exists
                        # Copy the executable to the repli folder
                        cp "qdd" "$repli_folder/"
                    done

                    # Navigate back to the script directory
                    cd - || exit
                done
            done
        done
    done
done

