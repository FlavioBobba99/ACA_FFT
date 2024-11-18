#!/bin/bash

# Output name for the compiled program
OUTPUT_NAME="../parallel"

# Source files to compile
SRC_FILES="../main_parallel.c ../PARALLEL_PROGRAM_LIBS/matrix_utilities.c ../PARALLEL_PROGRAM_LIBS/FFTs.c ../PARALLEL_PROGRAM_LIBS/images_handling.c"

# GCC compiler flags
CFLAGS="-Wall"

LIBFLAG="-lm"

# Compile the program
echo "Compiling the program..."
mpicc $CFLAGS $SRC_FILES -o $OUTPUT_NAME $LIBFLAG

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful."
    echo "Output: $OUTPUT_NAME"
else
    echo "Compilation failed."
fi
