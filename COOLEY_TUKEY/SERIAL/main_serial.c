#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>
#include <time.h>  // For clock_gettime

#include "./SERIAL_PROGRAM_LIBS/images_handling.h"
#include "./SERIAL_PROGRAM_LIBS/matrix_utilities.h"
#include "./SERIAL_PROGRAM_LIBS/FFTs.h"

int main(int argc, char *argv[]) {

    // Check if the correct number of arguments is provided
    if (argc != 4) {
        fprintf(stderr, "Correct usage: %s <image_path> <module_output_path> <phase_output_path>\n", argv[0]);
        printf("Arguments detected %d\n", argc);
        return 1;
    }
    
    printf("Process start...\n");

    const char *filename = argv[1];  // Input image file path
    Image *img = read_ppm(filename); // Read the input image in PPM format
    // Placeholder for the module and phase output
    Image *module;  
    Image *phase;   

    // Allocate memory for the module and phase images
    module = (Image *)malloc(sizeof(Image));
    phase = (Image *)malloc(sizeof(Image));

    // Pad the input image to prepare for FFT
    Image *padded_image = pad_image(img);

    // Perform FFT on the padded image and store results in module and phase images
    FFT_image(padded_image, module, phase);

    // Write the phase image to the specified output path with scaling
    writePPM(argv[3], phase, find_scale_factor(phase));

    // Apply logarithmic scaling to the module image for visualization and save it
    Image* module_log = log_scale(module);
    writePPM(argv[2], module_log, find_scale_factor(module_log));

    // Free dynamically allocated memory for images and complex matrices
    free_image(img);
    free_image(padded_image);
    free_image(module);
    free_image(phase);

    printf("...Process done.\n");

    return 0;
}
