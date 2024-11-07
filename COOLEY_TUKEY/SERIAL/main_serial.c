#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "./SERIAL_PROGRAM_LIBS/images_handling.h"
#include "./SERIAL_PROGRAM_LIBS/matrix_utilities.h"
#include "./SERIAL_PROGRAM_LIBS/FFTs.h"

int main(int argc, char *argv[]) {

    if (argc != 4) {
        fprintf(stderr, "Correct usage: %s <image_path> <module_output_path> <phase_output_path>\n", argv[0]);
        printf("Arguments detected %d\n", argc);
        return 1;
    }
    
    printf("Process start...\n");

    const char *filename = argv[1];
    Image *img = read_ppm(filename);
    Image *module;
    Image *phase;

    module = (Image *)malloc(sizeof(Image));
    phase = (Image *)malloc(sizeof(Image));

    Image *padded_image = pad_image(img);

    int rows = 4, cols = 4;

    // Dynamically allocate memory for the matrix
    double complex **matrix = malloc(rows * sizeof(double complex*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(double complex));
    }

    // Initialize the matrix with some values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 1 + 0 * I;
        }
    }

    //double complex **result_matrix = matrix_FFT(matrix, 4,4);
    
    //print_complex_matrix(result_matrix, 4, 4);

    FFT_image(padded_image,module,phase);

    writePPM(argv[3], phase, find_scale_factor(phase));

    Image* module_log = log_scale(module);
    writePPM(argv[2], module_log, find_scale_factor(module_log));

    //printf("Image trnasformed and saved!\n");

    //printf("Freeing original image\n");
    free_image(img);
    //printf("Freeing padded image\n");
    free_image(padded_image);

    //printf("Freeing module image\n");
    free_image(module);
    //printf("Freeing phase image\n");
    free_image(phase);
    printf("...Process done.\n");
    
    return 0;
}
