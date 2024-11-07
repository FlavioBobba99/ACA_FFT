#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "FFTs.h"

// local modules
#include "matrix_utilities.h"
#include "images_handling.h"

#define PI 3.14159265358979323846  // Define constant value for PI

/*
    Calculates the FFT of a given complex vector using the Cooley-Tukey FFT algorithm.
    Input:
        - provavett: pointer to the input complex vector
        - lengthvett: length of the vector
    Output:
        - Returns a pointer to the FFT-transformed complex vector
*/
double complex* FFT_complex(double complex *provavett, int lengthvett) {
    // Base case for recursion: if only one element, return it directly
    if (lengthvett == 1) {
        double complex *single_out = malloc(sizeof(double complex));
        if (single_out == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        single_out[0] = provavett[0];
        return single_out;
    }

    int split_length = lengthvett / 2;

    // Allocate memory for even and odd indexed arrays
    double complex *even = malloc(split_length * sizeof(double complex));
    double complex *odd = malloc(split_length * sizeof(double complex));
    if (even == NULL || odd == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Split the input vector into even and odd elements
    for (int i = 0; i < split_length; i++) {
        even[i] = provavett[2 * i];
        odd[i] = provavett[2 * i + 1];
    }

    // Recursively calculate FFT for both halves
    double complex *y_even = FFT_complex(even, split_length);
    double complex *y_odd = FFT_complex(odd, split_length);

    // Allocate memory for the final output vector
    double complex *out = malloc(lengthvett * sizeof(double complex));
    if (out == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Combine FFT results from even and odd vectors using the twiddle factor
    for (int i = 0; i < split_length; i++) {
        double complex twiddle_factor = cexp(-2.0 * I * PI * i / lengthvett);
        out[i] = y_even[i] + twiddle_factor * y_odd[i];
        out[i + split_length] = y_even[i] - twiddle_factor * y_odd[i];
    }

    // Free memory for even, odd, and partial FFT results
    free(even);
    free(odd);
    free(y_even);
    free(y_odd);

    return out;
}

/*
    Transposes a given matrix.
    Input:
        - matrix: pointer to the input complex matrix
        - result: pointer to the output matrix (result of the transposition)
        - width: width of the input matrix
        - height: height of the input matrix
*/
void transpose(complex double **matrix, complex double **result, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            result[j][i] = matrix[i][j];  // Transpose operation
        }
    }
}

/*
    Calculates the FFT of a matrix.
    Input:
        - matrix: pointer to the input complex matrix
        - width: width of the matrix
        - height: height of the matrix
    Output:
        - Returns a pointer to the FFT-transformed matrix
*/
double complex **matrix_FFT(double complex **matrix, int width, int height) { 
    // Temporary matrices for transposing and output
    double complex **temporary_transposed_matrix = allocate_complex_matrix(height, width);
    double complex **out_matrix = allocate_complex_matrix(width, height);
    
    // Perform FFT row-wise on the matrix
    for (int i = 0; i < height; i++) {
        double complex *out_vect = FFT_complex(matrix[i], width);
        for (int j = 0; j < width; j++) {
            temporary_transposed_matrix[j][i] = out_vect[j];
        }
        free(out_vect);
    }
    
    // Perform FFT column-wise on the transposed result
    for (int i = 0; i < width; i++) {
        double complex *out_vect = FFT_complex(temporary_transposed_matrix[i], height);
        for (int j = 0; j < height; j++) {
            out_matrix[j][i] = out_vect[j];
        }
        free(out_vect);
    }

    // Free temporary matrix memory
    free_complex_matrix(temporary_transposed_matrix, width);
    
    return out_matrix;
}

/*
    Performs FFT on each color channel of an image and calculates module and phase.
    Input:
        - in: pointer to the input image
        - module: pointer to the image to store the module result
        - phase: pointer to the image to store the phase result
*/
void FFT_image(Image *in, Image *module, Image *phase) {
    int height = in->height;
    int width = in->width;

    // Initialize module and phase image dimensions
    module->height = height;
    module->width = width;
    module->max_color = 255;
    phase->height = height;
    phase->width = width;
    phase->max_color = 255;

    // Check if input image has valid red channel data
    if (in->red == NULL) {
        fprintf(stderr, "Input image red channel is NULL\n");
        return;
    }

    // Convert each color channel of the image to complex format for FFT
    double complex **temp_red = convert_to_complex_matrix(in->red, height, width);
    double complex **temp_green = convert_to_complex_matrix(in->green, height, width);
    double complex **temp_blue = convert_to_complex_matrix(in->blue, height, width);

    // Perform FFT on each color channel
    double complex **complex_red = matrix_FFT(temp_red, width, height);
    double complex **complex_green = matrix_FFT(temp_green, width, height);
    double complex **complex_blue = matrix_FFT(temp_blue, width, height);

    // Free the temporary complex matrix memory for each color channel
    free_complex_matrix(temp_red, height);
    free_complex_matrix(temp_green, height);
    free_complex_matrix(temp_blue, height);

    // Allocate memory for storing module and phase values in each color channel
    module->red = allocate_matrix(width, height);
    phase->red = allocate_matrix(width, height);
    module->green = allocate_matrix(width, height);
    phase->green = allocate_matrix(width, height);
    module->blue = allocate_matrix(width, height);
    phase->blue = allocate_matrix(width, height);

    // Calculate the module and phase for each color channel and store results
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            module->red[i][j] = cabs(complex_red[i][j]); 
            phase->red[i][j] = carg(complex_red[i][j]);
            module->green[i][j] = cabs(complex_green[i][j]);
            phase->green[i][j] = carg(complex_green[i][j]);
            module->blue[i][j] = cabs(complex_blue[i][j]);
            phase->blue[i][j] = carg(complex_blue[i][j]);
        }
    }
}


// I DON?T KNOW IT IS BORING RN
