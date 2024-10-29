#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "FFTs.h"

// local modules

#include "matrix_utilities.h"
#include "images_handling.h"

#define PI 3.14159265358979323846

/*
    This function calclates the FFT of a given vector

    The input is a complex vector of lenght n and the ouptput is the same
    This method implements the cooley tukey scheme for calculating the FFT
    it exploits recursive calls to iteslf in order to best handle the butterfly
    splits generated by the cooley tukey
*/

double complex* FFT_complex(double complex *provavett, int lengthvett) {
    // Base case: if there's only one element, return it
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

    // Allocate memory for even and odd arrays
    double complex *even = malloc(split_length * sizeof(double complex));
    double complex *odd = malloc(split_length * sizeof(double complex));
    if (even == NULL || odd == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Separate the input into even and odd indexed elements
    for (int i = 0; i < split_length; i++) {
        even[i] = provavett[2 * i];
        odd[i] = provavett[2 * i + 1];
    }

    // Recursively compute FFT for even and odd parts
    double complex *y_even = FFT_complex(even, split_length);
    double complex *y_odd = FFT_complex(odd, split_length);

    // Allocate memory for output
    double complex *out = malloc(lengthvett * sizeof(double complex));
    if (out == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Compute the FFT combining step
    for (int i = 0; i < split_length; i++) {
        double complex twiddle_factor = cexp(-2.0 * I * PI * i / lengthvett); // e^(-2*pi*i*k/N)
        out[i] = y_even[i] + twiddle_factor * y_odd[i];
        out[i + split_length] = y_even[i] - twiddle_factor * y_odd[i];
    }

    // Free allocated memory for even and odd arrays
    free(even);
    free(odd);
    free(y_even);
    free(y_odd);

    return out;
}

void FFT_complex_with_range(double complex *provavett, int starting_index, double complex *output_vector, int lengthvett) {
    // Base case: if there's only one element, return it
    if (lengthvett == 1) {
        printf("ERROR HOW DID YOU END UP HERE??? lenghtvect %d\n", lengthvett);
        double complex *single_out = malloc(sizeof(double complex));
        if (single_out == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        single_out[0] = provavett[0];
        //return single_out;
    }

    printf("lengthvett %d\n", lengthvett);

    int split_length = lengthvett / 2;

    // Allocate memory for even and odd arrays
    double complex *even = malloc(split_length * sizeof(double complex));
    double complex *odd = malloc(split_length * sizeof(double complex));
    if (even == NULL || odd == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    printf("Population of even and odd reached\n");

    // Separate the input into even and odd indexed elements
    for (int i = 0; i < split_length; i++) {
        printf("EVEN INDEX: %d\n", 2 * i + starting_index);
        printf("ODD INDEX: %d\n",2 * i + 1 + starting_index);
        even[i] = provavett[2 * i + starting_index];
        odd[i] = provavett[2 * i + 1 + starting_index];
    }

    print_complex_vector(even, split_length);
    print_complex_vector(odd, split_length);

    // Recursively compute FFT for even and odd parts
    double complex *y_even = FFT_complex(even, split_length);
    double complex *y_odd = FFT_complex(odd, split_length);

    print_complex_vector(y_even, split_length);
    print_complex_vector(y_odd, split_length);

    printf("FFT even and odd done\n");

    // Compute the FFT combining step
    for (int i = 0; i < split_length; i++) {
        double complex twiddle_factor = cexp(-2.0 * I * PI * i / lengthvett); // e^(-2*pi*i*k/N)
        output_vector[i + starting_index] = y_even[i] + twiddle_factor * y_odd[i];
        output_vector[i + split_length + starting_index] = y_even[i] - twiddle_factor * y_odd[i];
    }

    printf("Combination done\n");

    // Free allocated memory for even and odd arrays
    free(even);
    free(odd);
    free(y_even);
    free(y_odd);
}

void transpose(complex double **matrix, complex double **result, int widht, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            result[j][i] = matrix[i][j];  // Transpose operation
        }
    }
}

/*
    This functions peforms the FFT of a matrix. Input is a matrix the output is a complex
    matrix heigth and widdth are required to make this work

    it calculatest the fft first row by row then transposes the result and then FFT and transposition
    again in orer to have the correct output
*/

double complex **matrix_FFT (double complex **matrix, int width, int height) { 

    double complex **temporary_transposed_matrix = allocate_complex_matrix(height,width);
    double complex **out_matrix = allocate_complex_matrix(width,height);
    
	for (int i = 0; i < height; i++) {
            double complex *out_vect = FFT_complex(matrix[i], width);
            for(int j = 0; j< width; j++){
                temporary_transposed_matrix[j][i] = out_vect[j];
            }
            free(out_vect);
		}
    printf("----------------- DEBUG ----------------\n");
    
    for (int i = 0; i < width; i++) {
            double complex *out_vect = FFT_complex(temporary_transposed_matrix[i], height);
            for(int j = 0; j < height; j++){
                out_matrix[j][i] = out_vect[j];
            }
            free(out_vect);
		}

    free_complex_matrix(temporary_transposed_matrix, width);
	
	return out_matrix;
}

void FFT_image(Image *in, Image *module, Image *phase){

    int height = in->height;
    int width = in->width;

    module->height = height;
    module->width = width;
    module->max_color = 255;
    phase->height = height;
    phase->width = width;
    phase->max_color = 255;

    // Check if 'in' has valid pointers for 'red'
    if (in->red == NULL) {
        fprintf(stderr, "Input image red channel is NULL\n");
        return;
    }

    double complex **temp_red = convert_to_complex_matrix(in->red, height, width);
    double complex **temp_green = convert_to_complex_matrix(in->green, height, width);
    double complex **temp_blue = convert_to_complex_matrix(in->blue, height, width);

    printf("Matrices converted :)\n");

    // Check if matrix_FFT handles memory properly
    double complex **complex_red = matrix_FFT(temp_red, width, height);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT RED computation failed\n");
        return;
    }

    double complex **complex_green = matrix_FFT(temp_green, width, height);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT GREEN computation failed\n");
        return;
    }

    double complex **complex_blue = matrix_FFT(temp_blue, width, height);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT BLUE computation failed\n");
        return;
    }


    // Free allocated memory (example, modify as necessary)
    free_complex_matrix(temp_red, height);
    free_complex_matrix(temp_green, height);
    free_complex_matrix(temp_blue, height);

    module->red = allocate_matrix(width, height);
    phase->red = allocate_matrix(width, height);
    module->green = allocate_matrix(width, height);
    phase->green = allocate_matrix(width, height);
    module->blue = allocate_matrix(width, height);
    phase->blue = allocate_matrix(width, height);

    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            module->red[i][j] =  cabs(complex_red[i][j]); 
            phase->red[i][j] = carg(complex_red[i][j]);
            module->green[i][j] = cabs(complex_green[i][j]);
            phase->green[i][j] = carg(complex_green[i][j]);
            module->blue[i][j] = cabs(complex_blue[i][j]);
            phase->blue[i][j] =  carg(complex_blue[i][j]);
            
        }
    }

}
