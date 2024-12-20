#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "matrix_utilities.h"

/*
    This function allocates a float matrix and returns the pointer
    to the matrix as an output.

    It takes the width(NUMBER OF COLUMNS) and the height(NUMBER OF ROWS)
    as an input. a call with 5,20 produces a 20x5 matrix
*/

double **allocate_matrix(int width, int height) {
    double **matrix = (double **)malloc(height * sizeof(double *));
    for (int i = 0; i < height; i++) {
        matrix[i] = (double *)malloc(width * sizeof(double));
    }
    return matrix;
}

/*
    This function allocates a complex matrix, returns a pointer as 
    the one above and also uses the same convention to talk about
    rows and columns
*/

double complex **allocate_complex_matrix(int width, int height) {
    double complex **matrix = (double complex **)malloc(height * sizeof(double complex *));
    for (int i = 0; i < height; i++) {
        matrix[i] = (double complex *)malloc(width * sizeof(double complex));
    }
    return matrix;
}

/*
    This function frees a double matrix. It sets
    the pointer to NULL for safety reasons in order
    to perform a check on the pointer for correct freeing

    the height information is needed to free
*/

void free_matrix(double **matrix, int height) {
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
        matrix[i] = NULL;  // Set to NULL after freeing
    }
    free(matrix);
    matrix = NULL;
}

/*
    Same as above but for complex matrix
*/

void free_complex_matrix(double complex **matrix, int height) {
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
        matrix[i] = NULL;
    }
    free(matrix);
    matrix = NULL;
}



double complex **convert_to_complex_matrix(double **matrix, int rows, int cols) {
    // Allocate memory for the complex matrix
    double complex **complexMatrix = (double complex **)malloc(rows * sizeof(double complex *));
    if (complexMatrix == NULL) {
        perror("Failed to allocate memory for complex matrix");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        complexMatrix[i] = (double complex *)malloc(cols * sizeof(double complex));
        if (complexMatrix[i] == NULL) {
            perror("Failed to allocate memory for complex matrix row");
            exit(EXIT_FAILURE);
        }
    }

    // Convert double matrix to complex matrix (real part is the original value, imaginary part is 0)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            complexMatrix[i][j] = matrix[i][j] + 0 * I; // Real part is from the input matrix, Imaginary part is 0
        }
    }

    return complexMatrix;
}

/*
    Function used to print a complex matrix

    It takes as an input the matrix pointer and the number
    of rows and columns.
*/

void print_complex_matrix(double complex **matrix, int rows, int cols) {
    printf("Complex Matrix:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("(%6.2f + %6.2fi)  ", creal(matrix[i][j]), cimag(matrix[i][j]));
        }
        printf("\n");
    }
}

/*
    Function used to print a vector of integers

    it takes as an input the vector pointer and the number
    of elements contained in it.

    Might be implement a version that uses size() to calculate
    the number of elements
*/

void print_int_vector(int *vect, int lenght){
    printf("[");
    int i;
    for(i = 0; i < lenght-1; i++){
        printf("%d ",vect[i]);
    }
    printf("%d]\n",vect[i]);
}
