#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

extern double **allocate_matrix(int width, int height);

extern double complex **allocate_complex_matrix(int width, int height);

extern void free_matrix(double **matrix, int height);

extern void free_complex_matrix(double complex **matrix, int height);

extern double complex **convert_to_complex_matrix(double **matrix, int rows, int cols);

extern void print_complex_matrix(double complex **matrix, int rows, int cols);

extern void print_int_vector(int *vect, int length);

extern void print_double_vector(double *vect, int lenght);

extern void print_double_matrix(double **matrix, int rows, int cols);

#endif