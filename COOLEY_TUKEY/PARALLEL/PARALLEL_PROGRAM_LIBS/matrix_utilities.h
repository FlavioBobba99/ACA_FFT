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

extern void print_complex_vector(double complex *input_vector, int lenght_vector);

extern void print_int_vector(int *vect, int length);

extern void print_double_vector(double *vect, int lenght);

extern void print_double_matrix(double **matrix, int rows, int cols);

extern double complex **unflatten_complex_matrix(double complex *in, int width, int height);

extern double **unflatten_double_matrix(double *in, int width, int height);

extern double complex **transpose_complex_matrix(double complex **in, int width, int height);

extern double **transpose_double_matrix(double **in, int width, int height);

extern double complex *flatten_complex_matrix(double complex **matrix, int width, int height);

extern double *flatten_double_matrix(double **matrix, int heigth, int width);

extern double complex *double_to_complex_vector(double *input_vector, int lenght_vector);

extern double find_max_in_double_vector(double *in, int lenght);

extern double find_max_in_double_vector_and_logscale(double *in, int lenght);

#endif
