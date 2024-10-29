#ifndef FFTs_H
#define FFTs_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "images_handling.h"

extern double complex* FFT_complex(double complex *provavett, int lengthvett);
extern void transpose(complex double **matrix, complex double **result, int widht, int height);
extern double complex **matrix_FFT (double complex **matrix, int width, int height);
extern void FFT_image(Image *in, Image *module, Image *phase);
extern void FFT_complex_with_range(double complex *provavett, int starting_index, double complex *output_vector, int lenghtvett);

#endif