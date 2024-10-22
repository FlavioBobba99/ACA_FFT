#ifndef PRINTING_H
#define PRINTING_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

/*
    Structure to holdi the image data. This structre
    utilizes 3 matrices to memorize the RGB values of each image.
    This is doone to ensure contiguity in the storage of the differnet 
    vector that make up a matrix. This is done because when values are used 
    i avoid storig in cache irrelevant values to perform the transformate.
*/

typedef struct {
    int width;
    int height;
    int max_color;
    double **red;
    double **green;
    double **blue;
} Image;

extern Image *read_ppm(const char *filename);

extern void free_image(Image *img);

extern void print_image(Image *img);

extern void writePPM(const char *filename, Image *img, float scale_factor);

extern void fftshift(double** input, double** output, int height, int width);

extern double find_scale_factor(Image *img);

extern Image* pad_image(const Image* input_image);

extern Image* log_scale (Image* img);

#endif