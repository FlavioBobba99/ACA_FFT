#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "images_handling.h"
#include "matrix_utilities.h"
#include "FFTs.h"

/*
    This function reads ppm files and returns an image type data structure
*/

Image *read_ppm(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(1);
    }

    Image *img = (Image *)malloc(sizeof(Image));
    char format[3];
    fscanf(file, "%s", format);

    // Check if the file is a P6 PPM file
    if (format[0] != 'P' || format[1] != '6') {
        fprintf(stderr, "Invalid PPM file format (must be 'P6')\n");
        fclose(file);
        exit(1);
    }

    // Read image size and max color value
    fscanf(file, "%d %d", &img->width, &img->height);
    fscanf(file, "%d", &img->max_color);
    fgetc(file);  // Consume the newline character after the max color value

    // Allocate memory for the color matrices
    img->red = allocate_matrix(img->width, img->height);
    img->green = allocate_matrix(img->width, img->height);
    img->blue = allocate_matrix(img->width, img->height);

    // Read the pixel values into the matrices
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            img->red[i][j] = (double)fgetc(file);
            img->green[i][j] = (double)fgetc(file);
            img->blue[i][j] = (double)fgetc(file);
        }
    }

    fclose(file);
    return img;
}

/*
    Function used to free a IMAGE data structure
*/

void free_image(Image *img) {
    free_matrix(img->red, img->height);
    printf("FREE_MATRIX_IMG_DEBUG\n");
    free_matrix(img->green, img->height);
    printf("FREE_MATRIX_IMG_DEBUG\n");
    free_matrix(img->blue, img->height);
    printf("FREE_MATRIX_IMG_DEBUG\n");
    free(img);
    printf("FREE_IMG_DEBUG\n");
}

/*
    Function used to print to console a IMAGE data structure
*/

void print_image(Image *img) {
    printf("Width: %d, Height: %d, Max Color: %d\n", img->width, img->height, img->max_color);
    printf("Red Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->red[i][j]);
        }
        printf("\n");
    }

    printf("Green Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->green[i][j]);
        }
        printf("\n");
    }

    printf("Blue Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->blue[i][j]);
        }
        printf("\n");
    }
}

void fftshift(double** input, double** output, int height, int width) {
    int half_height = height / 2;
    int half_width = width / 2;

    // Shift each quadrant into its new position
    for (int i = 0; i < half_height; i++) {
        for (int j = 0; j < half_width; j++) {
            // Top-left to bottom-right
            output[i + half_height][j + half_width] = input[i][j];

            // Bottom-right to top-left
            output[i][j] = input[i + half_height][j + half_width];

            // Top-right to bottom-left
            output[i + half_height][j] = input[i][j + half_width];

            // Bottom-left to top-right
            output[i][j + half_width] = input[i + half_height][j];
        }
    }
}

/*
    The method must scale the int values before casting them
    ex 4000:255=current_value:X this gives a scale factor to
    which each pixel must be multiplied to in order to retain 
    scale. 
    The scale factor could be caluclated in the method or outside
    by passing the max value or the factori directly.
*/

void writePPM(const char *filename, Image *img, float scale_factor) {

    double **shifted_red = allocate_matrix(img->width, img->height);
    double **shifted_green= allocate_matrix(img->width, img->height);
    double **shifted_blue = allocate_matrix(img->width, img->height);

    fftshift(img->red, shifted_red, img->height, img->width);
    fftshift(img->green, shifted_green,img->height, img->width);
    fftshift(img->blue, shifted_blue, img->height, img->width);


    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    // Write the PPM header
    fprintf(fp, "P6\n%d %d\n%d\n", img->width, img->height, img->max_color);

    // Write pixel data
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            unsigned char pixel[3];
            pixel[0] = (unsigned char)(shifted_red[i][j]*scale_factor);
            pixel[1] = (unsigned char)(shifted_green[i][j]*scale_factor);
            //printf("PRE BLUE VALUE = %f \n", shifted_blue[i][j]*scale_factor);
            pixel[2] = (unsigned char)(shifted_blue[i][j]*scale_factor);
          //  printf("BLUE VALUE = %d \n", pixel[2]);
            fwrite(pixel, sizeof(unsigned char), 3, fp);
        }
    }

    fclose(fp);
}

int closest_square(int target){

    int square = 1;

    while(square<target){
        square *= 2;
    }
    return square;
}

double find_scale_factor(Image *img){

    double scale = 1.0;
    double max_value = 0;

    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->red[i][j]>max_value){
                max_value = img->red[i][j];
            }
        }
    }

     for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->green[i][j]>max_value){
                max_value = img->green[i][j];
            }
        }
    }

     for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->blue[i][j]>max_value){
                max_value = img->blue[i][j];
            }
        }
    }
    
    scale = (255 / max_value);
    printf("SCALE FACTOR = %f\n", scale);
    printf("MAX VALUE = %f\n", max_value);
    return scale;
}

Image* pad_image(const Image* input_image) {
    // Determine the new padded width and height (next power of 2)
    int new_width = closest_square(input_image->width);
    int new_height = closest_square(input_image->height);

    // Allocate memory for the new padded image
    Image* padded_image = (Image*)malloc(sizeof(Image));
    if (!padded_image) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // Set the dimensions of the padded image
    padded_image->width = new_width;
    padded_image->height = new_height;
    padded_image->max_color = input_image->max_color;

    // Allocate matrices for the color channels using your allocate_matrix function
    padded_image->red = allocate_matrix(new_width, new_height);
    padded_image->green = allocate_matrix(new_width, new_height);
    padded_image->blue = allocate_matrix(new_width, new_height);

    // Initialize the padded matrices with zero and copy the original image data
    for (int i = 0; i < new_height; i++) {
        for (int j = 0; j < new_width; j++) {
            if (i < input_image->height && j < input_image->width) {
                // Copy original pixel values
                padded_image->red[i][j] = input_image->red[i][j];
                padded_image->green[i][j] = input_image->green[i][j];
                padded_image->blue[i][j] = input_image->blue[i][j];
            } else {
                // Zero padding for the extra rows and columns
                padded_image->red[i][j] = 0.0;
                padded_image->green[i][j] = 0.0;
                padded_image->blue[i][j] = 0.0;
            }
        }
    }

    return padded_image;
}

Image* log_scale (Image* img){
    Image* log_out = (Image *)malloc(sizeof(Image));

    log_out->width = img->width;
    log_out->height = img->height;
    log_out->max_color = img->max_color;

    log_out->red = allocate_matrix(img->width, img->height);
    log_out->green = allocate_matrix(img->width, img->height);
    log_out->blue = allocate_matrix(img->width, img->height);

    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->red[i][j] = 20*log10(img->red[i][j]);
        }
    }

     for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->green[i][j] = 20*log10(img->green[i][j]);
        }
    }

     for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->blue[i][j] = 20*log10(img->blue[i][j]);
        }
    }
    printf("scale convertion done!\n");

    return log_out;
}