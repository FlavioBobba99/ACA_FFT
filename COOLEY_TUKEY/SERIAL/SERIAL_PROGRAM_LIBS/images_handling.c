#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>

#include "images_handling.h"
#include "matrix_utilities.h"
#include "FFTs.h"

/*
    This function reads a PPM (Portable Pixmap) file and returns an Image structure
    with its width, height, max color value, and color matrices for red, green, and blue channels.
*/

Image *read_ppm(const char *filename) {
    FILE *file = fopen(filename, "rb");  // Open the file in binary mode
    if (file == NULL) {  // Error handling if the file cannot be opened
        fprintf(stderr, "Unable to open file %s\n", filename);
        exit(1);
    }

    // Allocate memory for the Image structure
    Image *img = (Image *)malloc(sizeof(Image));
    char format[3];
    fscanf(file, "%s", format);  // Read the file format (should be "P6")

    // Check that the file is in the correct PPM format (P6)
    if (format[0] != 'P' || format[1] != '6') {
        fprintf(stderr, "Invalid PPM file format (must be 'P6')\n");
        fclose(file);
        exit(1);
    }

    // Read the image dimensions (width, height) and max color value
    fscanf(file, "%d %d", &img->width, &img->height);
    fscanf(file, "%d", &img->max_color);
    fgetc(file);  // Consume the newline character after the max color value

    // Allocate memory for the color channels (red, green, blue)
    img->red = allocate_matrix(img->width, img->height);
    img->green = allocate_matrix(img->width, img->height);
    img->blue = allocate_matrix(img->width, img->height);

    // Read the pixel values for the three color channels (RGB)
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            img->red[i][j] = (double)fgetc(file);  // Read red pixel value
            img->green[i][j] = (double)fgetc(file);  // Read green pixel value
            img->blue[i][j] = (double)fgetc(file);  // Read blue pixel value
        }
    }

    fclose(file);  // Close the file after reading
    return img;  // Return the populated Image structure
}

/*
    Function used to free the memory allocated for the Image structure and its matrices.
*/

void free_image(Image *img) {
    free_matrix(img->red, img->height);  // Free the red channel matrix
    free_matrix(img->green, img->height);  // Free the green channel matrix
    free_matrix(img->blue, img->height);  // Free the blue channel matrix
    free(img);  // Free the Image structure
}

/*
    Function used to print the details of an Image structure (dimensions, color matrices).
*/

void print_image(Image *img) {
    printf("Width: %d, Height: %d, Max Color: %d\n", img->width, img->height, img->max_color);
    
    // Print the red channel matrix
    printf("Red Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->red[i][j]);
        }
        printf("\n");
    }

    // Print the green channel matrix
    printf("Green Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->green[i][j]);
        }
        printf("\n");
    }

    // Print the blue channel matrix
    printf("Blue Matrix:\n");
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            printf("%3f ", img->blue[i][j]);
        }
        printf("\n");
    }
}

/*
    This function shifts the quadrants of the input matrix to re-center the image 
    for frequency analysis (FFT shifting).
*/

void fftshift(double** input, double** output, int height, int width) {
    int half_height = height / 2;
    int half_width = width / 2;

    // Swap quadrants of the matrix
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
    This function scales the pixel values in the image based on a scale factor
    and writes the resulting image to a PPM file.
*/

void writePPM(const char *filename, Image *img, float scale_factor) {
    // Allocate memory for the shifted color matrices
    double **shifted_red = allocate_matrix(img->width, img->height);
    double **shifted_green = allocate_matrix(img->width, img->height);
    double **shifted_blue = allocate_matrix(img->width, img->height);

    // Perform FFT shift on each color channel
    fftshift(img->red, shifted_red, img->height, img->width);
    fftshift(img->green, shifted_green, img->height, img->width);
    fftshift(img->blue, shifted_blue, img->height, img->width);

    FILE *fp = fopen(filename, "wb");  // Open the output file in binary mode
    if (!fp) {  // Error handling if the file cannot be opened
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    // Write the PPM header
    fprintf(fp, "P6\n%d %d\n%d\n", img->width, img->height, img->max_color);

    // Write the pixel data (scaled) to the file
    for (int i = 0; i < img->height; i++) {
        for (int j = 0; j < img->width; j++) {
            unsigned char pixel[3];
            pixel[0] = (unsigned char)(shifted_red[i][j] * scale_factor);  // Red channel
            pixel[1] = (unsigned char)(shifted_green[i][j] * scale_factor);  // Green channel
            pixel[2] = (unsigned char)(shifted_blue[i][j] * scale_factor);  // Blue channel
            fwrite(pixel, sizeof(unsigned char), 3, fp);  // Write the pixel to the file
        }
    }

    fclose(fp);  // Close the file after writing
}

/*
    This function calculates the next power of 2 greater than or equal to the target value.
    It finds the smallest power of 2 that can contain the target value.
*/

int closest_square(int target){
    int square = 1;

    // Keep doubling the square value until it exceeds or equals the target
    while(square < target){
        square *= 2;
    }

    return square;  // Return the next power of 2 that is greater than or equal to the target
}

/*
    This function finds the scale factor for an image.
    The scale factor is calculated as the ratio of 255 to the maximum pixel value in the image.
    This is useful for scaling the pixel values to the 8-bit range (0-255).
*/

double find_scale_factor(Image *img){
    double scale = 1.0;
    double max_value = 0;

    // Check the red channel for the maximum pixel value
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->red[i][j] > max_value){
                max_value = img->red[i][j];
            }
        }
    }

    // Check the green channel for the maximum pixel value
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->green[i][j] > max_value){
                max_value = img->green[i][j];
            }
        }
    }

    // Check the blue channel for the maximum pixel value
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            if(img->blue[i][j] > max_value){
                max_value = img->blue[i][j];
            }
        }
    }

    // Calculate the scale factor to fit the maximum pixel value into the 8-bit range
    scale = (255 / max_value);
    
    return scale;  // Return the calculated scale factor
}

/*
    This function pads the input image to the next power of 2 for both width and height.
    This is typically done to optimize performance for Fourier transforms (FFT).
    The padding is done by allocating a new image with zeros in the padded areas.
*/

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

    // Allocate matrices for the color channels (red, green, blue) in the padded image
    padded_image->red = allocate_matrix(new_width, new_height);
    padded_image->green = allocate_matrix(new_width, new_height);
    padded_image->blue = allocate_matrix(new_width, new_height);

    // Initialize the padded matrices with zero and copy the original image data
    for (int i = 0; i < new_height; i++) {
        for (int j = 0; j < new_width; j++) {
            if (i < input_image->height && j < input_image->width) {
                // Copy original pixel values from the input image to the padded image
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

    return padded_image;  // Return the padded image
}

/*
    This function applies a logarithmic scale to the pixel values of the image.
    It enhances low-magnitude pixel values and compresses high-magnitude values for better visualization.
    Typically used after performing an FFT or for contrast enhancement.
*/

Image* log_scale (Image* img){
    // Allocate memory for the new image after applying log scaling
    Image* log_out = (Image *)malloc(sizeof(Image));

    log_out->width = img->width;
    log_out->height = img->height;
    log_out->max_color = img->max_color;

    // Allocate matrices for the color channels (red, green, blue) in the log-scaled image
    log_out->red = allocate_matrix(img->width, img->height);
    log_out->green = allocate_matrix(img->width, img->height);
    log_out->blue = allocate_matrix(img->width, img->height);

    // Apply log scaling to the red channel
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->red[i][j] = 20 * log10(img->red[i][j]);
        }
    }

    // Apply log scaling to the green channel
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->green[i][j] = 20 * log10(img->green[i][j]);
        }
    }

    // Apply log scaling to the blue channel
    for(int i = 0; i < img->height; i++){
        for(int j = 0; j < img->width; j++){
            log_out->blue[i][j] = 20 * log10(img->blue[i][j]);
        }
    }

    return log_out;  // Return the log-scaled image
}
