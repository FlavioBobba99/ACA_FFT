#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#define PI 3.14159265358979323846

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

// Function to allocate memory for the matrices
double **allocate_matrix(int width, int height) {
    double **matrix = (double **)malloc(height * sizeof(double *));
    for (int i = 0; i < height; i++) {
        matrix[i] = (double *)malloc(width * sizeof(double));
    }
    return matrix;
}

// Function to allocate memory for the complex matrices
double complex **allocate_complex_matrix(int width, int height) {
    double complex **matrix = (double complex **)malloc(height * sizeof(double complex *));
    for (int i = 0; i < height; i++) {
        matrix[i] = (double complex *)malloc(width * sizeof(double complex));
    }
    return matrix;
}

void print_complex_matrix(double complex **matrix, int rows, int cols) {
    printf("Complex Matrix:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("(%6.2f + %6.2fi)  ", creal(matrix[i][j]), cimag(matrix[i][j]));
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

// TODO expand this method: 
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

// Function to free the allocated memory
void free_matrix(double **matrix, int height) {
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
        matrix[i] = NULL;  // Set to NULL after freeing
    }
    free(matrix);
    matrix = NULL;
}


// Function to free a 2D complex matrix
void free_complex_matrix(double complex **matrix, int height) {
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to read the PPM file and fill the Image structure
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


// Function to free the memory allocated for the Image structure
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

// Function to print the image data (for debugging)
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

double complex **matrix_FFT(double complex **matrix, int width, int height, int rank, int size) {
    // Allocate memory for the transposed and output matrices
    double complex **temporary_transposed_matrix = allocate_complex_matrix(height, width);
    double complex **out_matrix = allocate_complex_matrix(width, height);

    // Arrays for Scatterv and Gatherv
    int *sendcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));

    // Calculate sendcounts and displacements for row-wise scatter
    int base_row_count = height / size;  // Base number of rows per process
    int extra_rows = height % size;      // Extra rows for uneven distribution

    // Calculate sendcounts (in terms of rows) and displacements (in terms of rows)
    for (int i = 0; i < size; i++) {
        sendcounts[i] = ((i < extra_rows) ? (base_row_count + 1) : base_row_count) * width;  // Number of elements
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];  // Displacement in elements
    }

    // Allocate memory for local rows to be received by each process
    int local_row_count = sendcounts[rank]; // Number of rows each process will handle
    double complex **local_rows = allocate_complex_matrix(width, local_row_count);

    // Scatter full rows of the matrix to different processes
    MPI_Scatterv(&(matrix[0][0]), sendcounts, displs, MPI_C_DOUBLE_COMPLEX, &(local_rows[0][0]),
                 local_row_count * width, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform FFT on the local rows (vectors)
    for (int i = 0; i < local_row_count; i++) {
        //printf("Performing FFT of vector %d", i);
        double complex *out_vect = FFT_complex(local_rows[i], width);  // FFT on each row
        for (int j = 0; j < width; j++) {
            temporary_transposed_matrix[j][displs[rank] + i] = out_vect[j];  // Transpose and store the result
        }
        free(out_vect);  // Free the temporary vector
    }

    printf("NOW FINALIZING\n");

    // Gather the transposed matrix back from all processes
    MPI_Gatherv(&(temporary_transposed_matrix[0][0]), local_row_count * width, MPI_C_DOUBLE_COMPLEX,
                &(temporary_transposed_matrix[0][0]), sendcounts, displs, MPI_C_DOUBLE_COMPLEX,
                0, MPI_COMM_WORLD);

    printf("GATHER\n");

    // Recalculate sendcounts and displs for column-wise scatter
    for (int i = 0; i < size; i++) {
        sendcounts[i] = ((i < extra_rows) ? (base_row_count + 1) : base_row_count) * height;  // Number of elements
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];  // Displacement in elements
    }

    local_row_count = sendcounts[rank] / height;  // Update the number of rows to process for columns
    double complex **local_cols = allocate_complex_matrix(height, local_row_count);

    // Scatter columns (transposed matrix rows) to the processes for FFT
    MPI_Scatterv(&(temporary_transposed_matrix[0][0]), sendcounts, displs, MPI_C_DOUBLE_COMPLEX,
                 &(local_cols[0][0]), local_row_count * height, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform FFT on the transposed matrix columns (now scattered as rows)
    for (int i = 0; i < local_row_count; i++) {
        double complex *out_vect = FFT_complex(local_cols[i], height);  // FFT on each column
        for (int j = 0; j < height; j++) {
            out_matrix[displs[rank] / height + i][j] = out_vect[j];  // Store the result in the output matrix
        }
        free(out_vect);  // Free the temporary vector
    }

    printf("NOW FINALIZING\n");

    // Gather the final output matrix back from all processes
    MPI_Gatherv(&(out_matrix[0][0]), local_row_count * height, MPI_C_DOUBLE_COMPLEX,
                &(out_matrix[0][0]), sendcounts, displs, MPI_C_DOUBLE_COMPLEX,
                0, MPI_COMM_WORLD);
    printf("GATHER SUCCESSFUL\n");

    // Clean up
    free(sendcounts);
    free(displs);
    free_complex_matrix(local_rows, local_row_count);
    free_complex_matrix(local_cols, local_row_count);
    free_complex_matrix(temporary_transposed_matrix, height);

    return out_matrix;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void FFT_image(Image *in, Image *module, Image *phase, int rank, int size){

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
    double complex **complex_red = matrix_FFT(temp_red, width, height, rank, size);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT RED computation failed\n");
        return;
    }
    
    double complex **complex_green = matrix_FFT(temp_green, width, height, rank, size);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT GREEN computation failed\n");
        return;
    }

    double complex **complex_blue = matrix_FFT(temp_blue, width, height, rank, size);
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
    // Free complex_red if necessary (depends on the implementation of matrix_FFT)
}


int closest_square(int target){

    int square = 1;

    while(square<target){
        square *= 2;
    }
    return square;
}

void zero_padding(Image *in, Image *out){

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


int main(int argc, char *argv[]) {

    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <image_path>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    Image *img = read_ppm(filename);
    Image *module;
    Image *phase;

    module = (Image *)malloc(sizeof(Image));
    phase = (Image *)malloc(sizeof(Image));

    Image *padded_image = pad_image(img);
    
    char module_out[] = "/home/flavio/Desktop/ferretti_MPI/IMAGES/output/module_out_test_1.ppm";
    char phase_out[] = "/home/flavio/Desktop/ferretti_MPI/IMAGES/output/phase_out_test_1.ppm";

    FFT_image(padded_image,module,phase, rank, size);

    //writePPM(module_out, module, find_scale_factor(module));
    writePPM(phase_out, phase, find_scale_factor(phase));

    Image* module_log = log_scale(module);
    writePPM(module_out, module_log, find_scale_factor(module_log));

    printf("Image trnasformed and saved!\n");

    printf("Freeing original image\n");
    free_image(img);
    printf("Freeing padded image\n");
    free_image(padded_image);

    printf("Freeing module image\n");
    free_image(module);
    printf("Freeing phase image\n");
    free_image(phase);

    MPI_Finalize();
    
    return 0;
}
