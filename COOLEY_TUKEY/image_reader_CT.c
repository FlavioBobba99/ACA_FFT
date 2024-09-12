#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

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
            pixel[0] = (unsigned char)img->red[i][j]*scale_factor;
            pixel[1] = (unsigned char)img->green[i][j]*scale_factor;
            pixel[2] = (unsigned char)img->blue[i][j]*scale_factor;
            fwrite(pixel, sizeof(unsigned char), 3, fp);
        }
    }

    fclose(fp);
}

// Function to free the allocated memory
void free_matrix(double **matrix, int height) {
    for (int i = 0; i < height; i++) {
        free(matrix[i]);
    }
    free(matrix);
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
    free_matrix(img->green, img->height);
    free_matrix(img->blue, img->height);
    free(img);
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

/*double complex* FFT(double *provavett, int lengthvett) {
    // Base case: lengthvett == 1
    if (lengthvett == 1) {
        double complex *single_out = malloc(sizeof(double complex));
        if (single_out == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        single_out[0] = provavett[0] + 0.0 * I;
        return single_out;
    }

    
        Malloc must be used in this case because allocating 
        a variable like a pointer in the normal way, inside a nested
        call, when the function terminates the variable is deallocated 
        casing undefined behaviour
    

    int split_length = lengthvett / 2;
    double theta = (2 * PI) / lengthvett;
    double complex w = cos(theta) + I * sin(theta);

    // Allocate memory for even and odd arrays
    double *even = malloc(split_length * sizeof(double));
    double *odd = malloc(split_length * sizeof(double));
    if (even == NULL || odd == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < split_length; i++) {
        even[i] = provavett[2 * i];
        odd[i] = provavett[2 * i + 1];
    }

    // Recursively compute FFT for even and odd parts
    double complex *y_even = FFT(even, split_length);
    double complex *y_odd = FFT(odd, split_length);

    // Allocate memory for output
    double complex *out = malloc(lengthvett * sizeof(double complex));
    if (out == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < split_length; i++) {
        out[i] = y_even[i] + cpow(w, i) * y_odd[i];
        out[i + split_length] = y_even[i] - cpow(w, i) * y_odd[i];
    }

    // Free allocated memory for even and odd arrays
    free(even);
    free(odd);
    free(y_even);
    free(y_odd);

    return out;
}*/

// Unecessary, deactivated
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

void transpose(complex double **matrix, complex double **result, int widht, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            result[j][i] = matrix[i][j];  // Transpose operation
        }
    }
}

/*
    This function calculates the FFT of a complex matrix in input
    and returns a complex matrix. It uses the previously defined complex 
    matrix
*/

double complex **matrix_FFT (double complex **matrix, int width, int height) { 

	double complex **temporary_matrix = allocate_complex_matrix(width, height);
    double complex **temporary_transposed_matrix = allocate_complex_matrix(height,width);
    double complex **out_transposed = allocate_complex_matrix(height,width);
    double complex **out_matrix = allocate_complex_matrix(width,height);
    
	for (int i = 0; i < height; i++) {
            double complex *out_vect = FFT_complex(matrix[i], width);
            for(int j = 0; j< width; j++){
                temporary_matrix[i][j] = out_vect[j];
            }
            free(out_vect);
		}

    transpose(temporary_matrix, temporary_transposed_matrix, width, height);
    free_complex_matrix(temporary_matrix, height);


    for (int i = 0; i < width; i++) {
            double complex *out_vect = FFT_complex(temporary_transposed_matrix[i], height);
            for(int j = 0; j < height; j++){
                out_transposed[i][j] = out_vect[j];
            }
            free(out_vect);
		}

    transpose(out_transposed, out_matrix, height, width);

    free_complex_matrix(temporary_transposed_matrix, width);
    free_complex_matrix(out_transposed, width);
	
	return out_matrix;
	}

/*
double complex **convert_to_complex_matrix(double **matrix, int rows, int cols) {
    // Allocate memory for the complex matrix
    double complex **complexMatrix = (double complex **)malloc(rows * sizeof(double complex *));
    for (int i = 0; i < rows; i++) {
        complexMatrix[i] = (double complex *)malloc(cols * sizeof(double complex));
    }

    // Convert double matrix to complex matrix (real part is the original value, imaginary part is 0)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            complexMatrix[i][j] = matrix[i][j] + 0 * I; // Real part is from the input matrix, Imaginary part is 0
        }
    }

    return complexMatrix;
}

void FFT_image(Image *in, Image *module, Image *phase){
    
    int height = in->height;
    int width = in->width;
    double complex **temp_red = convert_to_complex_matrix(in->red,height,width);
    printf("Matrix converted\n");
    double complex **complex_red = matrix_FFT(temp_red, width, height);
    print_complex_matrix(complex_red, height, width);
}*/

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

void FFT_image(Image *in, Image *module, Image *phase){
    int height = in->height;
    int width = in->width;

    // Check if 'in' has valid pointers for 'red'
    if (in->red == NULL) {
        fprintf(stderr, "Input image red channel is NULL\n");
        return;
    }

    double complex **temp_red = convert_to_complex_matrix(in->red, height, width);
    printf("Matrix converted\n");

    // Check if matrix_FFT handles memory properly
    double complex **complex_red = matrix_FFT(temp_red, width, height);
    if (complex_red == NULL) {
        fprintf(stderr, "FFT computation failed\n");
        return;
    }

    print_complex_matrix(complex_red, height, width);

    // Free allocated memory (example, modify as necessary)
    for (int i = 0; i < height; i++) {
        free(temp_red[i]);
    }
    free(temp_red);

    // Free complex_red if necessary (depends on the implementation of matrix_FFT)
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <image_path>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    Image *img = read_ppm(filename);
    
    //char new_image[] = "/home/flavio/Desktop/ferretti_MPI/IMAGES/output/test_out1.ppm";
    //double scale_factor = 0.2;
    //writePPM(new_image, img, scale_factor);

   // print_image(img);  // For debugging, can be removed if not needed

    
  /*  int  n = 8;
    double data[] = {1, 2, 3, 4, 5, 6, 7, 8};
    double complex *result = FFT(data, n);

    for (int i = 0; i < n; i++) {
        printf("FFT[%d] = %f + %fi\n", i, creal(result[i]), cimag(result[i]));
    }

    free(result);*/

    int rows = 4, cols = 8;

    // Dynamically allocate memory for the matrix
    double complex **matrix = malloc(rows * sizeof(double complex*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(double complex));
    }

    // Initialize the matrix with some values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = i + j * I;
        }
    }

    double complex **result_matrix = matrix_FFT(matrix, 8,4);

    print_complex_matrix(result_matrix, 4, 8);

    FFT_image(img,img,img);
    free_image(img);
    
    return 0;
}
