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

double *FFT(double *provavett, int lenghtvett) {
	if (lenghtvett ==1){
		return provavett;
	}
	double r= 1.0;
	double theta= (2*PI)/lenghtvett; 
	
	double complex w= cpow(r, theta);
	
	printf("%f i%f",creal(w), cimag(w));
	
	return provavett;
	
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

    print_image(img);  // For debugging, can be removed if not needed

    free_image(img);
    double test[4] = {1.0, 3.3, 2.0, PI};
    FFT(test, 4);
    return 0;
}
