#include <stdio.h>
#include <stdlib.h>

/*
    Defining a new type using the Pixel data structure which is composed of 3 values
*/ 

typedef struct {
    unsigned short int r, g, b;
} Pixel;

/* 
    Defining a new type using the data structure of pixel created beforehand.
    This data structure contains the width the height and the highest intensity
    of color registred in the image.
    the *data pointer it is used to store the address of a matrix that will contain the 
    pixel image in an accessible format.
*/

typedef struct {
    int width;
    int height;
    int max_color;
    Pixel *data;
} Image;

Image* read_ppm(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open file %s\n", filename);
        return NULL;
    }

    // Read the file header
    char format[3];
    fscanf(fp, "%2s", format);

    if (format[0] != 'P' || format[1] != '6') {
        fprintf(stderr, "Invalid PPM file format (must be 'P6')\n");
        fclose(fp);
        return NULL;
    }

    // Read image size information
    Image *img = (Image *)malloc(sizeof(Image));
    if (!img) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return NULL;
    }

    // Skip comments
    char ch;
    while ((ch = fgetc(fp)) == '#') {
        while (fgetc(fp) != '\n');
    }
    ungetc(ch, fp);

    fscanf(fp, "%d %d", &img->width, &img->height);
    fscanf(fp, "%d", &img->max_color);
    fgetc(fp);  // Skip the single whitespace character after max_color

    // Allocate memory for pixel data
    img->data = (Pixel *)malloc(img->width * img->height * sizeof(Pixel));
    if (!img->data) {
        fprintf(stderr, "Memory allocation failed for image data\n");
        free(img);
        fclose(fp);
        return NULL;
    }

    // Read the pixel data from the file
    fread(img->data, sizeof(Pixel), img->width * img->height, fp);

    fclose(fp);
    return img;
}

void free_image(Image *img) {
    if (img) {
        free(img->data);
        free(img);
    }
}

void print_rgb_matrices(Image *img) {
    // Crea le matrici per R, G e B
    unsigned int **R = (unsigned int **)malloc(img->height * sizeof(unsigned int *));
    unsigned int **G = (unsigned int **)malloc(img->height * sizeof(unsigned int *));
    unsigned int **B = (unsigned int **)malloc(img->height * sizeof(unsigned int *));
    
    for (int i = 0; i < img->height; i++) {
        R[i] = (unsigned int *)malloc(img->width * sizeof(unsigned int));
        G[i] = (unsigned int *)malloc(img->width * sizeof(unsigned int));
        B[i] = (unsigned int *)malloc(img->width * sizeof(unsigned int));
    }

    // Popola le matrici con i valori R, G e B
    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            Pixel pixel = img->data[y * img->width + x];
            R[y][x] = pixel.r;
            G[y][x] = pixel.g;
            B[y][x] = pixel.b;
        }
    }

    // Stampa le matrici
    printf("Matrice R:\n");
    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            printf("%3u ", R[y][x]);
        }
        printf("\n");
    }

    printf("\nMatrice G:\n");
    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            printf("%3u ", G[y][x]);
        }
        printf("\n");
    }

    printf("\nMatrice B:\n");
    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            printf("%3u ", B[y][x]);
        }
        printf("\n");
    }

    // Libera la memoria delle matrici
    for (int i = 0; i < img->height; i++) {
        free(R[i]);
        free(G[i]);
        free(B[i]);
    }
    free(R);
    free(G);
    free(B);
}

int main(int argc, char *argv[]) {
    const char *filename = argv[1];

    Image *img = read_ppm(filename);

    if (img) {
        printf("Image loaded: %dx%d, max color value: %d\n", img->width, img->height, img->max_color);
        
        // Stampa le matrici R, G e B
        print_rgb_matrices(img);

        // Libera la memoria dell'immagine
        free_image(img);
    }

    return 0;
}