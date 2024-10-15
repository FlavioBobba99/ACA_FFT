
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#define PI 3.14159265358979323846

/*
    Structure to hold the image data. This structure
    utilizes 3 matrices to memorize the RGB values of each image.
    This is done to ensure contiguity in the storage of the different 
    vector that make up a matrix. This avoids storing irrelevant values 
    in cache while performing transformations.
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

// Function to free the matrix memory
void free_complex_matrix(double complex **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Example FFT function on a vector
double complex *FFT_complex(double complex *vector, int size) {
    // Placeholder FFT function; you should replace this with a real FFT implementation
    double complex *result = malloc(size * sizeof(double complex));
    for (int i = 0; i < size; i++) {
        result[i] = vector[i];  // A dummy transformation
    }
    return result;
}

// Parallel matrix FFT function using MPI
double complex **matrix_FFT(double complex **matrix, int width, int height, int rank, int size) {
    // Allocate memory for transposed and output matrices
    double complex **temporary_transposed_matrix = allocate_complex_matrix(height, width);
    double complex **out_matrix = allocate_complex_matrix(width, height);

    // Arrays for Scatterv and Gatherv
    int *sendcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));

    // Calculate sendcounts and displacements for row-wise scatter
    int base_row_count = height / size; // Base number of rows per process
    int extra_rows = height % size;     // Extra rows for uneven distribution

    for (int i = 0; i < size; i++) {
        sendcounts[i] = (i < extra_rows) ? (base_row_count + 1) * width : base_row_count * width;
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
    }

    // Allocate memory for local row chunks for each process
    int local_row_count = sendcounts[rank] / width;
    double complex *local_rows = malloc(sendcounts[rank] * sizeof(double complex));

    // Scatter the rows of the matrix to different processes
    MPI_Scatterv(&matrix[0][0], sendcounts, displs, MPI_C_DOUBLE_COMPLEX, local_rows,
                 sendcounts[rank], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform FFT on the local rows
    for (int i = 0; i < local_row_count; i++) {
        double complex *out_vect = FFT_complex(&local_rows[i * width], width);
        for (int j = 0; j < width; j++) {
            temporary_transposed_matrix[j][displs[rank] / width + i] = out_vect[j];
        }
        free(out_vect);
    }

    // Gather the transposed matrix back from all processes
    MPI_Gatherv(&temporary_transposed_matrix[0][0], sendcounts[rank], MPI_C_DOUBLE_COMPLEX,
                &temporary_transposed_matrix[0][0], sendcounts, displs, MPI_C_DOUBLE_COMPLEX,
                0, MPI_COMM_WORLD);

    // Now scatter columns for the second pass (columns FFT)
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (i < extra_rows) ? (base_row_count + 1) * height : base_row_count * height;
        displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
    }

    local_row_count = sendcounts[rank] / height;
    double complex *local_cols = malloc(sendcounts[rank] * sizeof(double complex));

    MPI_Scatterv(&temporary_transposed_matrix[0][0], sendcounts, displs, MPI_C_DOUBLE_COMPLEX, local_cols,
                 sendcounts[rank], MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Perform FFT on the columns
    for (int i = 0; i < local_row_count; i++) {
        double complex *out_vect = FFT_complex(&local_cols[i * height], height);
        for (int j = 0; j < height; j++) {
            out_matrix[displs[rank] / height + i][j] = out_vect[j];
        }
        free(out_vect);
    }

    // Gather the final output matrix back from all processes
    MPI_Gatherv(&out_matrix[0][0], sendcounts[rank], MPI_C_DOUBLE_COMPLEX,
                &out_matrix[0][0], sendcounts, displs, MPI_C_DOUBLE_COMPLEX,
                0, MPI_COMM_WORLD);

    // Clean up
    free(sendcounts);
    free(displs);
    free(local_rows);
    free(local_cols);
    free_complex_matrix(temporary_transposed_matrix, height);

    return out_matrix;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int width = 8;  // Example width
    int height = 6; // Example height

    double complex **matrix = NULL;
    if (rank == 0) {
        matrix = allocate_complex_matrix(width, height);
        // Initialize your matrix here
    }

    double complex **result = matrix_FFT(matrix, width, height, rank, size);

    // Only the root process will print the result
    if (rank == 0) {
        // Print the result matrix here
        free_complex_matrix(matrix, height);
    }

    MPI_Finalize();
    return 0;
}