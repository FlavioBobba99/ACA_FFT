#define _GNU_SOURCE
#define TEST_MATRIX_HEIGTH 4
#define TEST_MATRIX_WIDTH 4
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#include "./PARALLEL_PROGRAM_LIBS/images_handling.h"
#include "./PARALLEL_PROGRAM_LIBS/matrix_utilities.h"
#include "./PARALLEL_PROGRAM_LIBS/FFTs.h"

double *flatten_double_matrix(double **matrix, int heigth, int width){

    int v_index = 0;

    double *flat_matrix = malloc(heigth*width*sizeof(double));

    for(int i = 0; i < heigth; i++){
        for(int j = 0; j < width; j++){
           flat_matrix[v_index] = matrix[i][j];
           v_index++;
        }
    }
    return flat_matrix;
}

double complex *double_to_complex_vector(double *input_vector, int lenght_vector){
	double complex *output_vector = malloc(lenght_vector * sizeof(double complex));
	
	if (output_vector == NULL) {
        perror("Failed to allocate memory for complex vector");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < lenght_vector; i ++){
		output_vector[i] = input_vector[i] + 0 * I;
	}
	print_complex_vector(output_vector, lenght_vector);
	
	return output_vector;
	}

void scatter_and_flatten_double_matrix(double **matrix, int height, int width, int rank, int size, double *local_chunk){

    int* elements_per_process = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));

    int baseline_elements = height*width/size;
    int remaining_elements = height*width%size;

    int elements_accumulator = 0;

    double *flat_matrix;

    for(int i = 0; i < size; i++){
        //ERROR IN COUNTIING THE ELEMENTS 
        elements_per_process[i] = (i == size -1) ? baseline_elements + remaining_elements : baseline_elements;
        //POSSIBLE ERROR IN -1
        displacements[i] = (i == 0) ? 0 : elements_accumulator;
        elements_accumulator += elements_per_process[i]; 
    }

    
    if (rank == 0){
		
		printf("Elements per core vector:\n");
		print_int_vector(elements_per_process, size);

		printf("Displacements vector:\n");
		print_int_vector(displacements, size);
			
        printf("Process %d reached flattening stage\n", rank); 
        flat_matrix = flatten_double_matrix(matrix, height, width);
        printf("Debug flattened matrix:\n");
        print_double_vector(flat_matrix, height*width);
    }
    

    local_chunk = malloc(elements_per_process[rank]*sizeof(double));

    MPI_Scatterv(flat_matrix, elements_per_process, displacements, MPI_DOUBLE,
                 local_chunk, elements_per_process[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    
    printf("- - - - - - PROCESS %d - - - - - - - - - - - -\n", rank);
    printf("Elements of process %d has to handle: %d\n", rank, elements_per_process[rank]);
    print_double_vector(local_chunk, elements_per_process[rank]);
    printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  \n");
    
    double_to_complex_vector(local_chunk, elements_per_process[rank]);
    
    // Allocate space for the full matrix only on the root process
	double complex **gathered_matrix = NULL;
	
	if (rank == 0){
		gathered_matrix = allocate_complex_matrix(width, height);
		printf("gathered matrix allocated \n");
	}
	
	// Gather the data from all processes to process 0
	MPI_Gatherv(local_chunk, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
            gathered_matrix, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	
	if (rank == 0){
		printf("Gathered Matrix\n");
		print_complex_matrix(gathered_matrix, height, width);
	}
		
}

int main(int argc, char *argv[]) {

    if (argc != 4) {
        fprintf(stderr, "Correct usage: %s <image_path> <module_output_path> <phase_output_path>\n", argv[0]);
        printf("Arguments detected %d\n", argc);
        return 1;
    }

    int rank, size;

    Image *img = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    
    double **test_matrix = allocate_matrix(TEST_MATRIX_WIDTH,TEST_MATRIX_HEIGTH);

    int acc = 0;

    for(int i  = 0; i < TEST_MATRIX_HEIGTH; i++){
        for(int j = 0; j < TEST_MATRIX_WIDTH; j++){
            test_matrix[i][j] = acc;
            acc++;
        }
    }

    if(rank == 0){
        printf("Original double matrix:\n");
        print_double_matrix(test_matrix, TEST_MATRIX_HEIGTH, TEST_MATRIX_WIDTH);
    }

    double *local_chunk = NULL;
    scatter_and_flatten_double_matrix(test_matrix, TEST_MATRIX_HEIGTH, TEST_MATRIX_WIDTH, rank, size, local_chunk);
    

    MPI_Finalize();

    return 0;
    
    //fede
    //mpirun -n 3 ../parallel /home/fede/Documenti/ACA_FFT/IMAGES/blue16x16.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/out1.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/out2.ppm

}

