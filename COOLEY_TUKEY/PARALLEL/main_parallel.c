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

void FFT_pt2(double complex **matrix, int height, int width, int rank, int size, double **module_matrix, double **phase_matrix){

    int* elements_per_process = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));

    int baseline_rows = height/size;
    int baseline_elements = baseline_rows*width;
    int remaining_elements = height%size*width;

    int elements_accumulator = 0;

    for(int i = 0; i < size; i++){
        //ERROR IN COUNTIING THE ELEMENTS 
        elements_per_process[i] = (i == size -1) ? baseline_elements + remaining_elements : baseline_elements;
        //POSSIBLE ERROR IN -1
        displacements[i] = (i == 0) ? 0 : elements_accumulator;
        elements_accumulator += elements_per_process[i]; 
    }

    double complex *flat_complex_matrix = NULL;

    if (rank == 0){
		
		printf("Elements per core vector:\n");
		print_int_vector(elements_per_process, size);

		printf("Displacements vector:\n");
		print_int_vector(displacements, size);

        flat_complex_matrix = flatten_complex_matrix(matrix, width, height);
        printf("Flattened complex matrix from FFT PT2\n");
        print_complex_vector(flat_complex_matrix, width*height);
    }
    

    double complex *local_chunk = malloc(elements_per_process[rank]*sizeof(double complex));

    MPI_Scatterv(flat_complex_matrix, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX,
                 local_chunk, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
                 0, MPI_COMM_WORLD);
    
    printf("- - - - - - PROCESS %d - - - - - - - - - - - -\n", rank);
    printf("Elements of process %d has to handle: %d\n", rank, elements_per_process[rank]);
    print_complex_vector(local_chunk, elements_per_process[rank]);
    printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  \n");

    double complex *local_FFT_output = malloc(elements_per_process[rank] * sizeof(double complex));

    for(int i = 0; i < elements_per_process[rank]; i += width){
        printf("DEBUG INDEXES i = %d\n", i);
        FFT_complex_with_range(local_chunk, i, local_FFT_output, width);
    }
    printf("FFT OUPTUT TEST FROM CORE %d\n", rank);
    print_complex_vector(local_FFT_output, elements_per_process[rank]);
    
    // Allocate space for the full matrix only on the root process
	double complex *gathered_vector = NULL;
    double *gathered_module = NULL;
    double *gathered_phase = NULL;

	
	if (rank == 0){
		gathered_vector = malloc( width * height * sizeof(double complex)); 
		printf("gathered vector allocated \n");

        gathered_module = malloc( width * height * sizeof(double)); 
		printf("gathered module allocated \n");

        gathered_phase = malloc( width * height * sizeof(double)); 
		printf("gathered phase allocated \n");
	}

    double *local_module = malloc(elements_per_process[rank]*sizeof(double));
    double *local_phase = malloc(elements_per_process[rank]*sizeof(double));
    
    for(int i = 0; i < elements_per_process[rank]; i++){
        local_module[i] = cabs(local_FFT_output[i]);
        local_phase[i] = carg(local_FFT_output[i]);
    }


    printf("Local moudle test\n");
    print_double_vector(local_module, elements_per_process[rank]);
    printf("Local phase test\n");
    print_double_vector(local_phase, elements_per_process[rank]);

	// Gather the data from all processes to process 0
	MPI_Gatherv(local_FFT_output, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
            gathered_vector, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    MPI_Gatherv(local_module, elements_per_process[rank], MPI_DOUBLE,
            gathered_module, elements_per_process, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(local_phase, elements_per_process[rank], MPI_DOUBLE,
            gathered_phase, elements_per_process, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double complex **FFT_pt1 = NULL;
    double complex **FFT_pt1_transposed = NULL;
	
	if (rank == 0){
		printf("Gathered Vector\n");
		print_complex_vector(gathered_vector, height * width);
        printf("Gathered module\n");
        print_double_vector(gathered_module, width * height);
        printf("Gathered phase\n");
        print_double_vector(gathered_phase, width * height);

        module_matrix = unflatten_double_matrix(gathered_module, width, height);
        phase_matrix = unflatten_double_matrix(gathered_phase, width, height);

        print_double_matrix(module_matrix, height, width);
        print_double_matrix(phase_matrix, height, width);

        FFT_pt1 = unflatten_complex_matrix(gathered_vector, width, height);
        printf("FINAL RESULT\n");
        printf("UNFLATTENED MATRIX\n");
        print_complex_matrix(FFT_pt1, height, width);

        FFT_pt1_transposed = transpose_complex_matrix(FFT_pt1, width, height);

        print_complex_matrix(FFT_pt1_transposed, width, height);
        
	}
}

void PARALLEL_FFT_matrix(double **matrix, int height, int width, int rank, int size, double **matrix_module_out, double **matrix_phase_out){

    int* elements_per_process = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));

    int baseline_rows = height/size;
    int baseline_elements = baseline_rows*width;
    int remaining_elements = height%size*width;

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
    

    double *local_chunk = malloc(elements_per_process[rank]*sizeof(double));

    MPI_Scatterv(flat_matrix, elements_per_process, displacements, MPI_DOUBLE,
                 local_chunk, elements_per_process[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    
    printf("- - - - - - PROCESS %d - - - - - - - - - - - -\n", rank);
    printf("Elements of process %d has to handle: %d\n", rank, elements_per_process[rank]);
    print_double_vector(local_chunk, elements_per_process[rank]);
    printf("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  \n");
    
    double complex *local_complex_chunk = double_to_complex_vector(local_chunk, elements_per_process[rank]);
    double complex *local_FFT_output = malloc(elements_per_process[rank] * sizeof(double complex));

    for(int i = 0; i < elements_per_process[rank]; i += width){
        printf("DEBUG INDEXES i = %d\n", i);
        FFT_complex_with_range(local_complex_chunk, i, local_FFT_output, width);
    }
    printf("FFT OUPTUT TEST FROM CORE %d\n", rank);
    print_complex_vector(local_FFT_output, elements_per_process[rank]);
    
    // Allocate space for the full matrix only on the root process
	double complex *gathered_vector = NULL;
	
	if (rank == 0){
		gathered_vector = malloc( width * height * sizeof(double complex)); 
		//printf("gathered vector allocated \n");
	}
	
	// Gather the data from all processes to process 0
	MPI_Gatherv(local_FFT_output, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
            gathered_vector, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    double complex **FFT_pt1 = NULL;
    double complex **FFT_pt1_transposed = NULL;
	
	if (rank == 0){
		printf("Gathered Vector\n");
		print_complex_vector(gathered_vector, height * width);
        FFT_pt1 = unflatten_complex_matrix(gathered_vector, width, height);
        printf("FIRST FFT STEP\n\n");
        printf("UNFLATTENED MATRIX\n");
        print_complex_matrix(FFT_pt1, height, width);

        FFT_pt1_transposed = transpose_complex_matrix(FFT_pt1, width, height);
        printf("TRANSPOSED MATRIX\n");
        print_complex_matrix(FFT_pt1_transposed, width, height);

	}
    FFT_pt2(FFT_pt1_transposed, width, height, rank, size, matrix_module_out, matrix_phase_out);
}

void PARALLEL_image_FFT(Image *in, Image *module, Image *phase, int rank, int size){

    double **module_red = NULL;
    double **phase_red = NULL;

    PARALLEL_FFT_matrix(in->red, in->height,in->width, rank, size, module_red, phase_red);
    
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

    if (rank == 0){
        img = read_ppm(argv[1]);
    }
    
    
    double **test_matrix = allocate_matrix(TEST_MATRIX_WIDTH,TEST_MATRIX_HEIGTH);

    int acc = 0;

    for(int i  = 0; i < TEST_MATRIX_HEIGTH; i++){
        for(int j = 0; j < TEST_MATRIX_WIDTH; j++){
            test_matrix[i][j] = 255;
            acc++;
        }
    }

    if(rank == 0){
        printf("Original double matrix:\n");
        print_double_matrix(test_matrix, TEST_MATRIX_HEIGTH, TEST_MATRIX_WIDTH);
    }

    double *local_chunk = NULL;
    double **matrix_module_output = NULL;
    double **matrix_phase_output = NULL;
    PARALLEL_FFT_matrix(test_matrix, TEST_MATRIX_HEIGTH, TEST_MATRIX_WIDTH, rank, size, matrix_module_output, matrix_phase_output);
    

    MPI_Finalize();

    return 0;
    
    //fede
    //mpirun -n 3 ../parallel /home/fede/Documenti/ACA_FFT/IMAGES/blue16x16.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/out1.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/out2.ppm


    //flavio
    //mpirun -n 4 ../parallel /home/flavio/Desktop/ferretti_MPI/IMAGES/blue_4x4.ppm /home/flavio/Desktop/ferretti_MPI/IMAGES/output/module.ppm /home/flavio/Desktop/ferretti_MPI/IMAGES/output/phase.ppm
}

