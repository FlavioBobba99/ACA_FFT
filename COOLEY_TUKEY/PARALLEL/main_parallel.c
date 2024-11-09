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

// ATTENZIONE LA FUNIONE ACCETTA HEIGHT WIDTH PER MANTENERE L'ORDINE USATO NELLE CHIAMATE A FUNZIONE

void FFT_pt2(double complex **matrix, int height, int width, int rank, int size, double ***module_matrix, double ***phase_matrix, double *color_max){

	double func_start_time = MPI_Wtime();

    int* elements_per_process = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));

    int baseline_rows = height/size;
    int baseline_elements = baseline_rows*width;
    int remaining_elements = height%size*width;

    int elements_accumulator = 0;

    for(int i = 0; i < size; i++){
        elements_per_process[i] = (i == size -1) ? baseline_elements + remaining_elements : baseline_elements;
        displacements[i] = (i == 0) ? 0 : elements_accumulator;
        elements_accumulator += elements_per_process[i]; 
    }

    double complex *flat_complex_matrix = NULL;

    if (rank == 0){
        flat_complex_matrix = flatten_complex_matrix(matrix, width, height);
    }
    

    double complex *local_chunk = malloc(elements_per_process[rank]*sizeof(double complex));

    MPI_Scatterv(flat_complex_matrix, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX,
                 local_chunk, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
                 0, MPI_COMM_WORLD);

    double complex *local_FFT_output = malloc(elements_per_process[rank] * sizeof(double complex));

    for(int i = 0; i < elements_per_process[rank]; i += width){
       FFT_complex_with_range(local_chunk, i, local_FFT_output, width);
    }

    double *gathered_module = NULL;
    double *gathered_phase = NULL;
    double *gathered_max_values = NULL;

	
	if (rank == 0){

        gathered_module = malloc( width * height * sizeof(double)); 

        gathered_phase = malloc( width * height * sizeof(double)); 

        gathered_max_values = (double*)malloc(size * sizeof(double));
	}

    double *local_module = malloc(elements_per_process[rank]*sizeof(double));
    double *local_phase = malloc(elements_per_process[rank]*sizeof(double));
    
    for(int i = 0; i < elements_per_process[rank]; i++){
        local_module[i] = cabs(local_FFT_output[i]);
        local_phase[i] = carg(local_FFT_output[i]);
    }

    free(local_FFT_output);

    double local_max_module = find_max_in_double_vector_and_logscale(local_module, elements_per_process[rank]);

    MPI_Gatherv(local_module, elements_per_process[rank], MPI_DOUBLE,
            gathered_module, elements_per_process, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(local_phase, elements_per_process[rank], MPI_DOUBLE,
            gathered_phase, elements_per_process, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gather(&local_max_module, 1, MPI_DOUBLE, gathered_max_values, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);



    free(local_module);
    free(local_phase);
	
	if (rank == 0){

        *module_matrix = unflatten_double_matrix(gathered_module, width, height);
        *phase_matrix = unflatten_double_matrix(gathered_phase, width, height);

        *color_max = find_max_in_double_vector(gathered_max_values, size);

        *module_matrix = transpose_double_matrix(*module_matrix, width, height);
        *phase_matrix = transpose_double_matrix(*phase_matrix, width, height);

        free(gathered_module);
        free(gathered_phase);
        free(local_chunk);
        free(elements_per_process);
        free(displacements);
	}
	
	double func_end_time = MPI_Wtime();
    printf("Execution time of FFT_pt2: %f seconds\n", func_end_time - func_start_time);
}

void PARALLEL_FFT_matrix(double **matrix, int height, int width, int rank, int size, double ***matrix_module_out, double ***matrix_phase_out,double *color_max){

	double func_start_time = MPI_Wtime();

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
        flat_matrix = flatten_double_matrix(matrix, height, width);
    }
    

    double *local_chunk = malloc(elements_per_process[rank]*sizeof(double));

    MPI_Scatterv(flat_matrix, elements_per_process, displacements, MPI_DOUBLE,
                 local_chunk, elements_per_process[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    
    double complex *local_complex_chunk = double_to_complex_vector(local_chunk, elements_per_process[rank]);
    double complex *local_FFT_output = malloc(elements_per_process[rank] * sizeof(double complex));

    for(int i = 0; i < elements_per_process[rank]; i += width){
        FFT_complex_with_range(local_complex_chunk, i, local_FFT_output, width);
    }
    
	double complex *gathered_vector = NULL;
	
	if (rank == 0){
		gathered_vector = malloc( width * height * sizeof(double complex)); 
	}
	
	MPI_Gatherv(local_FFT_output, elements_per_process[rank], MPI_C_DOUBLE_COMPLEX,
            gathered_vector, elements_per_process, displacements, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    double complex **FFT_pt1 = NULL;
    double complex **FFT_pt1_transposed = NULL;
	
	if (rank == 0){
        FFT_pt1 = unflatten_complex_matrix(gathered_vector, width, height);

        FFT_pt1_transposed = transpose_complex_matrix(FFT_pt1, width, height);

        free_complex_matrix(FFT_pt1, height);
        free(local_complex_chunk);
        free(local_chunk);
        free(gathered_vector);
        free(local_FFT_output);
        free(elements_per_process);
        free(displacements);

	}
    FFT_pt2(FFT_pt1_transposed, width, height, rank, size, matrix_module_out, matrix_phase_out, color_max);
    
    double func_end_time = MPI_Wtime();
    printf("Execution time of PARALLEL_FFT_MATRIX: %f seconds\n", func_end_time - func_start_time);
    
}

void PARALLEL_image_FFT(Image *in, Image **module, Image **phase, int rank, int size, double *global_max){

	double func_start_time = MPI_Wtime();

    double **module_channel = NULL;
    double **phase_channel = NULL;
    double max_RGB[3] = {0};

    if(rank == 0){
        printf("PARALLEL FFT REACHED\n");
    }

    if(rank == 0){
        // Initialize the output module and phase images if rank is 0
        *module = (Image *)malloc(sizeof(Image));
        (*module)->width = in->width;
        (*module)->height = in->height;
        (*module)->max_color = 255;
        (*module)->red = NULL;   // Allocate only when needed
        (*module)->green = NULL;
        (*module)->blue = NULL;

        *phase = (Image *)malloc(sizeof(Image));
        (*phase)->width = in->width;
        (*phase)->height = in->height;
        (*phase)->max_color = 255;
        (*phase)->red = NULL;    // Allocate only when needed
        (*phase)->green = NULL;
        (*phase)->blue = NULL;
    }

    // Process RED channel
    PARALLEL_FFT_matrix(in->red, in->height, in->width, rank, size, &module_channel, &phase_channel, &max_RGB[0]);
    if(rank == 0) {
        (*module)->red = module_channel;
        (*phase)->red = phase_channel;
    }

    if(rank == 0){
        printf("RED CHANNEL DONE\n");
    }
    // Free module_channel and phase_channel after processing the red channel, to conserve memory
    module_channel = NULL;
    phase_channel = NULL;

    // Process GREEN channel
    PARALLEL_FFT_matrix(in->green, in->height, in->width, rank, size, &module_channel, &phase_channel, &max_RGB[1]);
    if(rank == 0) {
        (*module)->green = module_channel;
        (*phase)->green = phase_channel;
    }

    if(rank == 0){
        printf("GREEN CHANNEL DONE\n");
    }
    // Free module_channel and phase_channel after processing the green channel
    module_channel = NULL;
    phase_channel = NULL;

    // Process BLUE channel
    PARALLEL_FFT_matrix(in->blue, in->height, in->width, rank, size, &module_channel, &phase_channel, &max_RGB[2]);
    if(rank == 0) {
        (*module)->blue = module_channel;
        (*phase)->blue = phase_channel;
    }

    if(rank == 0){
        printf("BLUE CHANNEL DONE\n");
    }
    // No need to free here since we’re at the last color channel

    // Calculate global max RGB value
    if(rank == 0) {
        //print_double_vector(max_RGB, 3);
        double global_max_RGB = find_max_in_double_vector(max_RGB, 3);
        *global_max = global_max_RGB;
       // printf("The global RGB max is %f\n", global_max_RGB);
    }
    
    double func_end_time = MPI_Wtime();
    printf("Execution time of PARALLEL_image_FFT: %f seconds\n", func_end_time - func_start_time);
}


int main(int argc, char *argv[]) {

    if (argc != 4) {
        fprintf(stderr, "Correct usage: %s <image_path> <module_output_path> <phase_output_path>\n", argv[0]);
        printf("Arguments detected %d\n", argc);
        return 1;
    }

    int rank, size;

    Image *img = NULL;
    Image *phase = NULL;
    Image *module = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start_time = MPI_Wtime();

    int img_width, img_height;
    if (rank == 0){
        img = pad_image(read_ppm(argv[1]));
        if (img == NULL) {
            fprintf(stderr, "Error reading image.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        img_width = img->width;
        img_height = img->height;
    }

    MPI_Bcast(&img_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&img_width, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // If the image needs to be sent to other ranks, you might need to send img->data too
    if (rank != 0) {
        img = (Image *)malloc(sizeof(Image));
        img->width = img_width;
        img->height = img_height;
    }

    //printf("I am process %d and i have recieved width, height = %d, %d\n", rank, img->width, img->height);

    double max_RGB = 0;

    PARALLEL_image_FFT(img, &module, &phase, rank, size, &max_RGB);
   // printf("the max_RGB: %f\n",max_RGB);

    if(rank == 0){
        printf("SAVING IMAGE\n");
        double scale = (255 / max_RGB);
      //  printf("SCALE FACTOR = %f\n", scale);
		writePPM(argv[2], module, scale);

        scale = (255 / 2*M_PI);
        //printf("SCALE FACTOR = %f\n", scale);
		writePPM(argv[3], phase, scale);
    }
    
    double end_time = MPI_Wtime();

    MPI_Finalize();

    printf("RUN SUCCESSFUL\n");
    printf("Total execution time of main: %f seconds\n", end_time - start_time);
    return 0;
    
    //fede
    //mpirun -n 4 ../parallel /home/fede/Documenti/ACA_FFT/IMAGES/blue16x16.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/module.ppm /home/fede/Documenti/ACA_FFT/IMAGES/output/phase.ppm


    //flavio
    //mpirun -n 4 ../parallel /home/flavio/Desktop/ferretti_MPI/IMAGES/blue_4x4.ppm /home/flavio/Desktop/ferretti_MPI/IMAGES/output/module.ppm /home/flavio/Desktop/ferretti_MPI/IMAGES/output/phase.ppm
}

