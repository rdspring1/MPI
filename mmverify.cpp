#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "mpi_matrix.h"

int main( int argc, char **argv )
{
	MPI_Init (&argc, &argv); 

	const size_t MATRIX_SIZE = atoi(argv[1]);
	const size_t DEPTH = atoi(argv[2]);

	int num_procs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	const int gridDim = std::sqrt(num_procs / DEPTH);

	int* A = NULL;
	int* B = NULL;
	double t1, t2;
	if(!myrank)
	{
		printf("N=%lu P=%d C=%lu\n", MATRIX_SIZE, num_procs, DEPTH);
		t1 = MPI_Wtime();
		A = random_matrix(MATRIX_SIZE);
		B = random_matrix(MATRIX_SIZE);
		//A = count_matrix(MATRIX_SIZE, gridDim);
		//B = count_matrix(MATRIX_SIZE, gridDim);
	}

	int* result = mpi_matrix_multiplication(A, B, MATRIX_SIZE, DEPTH);
	if(result) 
	{
		t2 = MPI_Wtime();
		printf("Elapsed Time: %f\n", t2-t1);
		(check(A, B, result, MATRIX_SIZE, gridDim)) ? printf("PASS\n") : printf("FAIL\n"); 
	}

	MPI_Finalize();
	return 0;
}
