#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "mpi_matrix.h"

int main( int argc, char **argv )
{
	MPI_Init (&argc, &argv); 
	const size_t MATRIX_SIZE = 100;
	const size_t DEPTH = 2;
	int num_procs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	const int gridDim = std::sqrt(num_procs / DEPTH);

	int* A = NULL;
	int* B = NULL;
	if(!myrank)
	{
		A = random_matrix(MATRIX_SIZE);
		B = random_matrix(MATRIX_SIZE);
		//A = count_matrix(MATRIX_SIZE, gridDim);
		//B = count_matrix(MATRIX_SIZE, gridDim);
	}

	int* result = mpi_matrix_multiplication(A, B, MATRIX_SIZE, DEPTH);
	if(result) 
	{
		(check(A, B, result, MATRIX_SIZE, gridDim)) ? printf("\nPASS\n") : printf("\nFAIL\n"); 
	}

	MPI_Finalize();
	return 0;
}
