#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "mpi_matrix.h"

int main()
{
	size_t MATRIX_SIZE = 8;
	int** A = alloc_matrix(MATRIX_SIZE);
	int** B = alloc_matrix(MATRIX_SIZE);
	random_matrix(A, MATRIX_SIZE);
	random_matrix(B, MATRIX_SIZE);
	check(A, B, mpi_matrix_multiplication(A, B, MATRIX_SIZE), MATRIX_SIZE); 
	return 1;
}
