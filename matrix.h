#ifndef MATRIX_H_
#define MATRIX_H_

#include <ctime>
#include <random>
#include <algorithm>

void print_matrix(int** matrix, const size_t size) 
{
	printf("\n Printing the result \n");
	for(size_t i=0; i < size; ++i)
	{
		for(size_t j=0; j < size; ++j)
		{
			printf("%d\t", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n End of printing the result \n");
}

int** alloc_matrix(size_t N)
{
	int** matrix= (int **) malloc(sizeof(int*) * N);
	for(int i = 0; i < N; ++i)
	{
		matrix[i] = (int*) malloc(sizeof(int) * N);
	}
	return matrix;
}

int** random_matrix(int** matrix, size_t N)
{
	const int RND = 100;
  	std::default_random_engine generator(time(NULL));
  	std::uniform_int_distribution<int> distribution(-RND, RND);

	for(int idx = 0; idx < N; ++idx)
	{
		for(int jdx = 0; jdx < N; ++jdx)
		{
			matrix[idx][jdx] = distribution(generator); 
		}
	}
	return matrix;
}

bool check(int** A, int** B, int** result, size_t N) 
{
	print_matrix(A, N);
	print_matrix(B, N);
	print_matrix(result, N);

	for(size_t i = 0; i < N; ++i)
	{
		for(size_t j = 0; j < N; ++j)
		{
			int value = 0;
			for(size_t k = 0; k < N; ++k)
			{
				result += A[i][k] * B[k][j];
			}

			if(value != result[i][j])
			{
				printf("\n Difference at position %lu, %lu\n",i,j);
				printf("\n reference %d, test %d\n", value, result[i][j]);
				return false;
			}
		}
	}
	return true;
}
#endif /* MATRIX_H_ */
