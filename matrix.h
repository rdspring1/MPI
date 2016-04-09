#ifndef MATRIX_H_
#define MATRIX_H_

#include <ctime>
#include <random>
#include <algorithm>

#define OFFSET(idx, N) ((idx)*(N))
#define IDX(i,j, N) ((i)*(N)+(j))

void print_matrix(int* matrix, const size_t N, const size_t gridDim) 
{
	if(!matrix) { return; }

	const size_t blockDim = N / gridDim;
	printf("\n Printing the result \n");
	for(size_t i=0; i < N; ++i)
	{
		for(size_t j=0; j < N; ++j)
		{
			int gX = i / blockDim;
			int gY = j / blockDim;
			int idx = i % blockDim;
			int jdx = j % blockDim;
			printf("%d\t", matrix[OFFSET(IDX(gX, gY, gridDim), blockDim*blockDim) + IDX(idx, jdx, blockDim)]);
		}
		printf("\n");
	}
	printf("\n End of printing the result \n");
}

int* zero_matrix(size_t N)
{
	return new int[N*N] {};
}

int* count_matrix(const size_t N, const size_t gridDim)
{
	int* matrix = zero_matrix(N);
	const size_t blockDim = N / gridDim;
	int count = 0;
	for(size_t gX = 0; gX < gridDim; ++gX)
	{
		for(size_t idx = 0; idx < blockDim; ++idx)
		{
			for(size_t gY = 0; gY < gridDim; ++gY)
			{
				for(size_t jdx = 0; jdx < blockDim; ++jdx)
				{
					matrix[OFFSET(IDX(gX, gY, gridDim), blockDim*blockDim) + IDX(idx, jdx, blockDim)] = (++count);
				}
			}
		}
	}
	return matrix;
}

int* random_matrix(size_t N)
{
	const int RND = 100;
	std::default_random_engine generator(time(NULL));
	std::uniform_int_distribution<int> distribution(-RND, RND);

	int* matrix = zero_matrix(N);
	for(int idx = 0; idx < N*N; ++idx)
	{
		matrix[idx] = distribution(generator); 
	}
	return matrix;
}

bool check(int* A, int* B, int* result, size_t N, size_t gridDim) 
{
	print_matrix(A, N, gridDim);
	print_matrix(B, N, gridDim);
	print_matrix(result, N, gridDim);

	const size_t blockDim = N / gridDim;
	for(size_t i = 0; i < N; ++i)
	{
		int gI = i / blockDim;
		int idx = i % blockDim;

		for(size_t j = 0; j < N; ++j)
		{
			int gJ = j / blockDim;
			int jdx = j % blockDim;

			int index = OFFSET(IDX(gI, gJ, gridDim), blockDim*blockDim) + IDX(idx, jdx, blockDim);
			int value = 0;
			for(size_t k = 0; k < N; ++k)
			{
				int gK = k / blockDim;
				int kdx = k % blockDim;

				int A_index = OFFSET(IDX(gI, gK, gridDim), blockDim*blockDim) + IDX(idx, kdx, blockDim);
				int B_index = OFFSET(IDX(gK, gJ, gridDim), blockDim*blockDim) + IDX(kdx, jdx, blockDim);
				value += A[A_index] * B[B_index];
			}

			if(value != result[index])
			{
				printf("\n Difference at position %lu, %lu\n",i,j);
				printf("\n reference %d, test %d\n", value, result[index]);
				return false;
			}
		}
	}
	return true;
}
#endif /* MATRIX_H_ */
