#include <mpi.h>
#include <random>

void generateArray(int* vector, size_t size)
{
	const int RND = 10;
  	std::default_random_engine generator(time(NULL));
  	std::uniform_int_distribution<int> distribution(1, RND);

	for (int i = 0; i < size; i ++)
	{
		vector[i] = distribution(generator);
	}
}

int test_size = 4;

int main( int argc, char **argv )
{
	MPI_Init (&argc, &argv); 

	int comm_size, myrank;
	int my2drank, mycoords[2]; 
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


	MPI_Request requests[4];
	MPI_Status status[4];

	MPI_Comm comm_2d;
	int dims[2], periods[2];
	dims[0] = dims[1] = sqrt(comm_size);
	periods[0] = periods[1] = 1; 
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);

	int sizelocal = test_size/dims[0];
	int a[sizelocal*sizelocal], b[sizelocal*sizelocal];
	generateArray(a, sizelocal*sizelocal);                              
	generateArray(b, sizelocal*sizelocal);

	MPI_Comm_rank(comm_2d, &my2drank); 
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords); 

	int uprank, downrank, leftrank, rightrank; 
	MPI_Cart_shift(comm_2d, 0, -1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_2d, 1, -1, &downrank, &uprank); 

	int *buffersA[2], *buffersB[2];

	buffersA[0] = a; 
	buffersA[1] = (int *)malloc(sizelocal*sizelocal*sizeof(int)); 
	buffersB[0] = b; 
	buffersB[1] = (int *)malloc(sizelocal*sizelocal*sizeof(int));

	int shiftsource, shiftdest; 
	MPI_Cart_shift(comm_2d, 0, -mycoords[0], &shiftsource, &shiftdest); 
	MPI_Sendrecv_replace(buffersA[0], sizelocal*sizelocal, MPI_INT,shiftdest, 1, shiftsource, 1, comm_2d, &status[0]); 

	MPI_Cart_shift(comm_2d, 1, -mycoords[1], &shiftsource, &shiftdest); 
	MPI_Sendrecv_replace(buffersB[0], sizelocal*sizelocal, MPI_INT,shiftdest, 1, shiftsource, 1, comm_2d, &status[0]);

	printf("Rank %d at point 1\n", my2drank);    

	for (int i = 0; i < 2; i ++)
	{
		printf("Rank %d started iteration %d\n", my2drank, i);
		MPI_Isend(buffersA[i%2], sizelocal*sizelocal, MPI_INT,leftrank, 1, comm_2d, &requests[0]); 
		MPI_Isend(buffersB[i%2], sizelocal*sizelocal, MPI_INT,uprank, 1, comm_2d, &requests[1]); 
		MPI_Irecv(buffersA[(i+1)%2], sizelocal*sizelocal, MPI_INT,rightrank, 1, comm_2d, &requests[2]); 
		MPI_Irecv(buffersB[(i+1)%2], sizelocal*sizelocal, MPI_INT,downrank, 1, comm_2d, &requests[3]);
		MPI_Waitall(4, requests, status);
		printf("Rank %d stopped waiting in iteration %d\n", my2drank, i);
		MPI_Barrier(comm_2d);
	}

	printf("Rank %d at point 2\n", my2drank);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
