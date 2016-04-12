#ifndef MPI_MATRIX_
#define MPI_MATRIX_

#include <mpi.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"

void matrix_multiply(int* A, int* B, int*C, const int N)
{
	for(int idx = 0; idx < N; ++idx)
	{
		for(int jdx = 0; jdx < N; ++jdx)
		{
			for(int kdx = 0; kdx < N; ++kdx)
			{
				C[IDX(idx, jdx, N)] += A[IDX(idx, kdx, N)] * B[IDX(kdx, jdx, N)];
			}
		}
	}
}

void MPI_Check(int error_code)
{
	assert(error_code == MPI_SUCCESS);
}

void MPI_Group_delete(MPI_Group* grp)
{
	if(*grp != MPI_GROUP_NULL)
	{
		MPI_Group_free(grp);
	}
}

void MPI_Comm_delete(MPI_Comm* comm)
{
	if(*comm != MPI_COMM_NULL)
	{
		MPI_Comm_free(comm);
	}
}

int* mpi_matrix_multiplication(int* A, int* B, const int MATRIX_SIZE, const int depth)
{
	// Global MPI
	MPI_Group grp_world;
	MPI_Comm_group(MPI_COMM_WORLD, &grp_world);

	int num_procs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	const int ndims = 3;
	const int gridDim = std::sqrt(num_procs / depth);
	const int blockDim = MATRIX_SIZE / gridDim;
	const int rounds = std::sqrt(num_procs / std::pow(depth, 3));

	// 3D communication MPI
	MPI_Comm comm_3d;
	int dims[ndims], periods[ndims];
	dims[0] = dims[1] = gridDim;
	dims[2] = depth;
	periods[0] = periods[1] = periods[2] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, &comm_3d);

	int my_3d_rank, my_coords[ndims]; 
	MPI_Comm_rank(comm_3d, &my_3d_rank); 
	MPI_Cart_coords(comm_3d, my_3d_rank, ndims, my_coords); 

	// Allocate Sub-Matrices
	int* Asub = zero_matrix(blockDim);
	int* Bsub = zero_matrix(blockDim);
	int* Csub = zero_matrix(blockDim);

	// Initialize front face communication topology 
	int* front_face = new int[gridDim*gridDim];
	for(int idx = 0; idx < gridDim; ++idx)
	{
		for(int jdx = 0; jdx < gridDim; ++jdx)
		{
			MPI_Cart_rank(comm_3d, new int[ndims] {idx, jdx, 0}, &front_face[IDX(idx, jdx, gridDim)]);
		}
	}

	MPI_Group grp_front_face;
	MPI_Group_incl(grp_world, gridDim*gridDim, front_face, &grp_front_face);
	MPI_Comm comm_front_face;
	MPI_Comm_create(MPI_COMM_WORLD, grp_front_face, &comm_front_face);
	delete[] front_face;

	// Scatter Sub-Matrices across front face
	if(comm_front_face != MPI_COMM_NULL)
	{
		MPI_Scatter(A, blockDim*blockDim, MPI_INT, Asub, blockDim*blockDim, MPI_INT, 0, comm_front_face);
		MPI_Scatter(B, blockDim*blockDim, MPI_INT, Bsub, blockDim*blockDim, MPI_INT, 0, comm_front_face);
	}

	// Initialize depth commication topology
	std::vector<MPI_Group> grp_depth(gridDim*gridDim);
	std::vector<MPI_Comm> comm_depth(gridDim*gridDim);
	int* depth_ranks = new int[depth];
	for(int idx = 0; idx < gridDim; ++idx)
	{
		for(int jdx = 0; jdx < gridDim; ++jdx)
		{
			for(int kdx = 0; kdx < depth; ++kdx)
			{
				MPI_Cart_rank(comm_3d, new int[ndims] {idx, jdx, kdx}, &depth_ranks[kdx]);
			}
			const int index = IDX(idx, jdx, gridDim);
			MPI_Group_incl(grp_world, depth, depth_ranks, &grp_depth[index]);
			MPI_Comm_create(MPI_COMM_WORLD, grp_depth[index], &comm_depth[index]);

			// Broadcast Sub-Matrices along depth
			if(comm_depth[index] != MPI_COMM_NULL)
			{
				MPI_Bcast(Asub, blockDim*blockDim, MPI_INT, 0, comm_depth[index]);
				MPI_Bcast(Bsub, blockDim*blockDim, MPI_INT, 0, comm_depth[index]);
			}
		}
	}
	delete[] depth_ranks;

	// Initial Circular Shift
	int r = (my_coords[1] + my_coords[0] - my_coords[2] * rounds) % gridDim;
	int s = (my_coords[1] - my_coords[0] + my_coords[2] * rounds) % gridDim;
	int s1 = (my_coords[0] - my_coords[1] + my_coords[2] * rounds) % gridDim;

	int Asend, Adest, Bsend, Bdest;
	MPI_Cart_rank(comm_3d, new int[ndims] {my_coords[0], s, my_coords[2]}, &Asend);
	MPI_Cart_rank(comm_3d, new int[ndims] {my_coords[0], r, my_coords[2]}, &Adest);
	MPI_Cart_rank(comm_3d, new int[ndims] {s1, my_coords[1], my_coords[2]}, &Bsend);
	MPI_Cart_rank(comm_3d, new int[ndims] {r, my_coords[1], my_coords[2]}, &Bdest);

	MPI_Status status[2];
	MPI_Sendrecv_replace(Asub, blockDim*blockDim, MPI_INT, Asend, 1, Adest, 1, comm_3d, &status[0]); 
	MPI_Sendrecv_replace(Bsub, blockDim*blockDim, MPI_INT, Bsend, 1, Bdest, 1, comm_3d, &status[1]);

	//printf("Rank %d at point 1\n", my_3d_rank);    
	matrix_multiply(Asub, Bsub, Csub, blockDim);

	// Cannon's Algorithm
	MPI_Request requests[4];
	int uprank, downrank, leftrank, rightrank; 
        MPI_Cart_shift(comm_3d, 0, 1, &uprank, &downrank);
        MPI_Cart_shift(comm_3d, 1, 1, &leftrank, &rightrank);
	//printf("x:%d y:%d z:%d rank:%d Ashift:%d Bshift:%d\n", my_coords[0], my_coords[1], my_coords[2], my_3d_rank, leftrank, uprank);

	int *buffersA[2]; 
	int *buffersB[2];
	buffersA[0] = Asub; 
	buffersA[1] = zero_matrix(blockDim); 
	buffersB[0] = Bsub; 
	buffersB[1] = zero_matrix(blockDim);

	for (int idx = 0; idx < rounds-1; ++idx)
	{
		//printf("Rank %d started iteration %lu\n", my_3d_rank, idx);
		//MPI_Sendrecv_replace(Asub, blockDim*blockDim, MPI_INT, rightrank, 1, leftrank, 1, comm_3d, &status[0]); 
		//MPI_Sendrecv_replace(Bsub, blockDim*blockDim, MPI_INT, downrank, 1, uprank, 1, comm_3d, &status[1]);

		MPI_Isend(buffersA[idx % 2],     blockDim*blockDim, MPI_INT, rightrank, 1, comm_3d, &requests[0]); 
		MPI_Isend(buffersB[idx % 2],     blockDim*blockDim, MPI_INT, downrank,  1, comm_3d, &requests[1]); 
		MPI_Irecv(buffersA[(idx+1) % 2], blockDim*blockDim, MPI_INT, leftrank,  1, comm_3d, &requests[2]); 
		MPI_Irecv(buffersB[(idx+1) % 2], blockDim*blockDim, MPI_INT, uprank,    1, comm_3d, &requests[3]);
		MPI_Waitall(4, requests, status);
		//printf("Rank %d stopped waiting in iteration %lu\n", my_3d_rank, idx);

		matrix_multiply(Asub, Bsub, Csub, blockDim);

		MPI_Barrier(comm_3d);
	}
	//printf("Rank %d at point 2\n", my_3d_rank);

	// Sum-Reduction along depth
	int* result_sub = zero_matrix(blockDim);
	for(int idx = 0; idx < gridDim; ++idx)
	{
		for(int jdx = 0; jdx < gridDim; ++jdx)
		{
			const int index = IDX(idx, jdx, gridDim);
			if(comm_depth[index] != MPI_COMM_NULL)
			{
				MPI_Reduce(Csub, result_sub, blockDim*blockDim, MPI_INT, MPI_SUM, 0, comm_depth[index]);
			}
		}
	}
	delete[] Asub;
	delete[] Bsub;
	delete[] Csub;

	// Delete MPI Depth
	for(size_t idx = 0; idx < comm_depth.size(); ++idx)
	{
		MPI_Comm_delete(&comm_depth[idx]);
	}

	for(size_t idx = 0; idx < grp_depth.size(); ++idx)
	{
		MPI_Group_delete(&grp_depth[idx]);
	}

	int* result = NULL;
	if(myrank == 0)
	{
		result = zero_matrix(MATRIX_SIZE);
	}
	// Gather Sub-Matrices along front face
	if(comm_front_face != MPI_COMM_NULL)
	{
		MPI_Gather(result_sub, blockDim*blockDim, MPI_INT, result, blockDim*blockDim, MPI_INT, 0, comm_front_face);
	}
	delete[] result_sub;

	// Delete MPI Front Face 
	MPI_Comm_delete(&comm_front_face);
	MPI_Group_delete(&grp_front_face);

	MPI_Comm_free(&comm_3d);
	MPI_Barrier(MPI_COMM_WORLD);

	return result;
}

#endif /* MPI_MATRIX_ */
