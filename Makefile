all: build
	mpirun -np 8 ./mpi

build:
	mpicxx -O2 -std=c++11 mmverify.cpp -o mpi

clean:
	rm -rf mpi
