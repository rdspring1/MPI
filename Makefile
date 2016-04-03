all:
	g++ -std=c++11 mmverify.cpp -o mpi

clean:
	rm -rf mpi
