#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#define MASTER 0

using namespace std;

void multiply_matrix_by_vector(int matr_size, double** matrix, double* vector, double* result_vector);
void multiply_vector_by_scalar(int vector_size, double* vector, double scalar, double* result_vector);
void add_vectors(int vector_size, double* vector1, double* vector2, double* result_vector);
double scalar_product(int vector_size, double* vector1, double* vector2);
double norm(double vector_size, double* vector);

void main(int argc, char* argv[])
{
	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	MPI_Init(&argc, &argv);	/* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	cout << "Hello world from process " << rank << " of " << world_size << endl;

	MPI_Finalize();
}

