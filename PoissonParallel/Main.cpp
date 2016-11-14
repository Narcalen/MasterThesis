#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "..\MasterThesis\EquationDefinitions.h"

#define MASTER 0

using namespace std;

void multiply_matrix_by_vector(int matr_size, double** matrix, double* vector, double* result_vector);
void multiply_vector_by_scalar(int vector_size, double* vector, double scalar, double* result_vector);
void add_vectors(int vector_size, double* vector1, double* vector2, double* result_vector);
double scalar_product(int vector_size, double* vector1, double* vector2);
double norm(double vector_size, double* vector);

void main(int argc, char* argv[])
{
	const int coord_x_steps = (COORD_UPPER_BOUND_X - COORD_LOWER_BOUND_X) / delta_x + 1;
	const int coord_y_steps = (COORD_UPPER_BOUND_Y - COORD_LOWER_BOUND_Y) / delta_y + 1;
	const int matr_size = (coord_x_steps - 2) * (coord_y_steps - 2);

	int rank, world_size, error = 0;
	double** coeff = NULL;
	double* freeElements = NULL;
	double** u = NULL;

	double* result = NULL;
	double* z = NULL;
	double* oldR = NULL;
	double* Az = NULL;
	double* r = NULL;
	double alpha, beta, norm_free_elem, norm_r;

	double* tmp_vector = NULL;
	unsigned done = 0;

	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	MPI_Init(&argc, &argv);	/* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	cout << "Hello world from process " << rank << " of " << world_size << endl;

	if (matr_size % world_size != 0) {
		error = -1;
		MPI_Bcast(&error, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	}
	if (error != 0) {
		if (rank == MASTER) {
			cout << "Error: The number of processes is not matrix size divisor" << endl;
		}
		MPI_Finalize();
		exit(error);
	}


	int elements_per_processor = matr_size / world_size;

	if (rank == MASTER) {
		u = new double*[coord_x_steps];
		for (int i = 0; i < coord_x_steps; i++) {
			u[i] = new double[coord_y_steps];
			for (int j = 0; j < coord_y_steps; j++) {
				u[i][j] = -1;
			}
		}
		//set boundary values:
		for (int x = 0; x < coord_x_steps; x++) {
			u[x][0] = def::boundaryValue(COORD_LOWER_BOUND_X + x * delta_x, COORD_LOWER_BOUND_Y);
			u[x][coord_y_steps - 1] = def::boundaryValue(COORD_LOWER_BOUND_X + x * delta_x, COORD_LOWER_BOUND_Y + (coord_y_steps - 1) * delta_y);
		}
		for (int y = 0; y < coord_y_steps; y++) {
			u[0][y] = def::boundaryValue(COORD_LOWER_BOUND_X, COORD_LOWER_BOUND_Y + y * delta_y);
			u[coord_x_steps - 1][y] = def::boundaryValue(COORD_LOWER_BOUND_X + (coord_x_steps - 1) * delta_x, COORD_LOWER_BOUND_Y + y * delta_y);
		}

		coeff = new double*[matr_size];
		for (int i = 0; i < matr_size; i++) {
			coeff[i] = new double[matr_size];
			for (int j = 0; j < matr_size; j++) {
				coeff[i][j] = 0;
			}
			if (i >= coord_y_steps - 2) {
				coeff[i][i - (coord_y_steps - 2)] = 1 / (delta_x*delta_x);
			}
			if (i >= 1 && i % (coord_y_steps - 2) > 0) {
				coeff[i][i - 1] = 1 / (delta_y*delta_y);
			}

			coeff[i][i] = -2 / (delta_y*delta_y) - 2 / (delta_x*delta_x);

			if (i < matr_size - 1 && (i + 1) % (coord_y_steps - 2) > 0) {
				coeff[i][i + 1] = 1 / (delta_y*delta_y);
			}
			if (i < matr_size - (coord_y_steps - 2)) {
				coeff[i][i + (coord_y_steps - 2)] = 1 / (delta_x*delta_x);
			}

		}
		cout << "Coefficients: " << endl;
		for (int i = 0; i < matr_size; i++) {
			for (int j = 0; j < matr_size; j++) {
				cout << coeff[i][j] << "\t";
			}
			cout << endl;
		}

		freeElements = new double[matr_size];
		for (int i = 1; i < coord_x_steps - 1; i++) {
			for (int j = 1; j < coord_y_steps - 1; j++) {
				int coord = (i - 1)*(coord_y_steps - 2) + (j - 1);
				freeElements[coord] = -def::sigma(COORD_LOWER_BOUND_X + i * delta_x, COORD_LOWER_BOUND_Y + j * delta_y);
				if (j == 1) {
					freeElements[coord] += -u[i][j - 1] / (delta_y * delta_y);
				}
				else if (j == coord_y_steps - 2) {
					freeElements[coord] += -u[i][j + 1] / (delta_y * delta_y);
				}
				if (i == 1) {
					freeElements[coord] += -u[i - 1][j] / (delta_x * delta_x);
				}
				else if (i == coord_x_steps - 2) {
					freeElements[coord] += -u[i + 1][j] / (delta_x * delta_x);
				}
			}
		}
		cout << "Free elements: " << endl;
		for (int i = 0; i < matr_size; i++) {
			cout << freeElements[i] << "\t";
		}
		cout << endl;

		result = new double[matr_size];
		for (int i = 0; i < matr_size; i++) {
			result[i] = 1;
		}
		z = new double[matr_size];
		oldR = new double[matr_size];
		Az = new double[matr_size];
		r = new double[matr_size];
		tmp_vector = new double[matr_size];
	}


	//BEGIN: Conjugate gradients method

	//r = add(freeElements, multiply((T)-1, multiply(matrix, result)));
	multiply_matrix_by_vector(matr_size, coeff, result, tmp_vector);
	multiply_vector_by_scalar(matr_size, tmp_vector, -1, tmp_vector);
	add_vectors(matr_size, freeElements, tmp_vector, r);

	//z.assign(r.begin(), r.end());
	add_vectors(matr_size, freeElements, result, z);
	//if (rank == MASTER) {
	//	cout << "z: " << endl;
	//	for (int i = 0; i < matr_size; i++) {
	//		cout << z[i] << "\t";
	//	}
	//	cout << endl;
	//}

	//for (int i = 1; i <= matr_size; i++) {
	int i = 1;
	while (!done) {
		if (rank == MASTER) {
			cout << endl << "Iteration " << i++ << endl;
		}

		//Az = multiply(matrix, z);
		multiply_matrix_by_vector(matr_size, coeff, z, Az);
		//if (rank == MASTER) {
		//	cout << "Az: " << endl;
		//	for (int i = 0; i < matr_size; i++) {
		//		cout << Az[i] << "\t";
		//	}
		//	cout << endl;
		//}

		//alpha = scalar(r, r) / scalar(Az, z);
		alpha = scalar_product(matr_size, z, z) / scalar_product(matr_size, Az, z);
		//if (rank == MASTER) {
		//	cout << "alpha = " << alpha << endl;
		//}

		//result = add(result, multiply(alpha, z));
		multiply_vector_by_scalar(matr_size, z, alpha, tmp_vector);
		add_vectors(matr_size, result, tmp_vector, result);

		//oldR.assign(r.begin(), r.end());
		if (rank == MASTER) {
			for (int i = 0; i < matr_size; i++) {
				//cout << result[i] << "\t";
				oldR[i] = r[i];
			}
			//cout << endl;
		}

		//r = add(r, multiply(-1 * alpha, Az));
		multiply_vector_by_scalar(matr_size, Az, -alpha, tmp_vector);
		add_vectors(matr_size, r, tmp_vector, r);

		//beta = scalar(r, r) / scalar(oldR, oldR);
		beta = scalar_product(matr_size, r, r) / scalar_product(matr_size, oldR, oldR);
		//if (rank == MASTER) {
		//	cout << "beta = " << beta << endl;
		//}

		//z = add(r, multiply(beta, z));
		multiply_vector_by_scalar(matr_size, z, beta, tmp_vector);
		add_vectors(matr_size, r, tmp_vector, z);
		if (rank == MASTER && norm(matr_size, r) / norm(matr_size, freeElements) < EPSILON) {
			done = 1;
		}
		MPI_Bcast(&done, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	}
	//END: Conjugate gradients method

	if (rank == MASTER) {
		int j = 0;
		for (int x = 1; x < coord_x_steps - 1; x++) {
			for (int y = 1; y < coord_y_steps - 1; y++) {
				u[x][y] = result[j++];
			}
		}
		cout << "Result: " << endl;
		//util::print_matrix(u);

		cout << "\t";
		for (int j = 0; j < coord_y_steps; j++) {
			cout << COORD_LOWER_BOUND_Y + j * delta_y << "\t";
		}
		cout << endl;

		for (int i = 0; i < coord_x_steps; i++) {
			cout << COORD_LOWER_BOUND_X + i * delta_x << "\t";
			for (int j = 0; j < coord_y_steps; j++) {
				cout << u[i][j] << "\t";
			}
			cout << endl;
		}
		//free global data
		delete[] freeElements;
		for (int i = 0; i < coord_x_steps; i++) {
			delete[] u[i];
		}
		delete[] u;
		for (int i = 0; i < matr_size; i++) {
			delete[] coeff[i];
		}
		delete[] coeff;

		cout << "The end!" << endl;
		getchar();
	}

	MPI_Finalize();
}

void multiply_matrix_by_vector(int matr_size, double** matrix, double* vector, double* result_vector) {
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	int elements_per_processor = matr_size / world_size;

	double* local_matr_row = new double[elements_per_processor * matr_size];
	double* local_vector = new double[matr_size];
	double* local_result = new double[elements_per_processor];

	for (int i = 0; i < elements_per_processor; i++) {
		local_result[i] = 0;
	}
	double* matrix_to_send = NULL;

	if (rank == MASTER) {
		local_vector = vector;
		matrix_to_send = new double[matr_size*matr_size];
		for (int i = 0; i < matr_size; i++) {
			for (int j = 0; j < matr_size; j++) {
				matrix_to_send[i*matr_size + j] = matrix[i][j];
			}
		}
	}
	MPI_Bcast(local_vector, matr_size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Scatter(matrix_to_send, matr_size * elements_per_processor, MPI_DOUBLE, local_matr_row, matr_size * elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	for (int i = 0; i < elements_per_processor; i++) {
		for (int j = 0; j < matr_size; j++) {
			local_result[i] += local_matr_row[i * matr_size + j] * local_vector[j];
		}
	}

	MPI_Gather(local_result, elements_per_processor, MPI_DOUBLE, result_vector, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	//free local data
	delete[] local_result;
	delete[] local_matr_row;
}

void multiply_vector_by_scalar(int vector_size, double* vector, double scalar, double* result_vector) {
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	int elements_per_processor = vector_size / world_size;

	double* local_vector = new double[elements_per_processor];

	MPI_Bcast(&scalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Scatter(vector, elements_per_processor, MPI_DOUBLE, local_vector, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	for (int i = 0; i < elements_per_processor; i++) {
		local_vector[i] *= scalar;
	}
	MPI_Gather(local_vector, elements_per_processor, MPI_DOUBLE, result_vector, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	delete[] local_vector;
}

void add_vectors(int vector_size, double* vector1, double* vector2, double* result_vector) {
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	int elements_per_processor = vector_size / world_size;

	double* local_vector1 = new double[elements_per_processor];
	double* local_vector2 = new double[elements_per_processor];

	MPI_Scatter(vector1, elements_per_processor, MPI_DOUBLE, local_vector1, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Scatter(vector2, elements_per_processor, MPI_DOUBLE, local_vector2, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	for (int i = 0; i < elements_per_processor; i++) {
		local_vector1[i] += local_vector2[i];
	}
	MPI_Gather(local_vector1, elements_per_processor, MPI_DOUBLE, result_vector, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	delete[] local_vector1;
	delete[] local_vector2;
}

double scalar_product(int vector_size, double* vector1, double* vector2) {
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	int elements_per_processor = vector_size / world_size;

	double* local_vector1 = new double[elements_per_processor];
	double* local_vector2 = new double[elements_per_processor];
	double local_sum = 0;
	double* local_sums = NULL;
	if (rank == MASTER) {
		local_sums = new double[world_size];
	}

	MPI_Scatter(vector1, elements_per_processor, MPI_DOUBLE, local_vector1, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Scatter(vector2, elements_per_processor, MPI_DOUBLE, local_vector2, elements_per_processor, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	for (int i = 0; i < elements_per_processor; i++) {
		local_sum += local_vector1[i] * local_vector2[i];
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&local_sum, 1, MPI_DOUBLE, local_sums, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

	delete[] local_vector1;
	delete[] local_vector2;

	if (rank == MASTER) {
		double result = 0;
		for (int i = 0; i < world_size; i++) {
			result += local_sums[i];
		}

		delete[] local_sums;

		return result;
	}

	//int received_num = 0;
	//for (int k = log(world_size) / log(2) - 1; k >= 0; k--) {
	//	if (rank < pow(2, k)) {
	//		MPI_Recv(&received_num, 1, MPI_DOUBLE, rank + pow(2, k), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//		cout << "Processor " << rank << " received number " << received_num << " from processor " << rank + pow(2, k) << endl;
	//		local_sum += received_num;
	//		cout << "Processor " << rank << " has calculated subsum " << local_sum << endl;
	//	}
	//	else if (rank < pow(2, k + 1))
	//	{
	//		MPI_Send(&local_sum, 1, MPI_DOUBLE, rank - pow(2, k), 0, MPI_COMM_WORLD);
	//		cout << "Processor " << rank << " sent number " << local_sum << " to processor " << rank - pow(2, k) << endl;
	//	}
	//}

	//if (rank == MASTER) {
	//	cout << "Total sum calculated in the function: " << local_sum << endl;
	//	result = &local_sum;
	//}
}

double norm(double vector_size, double* vector)
{
	double norm = 0;
	for (int i = 0; i < vector_size; i++) {
		norm += vector[i] * vector[i];
	}

	return sqrt(norm);
}