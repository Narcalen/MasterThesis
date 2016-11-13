#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "..\MasterThesis\EquationDefinitions.h"

#define MASTER 0

using namespace std;

void multiply_matrix_by_vector(int matr_size, double** matrix, int vector_size, double* vector, double* result_vector);

void main(int argc, char* argv[])
{
	const int coord_x_steps = (COORD_UPPER_BOUND_X - COORD_LOWER_BOUND_X) / delta_x + 1;
	const int coord_y_steps = (COORD_UPPER_BOUND_Y - COORD_LOWER_BOUND_Y) / delta_y + 1;
	const int matr_size = (coord_x_steps - 2) * (coord_y_steps - 2);

	int rank, world_size;
	double** coeff = NULL;
	double* coeff_to_send = NULL;
	double* freeElements = new double[matr_size];
	double** u = NULL;
	double* result_vector = NULL;
	//double* local_matr_row = NULL;
	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	MPI_Init(&argc, &argv);	/* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	cout << "Hello world from process " << rank << " of " << world_size << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == MASTER) {
		if (matr_size % world_size != 0) {
			MPI_Finalize();
			exit(-1);
		}
	}
	int elements_per_processor = matr_size / world_size;

	//local_matr_row = new double[elements_per_processor * matr_size];
	////for (int i = 0; i < elements_per_processor; i++) {
	////	local_matr_row[i] = new double[matr_size];
	////}
	//double* local_vector = new double[matr_size];
	//double* gathered_result = new double[matr_size];
	//double* local_result = new double[elements_per_processor];
	//for (int i = 0; i < elements_per_processor; i++) {
	//	local_result[i] = 0;
	//}

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
		//cout << "Function: " << endl;
		//for (int i = 0; i < coord_x_steps; i++) {
		//	for (int j = 0; j < coord_y_steps; j++) {
		//		cout << u[i][j] << "\t";
		//	}
		//	cout << endl;
		//}

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

		//coeff_to_send = new double[matr_size*matr_size];
		//for (int i = 0; i < matr_size; i++) {
		//	for (int j = 0; j < matr_size; j++) {
		//		coeff_to_send[i*matr_size + j] = coeff[i][j];
		//	}
		//}

		//cout << "Coefficients: " << endl;
		//for (int i = 0; i < matr_size; i++) {
		//	for (int j = 0; j < matr_size; j++) {
		//		cout << coeff_to_send[i*matr_size + j] << "\t";
		//	}
		//	cout << endl;
		//}

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
		//cout << "Free elements: " << endl;
		//for (int j = 0; j < matr_size; j++) {
		//	cout << freeElements[j] << "\t";
		//}
		//cout << endl;

		result_vector = new double[matr_size];
	}


	//MPI_Bcast(freeElements, matr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Scatter(coeff_to_send, matr_size * elements_per_processor, MPI_DOUBLE, local_matr_row, matr_size * elements_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	//for (int i = 0; i < world_size; i++) {
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	if (rank == i) {
	//		cout << "Processor " << i << " received data:" << endl;
	//		cout << "Matrix row(s): " << endl;
	//		for (int row = 0; row < elements_per_processor; row++) {
	//			for (int col = 0; col < matr_size; col++) {
	//				cout << local_matr_row[row * matr_size + col] << "\t";
	//			}
	//			cout << endl;
	//		}
	//		cout << "Vector: " << endl;
	//		for (int j = 0; j < matr_size; j++) {
	//			cout << freeElements[j] << "\t";
	//		}
	//		cout << endl;
	//	}
	//}

	//for (int i = 0; i < elements_per_processor; i++) {
	//	for (int j = 0; j < matr_size; j++) {
	//		local_result[i] += local_matr_row[i * matr_size + j] * freeElements[j];
	//	}
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//cout << "processor " << rank << " - result: " << local_result[0] << endl;

	//MPI_Gather(local_result, elements_per_processor, MPI_DOUBLE, gathered_result, elements_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//if (rank == MASTER) {
	//	double res = 0;
	//	cout << "gathered vector:" << endl;
	//	for (int i = 0; i < matr_size; i++) {
	//		cout << gathered_result[i] << "\t";
	//		res += gathered_result[i];
	//	}
	//	cout << endl << res << endl;

	//	cout << "expected result:" << endl;
	//	for (int i = 0; i < matr_size; i++) {
	//		double res = 0;
	//		for (int j = 0; j < matr_size; j++) {
	//			res += coeff[i][j] * freeElements[j];
	//		}
	//		cout << res << "\t";
	//	}
	//	cout << endl;
	//}

	//delete[] local_vector;
	//delete[] gathered_result;
	//delete[] local_result;
	//delete[] local_matr_row;

	multiply_matrix_by_vector(matr_size, coeff, matr_size, freeElements, result_vector);

	if (rank == MASTER) {
		cout << "Received vector:" << endl;
		for (int i = 0; i < matr_size; i++) {
			cout << result_vector[i] << "\t";
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
		delete[] coeff_to_send;

		cout << "The end!" << endl;
		getchar();
	}

	MPI_Finalize();
}

void multiply_matrix_by_vector(int matr_size, double** matrix, int vector_size, double* vector, double* result_vector) {
	int rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	/* get number of processes */
	int elements_per_processor = matr_size / world_size;

	cout << "multiply function called" << endl;
	cout << "Hello world from process " << rank << " of " << world_size << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	double* local_matr_row = new double[elements_per_processor * matr_size];
	double* local_vector = new double[matr_size];
	double* local_result = new double[elements_per_processor];
	for (int i = 0; i < elements_per_processor; i++) {
		local_result[i] = 0;
	}
	double* matrix_to_send = NULL;
	cout << "Variables initialized" << endl;

	if (rank == MASTER) {
		matrix_to_send = new double[matr_size*matr_size];
		for (int i = 0; i < matr_size; i++) {
			for (int j = 0; j < matr_size; j++) {
				matrix_to_send[i*matr_size + j] = matrix[i][j];
			}
		}
	}
	cout << "Master process - Broadcasting data" << endl;
	MPI_Bcast(vector, matr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix_to_send, matr_size * elements_per_processor, MPI_DOUBLE, local_matr_row, matr_size * elements_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	for (int i = 0; i < world_size; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == i) {
			cout << "Processor " << i << " received data:" << endl;
			cout << "Matrix row(s): " << endl;
			for (int row = 0; row < elements_per_processor; row++) {
				for (int col = 0; col < matr_size; col++) {
					cout << local_matr_row[row * matr_size + col] << "\t";
				}
				cout << endl;
			}
			cout << "Vector: " << endl;
			for (int j = 0; j < matr_size; j++) {
				cout << vector[j] << "\t";
			}
			cout << endl;
		}
	}

	for (int i = 0; i < elements_per_processor; i++) {
		for (int j = 0; j < matr_size; j++) {
			local_result[i] += local_matr_row[i * matr_size + j] * vector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "processor " << rank << " - result: " << local_result[0] << endl;

	MPI_Gather(local_result, elements_per_processor, MPI_DOUBLE, result_vector, elements_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == MASTER) {
		cout << "gathered vector:" << endl;
		for (int i = 0; i < matr_size; i++) {
			cout << result_vector[i] << "\t";
		}
		cout << endl;
	}

	//free local data
	delete[] local_vector;
	delete[] local_result;
	delete[] local_matr_row;
}
