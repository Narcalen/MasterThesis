#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "..\MasterThesis\EquationDefinitions.h"
#include "..\FiniteVolumesConsecutive\Element.h"
#include "..\FiniteVolumesConsecutive\Node.h"
#include "..\FiniteVolumesConsecutive\FVMUtil.h"

#define MASTER 0

using namespace std;

void multiply_matrix_by_vector(int matr_size, double** matrix, double* vector, double* result_vector);
void multiply_vector_by_scalar(int vector_size, double* vector, double scalar, double* result_vector);
void add_vectors(int vector_size, double* vector1, double* vector2, double* result_vector);
double scalar_product(int vector_size, double* vector1, double* vector2);
double norm(double vector_size, double* vector);

void main(int argc, char* argv[])
{
	int rank, world_size, error = 0;
	double** coeff = NULL;
	double* freeElements = NULL;

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

	int a, b, c;
	ifstream infile;

	vector<Element> elements;
	vector<Node> nodes;

	//read nodes data
	infile.open("distmesh_nodes.txt", ios::in);
	if (!infile) {
		cout << "Cannot open input file distmesh_nodes.txt" << endl;
	}

	Node node(2);
	while (infile >> a >> b)
	{
		node.set(0, a);
		node.set(1, b);
		nodes.push_back(node);
	}
	infile.close();
	//cout << "Nodes: " << endl;
	//for (Node node: nodes) {
	//	cout << node.get(0) << " " << node.get(1) << endl;
	//}

	//read elements data
	infile.open("distmesh_elements.txt", ios::in);
	if (!infile) {
		cout << "Cannot open input file distmesh_elements.txt" << endl;
	}

	Element element(3);
	while (infile >> a >> b >> c)
	{
		element.setNode(0, a);
		element.setNode(1, b);
		element.setNode(2, c);
		elements.push_back(element);
	}
	infile.close();

	const int matr_size = nodes.size();
	cout << "Number of elements: " << matr_size << endl;


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
		freeElements = new double[matr_size];
		coeff = new double*[matr_size];
		for (int i = 0; i < matr_size; i++) {
			coeff[i] = new double[matr_size];
			freeElements[i] = 0;
			for (int j = 0; j < matr_size; j++) {
				coeff[i][j] = 0;
			}
		}

		for (Element element : elements) {
			int node1Num = element.getNode(0);
			Node node1 = nodes[node1Num];
			int node2Num = element.getNode(1);
			Node node2 = nodes[node2Num];
			int node3Num = element.getNode(2);
			Node node3 = nodes[node3Num];
			double multiplier = 1 / (2 * abs(detD(node1, node2, node3)));

			//add local main diagonal to global main diagonal
			double val = multiplier * (deltaX(node2, node3)*deltaX(node2, node3) + deltaY(node2, node3) * deltaY(node2, node3));
			coeff[node1Num][node1Num] += val;

			val = multiplier * (deltaX(node3, node1)*deltaX(node3, node1) + deltaY(node3, node1) * deltaY(node3, node1));
			coeff[node2Num][node2Num] += val;

			val = multiplier * (deltaX(node1, node2)*deltaX(node1, node2) + deltaY(node1, node2) * deltaY(node1, node2));
			coeff[node3Num][node3Num] += val;

			//add local elements above and below main diagonal to global diagonals
			val = multiplier * (deltaX(node3, node1)*deltaX(node2, node3) + deltaY(node3, node1) * deltaY(node2, node3));
			coeff[node1Num][node2Num] += val;
			coeff[node2Num][node1Num] += val;

			val = multiplier * (deltaX(node1, node2)*deltaX(node2, node3) + deltaY(node1, node2) * deltaY(node2, node3));
			coeff[node1Num][node3Num] += val;
			coeff[node3Num][node1Num] += val;

			val = multiplier * (deltaX(node1, node2)*deltaX(node3, node1) + deltaY(node1, node2) * deltaY(node3, node1));
			coeff[node3Num][node2Num] += val;
			coeff[node2Num][node3Num] += val;

			//set values for the right side of the equations system
			double f1 = -def::sigma(node1.get(0), node1.get(1));
			double f2 = -def::sigma(node2.get(0), node2.get(1));
			double f3 = -def::sigma(node3.get(0), node3.get(1));
			multiplier = abs(detD(node1, node2, node3)) / 216;

			freeElements[node1Num] += multiplier * (22 * f1 + 7 * f2 + 7 * f3);
			freeElements[node2Num] += multiplier * (7 * f1 + 22 * f2 + 7 * f3);
			freeElements[node3Num] += multiplier * (7 * f1 + 7 * f2 + 22 * f3);
		}

		//boundary conditions 1st type
		for (unsigned i = 0; i < nodes.size(); i++) {
			Node node = nodes[i];
			if (node.get(0) == COORD_LOWER_BOUND_X || node.get(0) == COORD_UPPER_BOUND_X
				|| node.get(1) == COORD_LOWER_BOUND_Y || node.get(1) == COORD_UPPER_BOUND_Y) {
				for (unsigned j = 0; j < matr_size; j++) {
					coeff[i][j] = 0;
					coeff[j][i] = 0;
				}
				coeff[i][i] = 1;
				freeElements[i] = def::boundaryValue(node.get(0), node.get(1));
			}
		}

		cout << "Coefficients: " << endl;
		for (int i = 0; i < matr_size; i++) {
			for (int j = 0; j < matr_size; j++) {
				cout << coeff[i][j] << "\t";
			}
			cout << endl;
		}

		cout << "Free elements: " << endl;
		for (int i = 0; i < matr_size; i++) {
			cout << freeElements[i] << "\t";
		}
		cout << endl;

		result = new double[matr_size];
		for (int i = 0; i < matr_size; i++) {
			result[i] = 0;
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
		i++;
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
				cout << result[i] << "\t";
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
		if (rank == MASTER && norm(matr_size, r) / norm(matr_size, freeElements) < EPSILON || i >= matr_size  /** matr_size  * matr_size*/) {
			done = 1;
		}
		MPI_Bcast(&done, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	}
	//END: Conjugate gradients method

	if (rank == MASTER) {

		cout << "Result: " << endl;
		cout << "x\ty\tvalue" << endl;
		for (unsigned i = 0; i < nodes.size(); i++) {
			Node node = nodes[i];
			cout << node.get(0) << "\t" << node.get(1) << "\t" << result[i] << endl;
		}

		//free global data
		delete[] freeElements;
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
	delete[] matrix_to_send;
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