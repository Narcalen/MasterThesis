#include <vector>
#include <iostream>
#include "SimpleIteration.h"
#include "math.h"
#include "Util.h"

using namespace std;


template<typename T> vector<T>* JacobiSolver(SparseMatrix<T> matrix, vector<T>* freeElements){
	int N = freeElements -> size();
	vector<T>* result = new vector<T>(N, 0.5);
	vector<T>* temp = new vector<T>(N);
	double norm;

	int j = 0;
	do {
		for (int i = 0; i < N; i++){
			temp->at(i) = freeElements->at(i);
			for (int j = 0; j < N; j++){
				if (i != j){
					temp->at(i) -= matrix(i,j) * result->at(j);
				}
			}
			temp->at(i) /= matrix(i,i);
		}
		norm = distance(result, temp);
		result -> assign(temp -> begin(), temp -> end());
		cout << "iteration " << j++ << endl;
		print(cout, result);
	}
	while (norm > EPSILON);

	return result;
}

template vector<float>* JacobiSolver(SparseMatrix<float> matrix, vector<float>* freeElements);
template vector<double>* JacobiSolver(SparseMatrix<double> matrix, vector<double>* freeElements);