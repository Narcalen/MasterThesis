#include <vector>
#include <iostream>
#include "SimpleIteration.h"
#include "math.h"
#include "Util.h"

using namespace std;


template<typename T> vector<T> JacobiSolver(SparseMatrix<T> matrix, vector<T> freeElements){
	int N = freeElements.size();
	vector<T> result(N, 0.5);
	vector<T> temp(N);
	double norm;

	int j = 0;
	do {
		for (int i = 0; i < N; i++){
			temp[i] = freeElements[i];
			for (int j = 0; j < N; j++){
				if (i != j){
					temp[i] -= matrix(i,j) * result[j];
				}
			}
			temp[i] /= matrix(i,i);
		}
		norm = util::distance(result, temp);
		result.assign(temp.begin(), temp.end());
		cout << "iteration " << j++ << endl;
		util::print(result);
	}
	while (norm > EPSILON);

	return result;
}

template vector<float> JacobiSolver(SparseMatrix<float> matrix, vector<float> freeElements);
template vector<double> JacobiSolver(SparseMatrix<double> matrix, vector<double> freeElements);