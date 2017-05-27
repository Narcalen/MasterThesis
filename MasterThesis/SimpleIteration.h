#pragma once
#include "..\Sparse-Matrix\src\SparseMatrix\SparseMatrix.h"
#include "EquationDefinitions.h"
#include <vector>
#include <assert.h>
//#include <iostream>
//#include <iomanip>

using namespace std;

template <typename T>
vector<T> JacobiSolver(const SparseMatrix<T>& matrix, const vector<T>& freeElements) {
	//	assert(matrix.size() == freeElements.size());
	int N = freeElements.size();
	vector<T> result(N, 0.5);
	vector<T> temp(N);
	double norm;

	int j = 0;
	do {
		for (int i = 0; i < N; i++) {
			temp[i] = freeElements[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					temp[i] -= matrix.get(i + 1, j + 1) * result[j];
				}
			}
			temp[i] /= matrix.get(i + 1, i + 1);
		}
		norm = util::norm(result, temp);
		result.assign(temp.begin(), temp.end());
		//cout << "iteration " << j++ << endl;
		//util::print(result);
		//getchar();
	} while (norm > EPSILON);

	return result;
}