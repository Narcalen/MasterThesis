#pragma once
#include "Util.h"
#include "SparseMatrix\SparseMatrix.h"
#include "EquationDefinitions.h"
#include <vector>

using namespace std;
using namespace util;

template <typename T> 
vector<T> ConjugateGradientSolver(const SparseMatrix<T>& matrix, const vector<T>& freeElements){
	assert(matrix.size() == freeElements.size());
	int N = matrix.size();
	vector<T> result(N, 1), z(N), oldR(N), Az(N), r(N);
	T alpha, beta;

	r = add(freeElements, multiply((T)-1, multiply(matrix, result)));
	z.assign(r.begin(), r.end());

	for (int i = 1; i <=N; i++){
		Az = multiply(matrix, z);
		alpha = scalar(r, r) / scalar(Az, z);
		result = add(result, multiply(alpha, z));
		oldR.assign(r.begin(), r.end());
		r = add(r, multiply(-1 * alpha, Az));
		beta = scalar(r, r) / scalar(oldR, oldR);
		z = add(r, multiply(beta, z));
	}

	return result;
}