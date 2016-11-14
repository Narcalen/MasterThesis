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
	cout << "z: " << endl;
	util::print(z);

	double norm_free_elem = norm(freeElements);

	int i = 0;

	for (i = 1; i <=N; i++){
		//while(norm(r) / norm(freeElements) >= EPSILON) {
		//cout << endl << "Iteration " << i << endl;
		Az = multiply(matrix, z);
		//cout << "Az: " << endl;
		//util::print(Az);
		alpha = scalar(r, r) / scalar(Az, z);
		//cout << "alpha = " << alpha << endl;
		result = add(result, multiply(alpha, z));
		//util::print(result);
		oldR.assign(r.begin(), r.end());
		r = add(r, multiply(-1 * alpha, Az));
		beta = scalar(r, r) / scalar(oldR, oldR);
		//cout << "beta = " << beta << endl;
		z = add(r, multiply(beta, z));
		//i++;
	}

	cout << "SLE solved in " << i << " iterations" << endl;
	return result;
}