#pragma once
#include "Util.h"
#include "SparseMatrix\SparseMatrix.h"
#include "EquationDefinitions.h"
#include <vector>

using namespace std;

template <typename T> 
vector<T> ConjugateGradientSolver(const SparseMatrix<T>& matrix, const vector<T>& freeElements){
	assert(matrix.size() == freeElements.size());
	int N = matrix.size();
	vector<T> result(N,0), d(N,0), prevG(N,0), curG(N,0);
	T s;

	prevG.assign(freeElements.begin(), freeElements.end());
	prevG *= -1;

	for (i = 1; i <=N; i++){
		curG = matrix * result + (-1) * freeElements;
		d *= (-1) * curG  + util::scalar(curG,curG)/util::scalar(prevG,prevG);
		s = util::scalar(d, curG) / util::scalar(d, matrix * d);
		result += s * d;
		prevG.assign(curG.begin(), curG.end());
	}

	return result;
}