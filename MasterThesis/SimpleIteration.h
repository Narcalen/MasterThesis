#pragma once
#include "SparseMatrix\SparseMatrix.h"
#include <vector>
#include <iostream>

using namespace std;

const double EPSILON = 0.0001;

template <typename T> 
vector<T>* JacobiSolver(SparseMatrix<T> matrix, vector<T>* freeElements);