#pragma once
#include <vector>
#include <iostream>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iterator>

using namespace std;

namespace util {

	template<typename T>
	void print(ostream& os, const vector<T>& vector) {
		for (unsigned i = 0; i < vector.size(); i++) {
			os << vector[i] << "\t";
		}

		os << endl;
	}

	template<typename T>
	void print(const vector<T>& vector) {
		print(cout, vector);
	}

	template<typename T>
	void print_matrix(const vector<vector<T>>& matrix) {
		int rows = matrix.size();
		int cols = matrix[0].size();

		for (int x = 0; x < rows; x++) {
			for (int y = 0; y < cols; y++) {
				cout << matrix[x][y] << "\t";
			}
			cout << endl;
		}
	}

	template<typename T>
	T distance(const vector<T>& vector1, const vector<T>& vector2) {
		assert(vector1.size() == vector2.size());
		T distance = 0;
		for (int i = 0; i < vector1.size(); i++) {
			distance += pow(vector1[i] - vector2[i], 2);
		}
		return sqrt(fabs(distance));
	}

	template<typename T>
	T norm(const vector<T>& vector1, const vector<T>& vector2) {
		assert(vector1.size() == vector2.size());
		T norm = fabs(vector1[0] - vector2[0]);
		for (int i = 1; i < vector1.size(); i++) {
			if (fabs(vector1[i] - vector2[i]) > norm) {
				norm = fabs(vector1[i] - vector2[i]);
			}
		}
		return norm;
	}

	template<typename T>
	T sum(const vector<T>& vector1, const vector<T>& vector2) {
		assert(vector1.size() == vector2.size());
		vector<T> result(vector1.size());
		for (int i = 0; i < vector1.size(); i++) {
			result[i] = vector1[i] + vector2[i];
		}

		return result;
	}

	template<typename T>
	T scalar(const vector<T>& vector1, const vector<T>& vector2) {
		assert(vector1.size() == vector2.size());
		T result = 0;
		for (int i = 0; i < vector1.size(); i++) {
			result += vector1[i] * vector2[i];
		}

		return result;
	}

	template<typename T>
	vector<T> multiply (const SparseMatrix<T>& matrix, const vector<T>& vector) {
		int rows = matrix.size();
		int cols = matrix.size();
		assert(vector.size() == cols);
		std::vector<T> result(rows, 0);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				result[i] += matrix.get(i,j) * vector[j];
			}
		}

		return result;
	}


	template <typename T>
	vector<T> multiply(const T number, const std::vector<T>& vector)
	{
		std::vector<T> result;
		result.reserve(vector.size());

		std::transform(vector.begin(), vector.end(), std::back_inserter(result), bind1st(multiplies<T>(), number));
		return result;
	}

	template <typename T>
	vector<T> add (const std::vector<T>& a, const std::vector<T>& b)
	{
		assert(a.size() == b.size());

		std::vector<T> result;
		result.reserve(a.size());

		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
		return result;
	}


}