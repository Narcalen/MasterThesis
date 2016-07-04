#include "Util.h"
#include <iostream>
#include <math.h>

using namespace std;

template<typename T> void print (ostream & os, vector<T>* vector){
	for (int i = 0; i < vector -> size(); i++){
		os << vector -> at(i) << "\t";
	}

	os << endl;
}

template<typename T> T distance (vector <T>* vector1, vector <T>* vector2){
	T distance = 0;
	for (int i = 0; i < vector1 -> size(); i++){
		distance += pow(vector1 -> at(i) - vector2 -> at(i), 2);
	}
	return sqrt(fabs(distance));
}

template void print (ostream & os, vector<int>* vector);
template void print (ostream & os, vector<float>* vector);
template void print (ostream & os, vector<double>* vector);

template float distance (vector <float>* vector1, vector <float>* vector2);
template double distance (vector <double>* vector1, vector <double>* vector2);