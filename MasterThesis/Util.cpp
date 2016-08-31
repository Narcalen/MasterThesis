#include "Util.h"
#include <iostream>
#include <math.h>

using namespace std;




template void util::print (ostream & os, vector<int> vector);
template void util::print (ostream & os, vector<float> vector);
template void util::print (ostream & os, vector<double> vector);

template void util::print (vector<int> vector);
template void util::print (vector<float> vector);
template void util::print (vector<double> vector);

template void util::print_matrix (vector <vector<int>> vector);
template void util::print_matrix (vector <vector<float>> vector);
template void util::print_matrix (vector <vector<double>> vector);

template float util::distance (vector <float> vector1, vector <float> vector2);
template double util::distance (vector <double> vector1, vector <double> vector2);

