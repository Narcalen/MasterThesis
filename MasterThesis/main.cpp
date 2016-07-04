#include "iostream"
#include "iomanip"
#include "SparseMatrix\SparseMatrix.h"
#include "SimpleIteration.h"
#include "Util.h"

using namespace std;

void main(){

	cout << fixed << showpoint;
    cout << setprecision(4);

	cout << "Hello world" << endl;
	SparseMatrix<double> m(3);
	m.set(10,0,0).
		set(1,0,1).
		set(-1,0,2).
		set(1,1,0).
		set(10,1,1).
		set(-1,1,2).
		set(-1,2,0).
		set(1,2,1).
		set(10,2,2);

	vector<double>* freeElements = new vector<double>(3);
	freeElements -> at(0) = 11;
	freeElements -> at(1) = 10;
	freeElements -> at(2) = 10;

	//cout << m(1,0) << endl;
	cout << m << endl;

	vector<double>* result = JacobiSolver(m, freeElements);
	print(cout, result);
	
	getchar();
}