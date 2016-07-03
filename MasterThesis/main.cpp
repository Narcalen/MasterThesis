#include "iostream"
#include "SparseMatrix\SparseMatrix.h"

using namespace std;
using namespace System;
using namespace System::Globalization;


void main(){
	Console::WriteLine("Hello world");
	SparseMatrix<int> m(3);
	m.set(1, 0, 2)
		.set(2, 1, 0)
		.set(3, 1, 2)
		.set(4, 2, 1)
		.set(7, 2, 2);

	m.print();
	Console::ReadKey();
}