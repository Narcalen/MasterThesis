#include "iostream"
#include "SparseMatrix\SparseMatrix.h"

using namespace std;
using namespace System;
using namespace System::Globalization;


void main(){
	Console::WriteLine("Hello world");
	SparseMatrix<int> m(3);
	m.set(1, 1, 3)
		.set(2, 2, 1)
		.set(3, 2, 3)
		.set(4, 3, 2)
		.set(7, 3, 3);

	m.print();
	Console::ReadKey();
}