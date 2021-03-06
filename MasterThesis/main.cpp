#include <iostream>
#include <iomanip>
#include "..\Sparse-Matrix\src\SparseMatrix\SparseMatrix.h"
#include "SimpleIteration.h"
#include "EquationDefinitions.h"
#include "ConjugateGradient.h"
#include "Util.h"

using namespace std;

void main(){

	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	const int coord_x_steps = (COORD_UPPER_BOUND_X - COORD_LOWER_BOUND_X) / delta_x + 1;
	const int coord_y_steps = (COORD_UPPER_BOUND_Y - COORD_LOWER_BOUND_Y) / delta_y + 1;
	const int matr_size = (coord_x_steps - 2) * (coord_y_steps - 2);

	//cout << "x steps: " << coord_x_steps << endl;
	//cout << "y steps: " << coord_y_steps << endl;

	//initialize u
	vector<vector<double>> u;
	u.resize(coord_x_steps, vector<double>(coord_y_steps, -1));

	//set boundary values:
	for (int x = 0; x < coord_x_steps; x++){
		u[x][0] = def::boundaryValue(COORD_LOWER_BOUND_X + x * delta_x, COORD_LOWER_BOUND_Y);
		u[x][coord_y_steps-1] = def::boundaryValue(COORD_LOWER_BOUND_X + x * delta_x, COORD_LOWER_BOUND_Y + (coord_y_steps - 1) * delta_y);
	}
	for (int y = 0; y < coord_y_steps; y++){
		u[0][y] = def::boundaryValue(COORD_LOWER_BOUND_X, COORD_LOWER_BOUND_Y + y * delta_y);
		u[coord_x_steps-1][y]	= def::boundaryValue(COORD_LOWER_BOUND_X + (coord_x_steps - 1) * delta_x, COORD_LOWER_BOUND_Y + y * delta_y);
	}

	//initialize matrix with the coefficients
	SparseMatrix<double> coeff(matr_size);
	for (int i = 0; i < matr_size; i++){
		if (i >= coord_y_steps - 2){
			coeff.set( 1/(delta_x*delta_x), i + 1, i - (coord_y_steps - 2) + 1); 
		}
		if (i >= 1 && i % (coord_y_steps -2) > 0){
			coeff.set( 1/(delta_y*delta_y), i + 1, i);
		}
		
		coeff.set ((-2/(delta_y*delta_y) - 2/(delta_x*delta_x)), i + 1, i + 1);

		if (i < matr_size - 1 && (i+1) % (coord_y_steps -2) > 0){
			coeff.set( 1/(delta_y*delta_y), i + 1, i + 2);
		}
		if (i < matr_size - (coord_y_steps - 2)){
			coeff.set( 1/(delta_x*delta_x), i + 1, i + (coord_y_steps - 2) + 1);
		}

	}
	cout << "Coefficients: " << endl;
	cout << coeff << endl;

	vector<double> freeElements(matr_size);
	for (int i = 1; i < coord_x_steps - 1; i++){
		for (int j = 1; j < coord_y_steps - 1; j++){
			int coord = (i-1)*(coord_y_steps-2) + (j-1);
			freeElements[coord] = -def::sigma(COORD_LOWER_BOUND_X + i * delta_x, COORD_LOWER_BOUND_Y + j * delta_y);
			if (j == 1){
				freeElements[coord] += -u[i][j-1] / (delta_y * delta_y);
			}
			else if (j == coord_y_steps - 2){
				freeElements[coord] += -u[i][j+1] / (delta_y * delta_y);
			}
			if (i == 1){
				freeElements[coord] += -u[i-1][j] / (delta_x * delta_x);
			}
			else if (i == coord_x_steps - 2){
				freeElements[coord] += -u[i+1][j] / (delta_x * delta_x);
			}
		}
	}
	cout << "Free elements: " << endl;
	util::print(freeElements);

	//vector<double> result = JacobiSolver(coeff, freeElements);
	vector<double> result = ConjugateGradientSolver(coeff, freeElements);
	util::print(result);

	int j = 0;
	for (int x = 1; x < coord_x_steps-1; x++){
		for (int y = 1; y < coord_y_steps-1; y++){
			u[x][y] = result[j++];
		}
	}

	cout << "Result: " << endl;
	//util::print_matrix(u);

	cout << "\t";
	for (int j = 0; j < coord_y_steps; j++) {
		cout << COORD_LOWER_BOUND_Y + j * delta_y << "\t";
	}
	cout << endl;

	for (int i = 0; i < coord_x_steps; i++) {
		cout << COORD_LOWER_BOUND_X + i * delta_x << "\t";
		for (int j = 0; j < coord_y_steps; j++) {
			cout << u[i][j] << "\t";
		}
		cout << endl;
	}
	cout << "The end!" << endl;
	getchar();
}
