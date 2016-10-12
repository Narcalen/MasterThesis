#include <iostream>
#include <iomanip>
#include "SparseMatrix\SparseMatrix.h"
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
		u[x][0] = def::boundaryValue(x, 0);
		u[x][coord_y_steps-1] = def::boundaryValue(x, coord_y_steps-1);
	}
	for (int y = 0; y < coord_y_steps; y++){
		u[0][y] = def::boundaryValue(0, y);
		u[coord_x_steps-1][y]	= def::boundaryValue(coord_x_steps-1, y);
	}

	//initialize matrix with the coefficients
	SparseMatrix<double> coeff(matr_size);
	for (int i = 0; i < matr_size; i++){
		if (i >= coord_y_steps - 2){
			coeff.set( 1/(delta_x*delta_x), i, i - (coord_y_steps - 2)); 
		}
		if (i >= 1 && i % (coord_y_steps -2) > 0){
			coeff.set( 1/(delta_y*delta_y), i, i - 1); 
		}
		if (delta_x != delta_y){ //Temporary workaround as SparseMatrix class does not handle 0 insertion properly
			coeff.set ((2/(delta_y*delta_y) - 2/(delta_x*delta_x)), i, i);
		}
		if (i < matr_size - 1 && (i+1) % (coord_y_steps -2) > 0){
			coeff.set( 1/(delta_y*delta_y), i, i + 1); 
		}
		if (i < matr_size - (coord_y_steps - 2)){
			coeff.set( 1/(delta_x*delta_x), i, i + (coord_y_steps - 2)); 
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

	//Jacobi solver does not work as the matrix is not diagonally dominant
	//vector<double> result = JacobiSolver(coeff, freeElements);
	vector<double> result = ConjugateGradientSolver(coeff, freeElements);
	//util::print(result);

	int j = 0;
	for (int x = 1; x < coord_x_steps-1; x++){
		for (int y = 1; y < coord_y_steps-1; y++){
			u[x][y] = result[j++];
		}
	}

	cout << "Result: " << endl;
	util::print_matrix(u);
	cout << "The end!" << endl;
	getchar();
}