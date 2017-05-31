#include <iostream>
#include <iomanip>
#include "..\MasterThesis\EquationDefinitions.h"

using namespace std;

void main()
{
	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	const int coord_x_steps = (COORD_UPPER_BOUND_X - COORD_LOWER_BOUND_X) / delta_x + 1;
	const int coord_y_steps = (COORD_UPPER_BOUND_Y - COORD_LOWER_BOUND_Y) / delta_y + 1;

	/*cout << "\t";
	for (int j = 0; j < coord_y_steps; j++) {
		cout << COORD_LOWER_BOUND_Y + j * delta_y << "\t";
	}
	cout << endl;

	for (int i = 0; i < coord_x_steps; i++) {
		cout << COORD_LOWER_BOUND_X + i * delta_x << "\t";
		for (int j = 0; j < coord_y_steps; j++) {
			cout << def::boundaryValue(COORD_LOWER_BOUND_X + i * delta_x, COORD_LOWER_BOUND_Y + j * delta_y) << "\t";
		}
		cout << endl;
	}*/

	cout << "x\ty\tvalue" << endl;
	for (int i = 0; i < coord_x_steps; i++) {
		for (int j = 0; j < coord_y_steps; j++) {
			cout << COORD_LOWER_BOUND_X + i * delta_x  << "\t" << COORD_LOWER_BOUND_Y + j * delta_y << "\t" << def::boundaryValue(COORD_LOWER_BOUND_X + i * delta_x, COORD_LOWER_BOUND_Y + j * delta_y) << endl;
		}
	}

	cout << "The end!" << endl;
	getchar();
}

