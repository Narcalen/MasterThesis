#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iostream>

#include "..\MasterThesis\EquationDefinitions.h"
#include "..\MasterThesis\ConjugateGradient.h"
#include "..\MasterThesis\SimpleIteration.h"
#include "..\Sparse-Matrix\src\SparseMatrix\SparseMatrix.h"
#include "..\MasterThesis\Util.h"
#include "Element.h"
#include "Node.h"
#include "FVMUtil.h"


using namespace std;

void main()
{
	//output format
	cout << fixed << showpoint;
	cout << setprecision(2);

	int a, b, c;
	ifstream infile;

	vector<Element> elements;
	vector<Node> nodes;

	//read nodes data
	infile.open("distmesh_nodes.txt", ios::in);
	if (!infile) {
		cout << "Cannot open input file distmesh_nodes.txt" << endl;
	}

	Node node(2);
	while (infile >> a >> b)
	{
		node.set(0, a);
		node.set(1, b);
		nodes.push_back(node);
	}
	infile.close();
	//cout << "Nodes: " << endl;
	//for (Node node: nodes) {
	//	cout << node.get(0) << " " << node.get(1) << endl;
	//}

	//read elements data
	infile.open("distmesh_elements.txt", ios::in);
	if (!infile) {
		cout << "Cannot open input file distmesh_elements.txt" << endl;
	}

	Element element(3);
	while (infile >> a >> b >> c)
	{
		element.setNode(0, a);
		element.setNode(1, b);
		element.setNode(2, c);
		elements.push_back(element);
	}
	infile.close();
	//cout << "Elements: " << endl;
	//for (Element element : elements) {
	//	cout << element.getNode(0) << " " << element.getNode(1) << " " << element.getNode(2) << endl;
	//}

	//initialize u
	vector<double> u;
	const int matr_size = nodes.size();
	cout << "Number of elements: " << matr_size << endl;

	//initialize matrix with the coefficients
	SparseMatrix<double> coeff(matr_size);
	vector<double> freeElements(matr_size, 0);

	for (Element element : elements) {
		int node1Num = element.getNode(0);
		Node node1 = nodes[node1Num];
		int node2Num = element.getNode(1);
		Node node2 = nodes[node2Num];
		int node3Num = element.getNode(2);
		Node node3 = nodes[node3Num];
		double multiplier = 1 / (2 * abs(detD(node1, node2, node3)));

		node1Num++;
		node2Num++;
		node3Num++;

		//add local main diagonal to global main diagonal
		double val = multiplier * (deltaX(node2, node3)*deltaX(node2, node3) + deltaY(node2, node3) * deltaY(node2, node3));
		coeff.set(coeff.get(node1Num, node1Num) + val, node1Num, node1Num);

		val = multiplier * (deltaX(node3, node1)*deltaX(node3, node1) + deltaY(node3, node1) * deltaY(node3, node1));
		coeff.set(coeff.get(node2Num, node2Num) + val, node2Num, node2Num);

		val = multiplier * (deltaX(node1, node2)*deltaX(node1, node2) + deltaY(node1, node2) * deltaY(node1, node2));
		coeff.set(coeff.get(node3Num, node3Num) + val, node3Num, node3Num);

		//add local elements above and below main diagonal to global diagonals
		val = multiplier * (deltaX(node3, node1)*deltaX(node2, node3) + deltaY(node3, node1) * deltaY(node2, node3));
		coeff.set(coeff.get(node1Num, node2Num) + val, node1Num, node2Num);
		coeff.set(coeff.get(node2Num, node1Num) + val, node2Num, node1Num);

		val = multiplier * (deltaX(node1, node2)*deltaX(node2, node3) + deltaY(node1, node2) * deltaY(node2, node3));
		coeff.set(coeff.get(node1Num, node3Num) + val, node1Num, node3Num);
		coeff.set(coeff.get(node3Num, node1Num) + val, node3Num, node1Num);


		val = multiplier * (deltaX(node1, node2)*deltaX(node3, node1) + deltaY(node1, node2) * deltaY(node3, node1));
		coeff.set(coeff.get(node3Num, node2Num) + val, node3Num, node2Num);
		coeff.set(coeff.get(node2Num, node3Num) + val, node2Num, node3Num);

		//set values for the right side of the equations system
		double f1 = -def::sigma(node1.get(0), node1.get(1));
		double f2 = -def::sigma(node2.get(0), node2.get(1));
		double f3 = -def::sigma(node3.get(0), node3.get(1));
		multiplier = abs(detD(node1, node2, node3)) / 216;

		freeElements[node1Num-1] += multiplier * (22 * f1 + 7 * f2 + 7 * f3);
		freeElements[node2Num-1] += multiplier * (7 * f1 + 22 * f2 + 7 * f3);
		freeElements[node3Num-1] += multiplier * (7 * f1 + 7 * f2 + 22 * f3);
	}

	//boundary conditions 1st type
	for (unsigned i = 1; i <= nodes.size(); i++) {
		Node node = nodes[i-1];
		if (node.get(0) == COORD_LOWER_BOUND_X || node.get(0) == COORD_UPPER_BOUND_X
			|| node.get(1) == COORD_LOWER_BOUND_Y || node.get(1) == COORD_UPPER_BOUND_Y) {
			for (unsigned j = 1; j <= matr_size; j++) {
				if (coeff.get(i, j) != 0) {
					coeff.set(0, i, j);
				}
				if (coeff.get(j, i) != 0) {
					coeff.set(0, j, i);
				}
			}
			coeff.set(1, i, i);
			freeElements[i - 1] = def::boundaryValue(node.get(0), node.get(1));
		}
	}

	cout << "Coefficients: " << endl;
	cout << coeff << endl;

	cout << "Free elements: " << endl;
	util::print(freeElements);

	cout << "Result: " << endl;
	//vector<double> result = JacobiSolver(coeff, freeElements);
	vector<double> result = ConjugateGradientSolver(coeff, freeElements);

	cout << "x\ty\tvalue" << endl;
	for (unsigned i = 0; i < nodes.size(); i++) {
		Node node = nodes[i];
		cout << node.get(0) << "\t" << node.get(1) << "\t" << result[i] << endl;
	}

	cout << "The end!" << endl;
	getchar();
}