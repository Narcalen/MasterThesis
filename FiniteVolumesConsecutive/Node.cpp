#include "Node.h"
#include <cmath>


Node::Node()
{
	Node(2);
}

Node::Node(int dim)
{
	this->dim = dim;
	coordinates= std::vector<double>(dim);
}


Node::~Node()
{
	//delete coordinates;
}

double Node::get(int num) {
	if (num >= 0 && num < dim) {
		return coordinates[num];
	}
	else {
		return NAN;
	}
}

void Node::set(int num, int value) {
	if (num >= 0 && num < dim) {
		coordinates[num] = value;
	}
}