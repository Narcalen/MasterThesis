#include "Element.h"
#include <cmath>


Element::Element(int dim)
{
	this -> dim = dim;
	nodes = std::vector<int>(dim);
}


Element::~Element()
{
}

int Element::getNode(int number) {
	if (number >= 0 && number < dim) {
		return nodes[number];
	}
	else {
		return NAN;
	}
}

void Element::setNode(int number, int value) {
	if (number >= 0 && number < dim) {
		nodes[number] = value;
	}
}