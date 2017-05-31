#pragma once
#include "Node.h"

double detD(double x1, double y1, double x2, double y2, double x3, double y3) {
	return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
}

double detD(Node node1, Node node2, Node node3) {
	return detD(node1.get(0), node1.get(1), node2.get(0), node1.get(1), node3.get(0), node3.get(1));
}

double delta(double coord1, double coord2) {
	return coord1 - coord2;
}

double deltaX(Node node1, Node node2) {
	return delta(node1.get(0), node2.get(0));
}

double deltaY(Node node1, Node node2){
	return delta(node1.get(1), node2.get(1));
}

