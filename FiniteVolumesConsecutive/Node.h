#pragma once
#include <vector>

class Node
{
public:
	Node();
	Node(int dim);
	~Node();
	double get(int num);
	void set(int num, int value);

private:
	int dim;
	std::vector<double> coordinates;
};

