#pragma once
#include <vector>

class Element
{
public:
	Element(int dim);
	~Element();
	int getNode(int number);
	void setNode(int number, int value);

private:
	int dim;
	std::vector<int> nodes;
};

