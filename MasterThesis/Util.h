#pragma once
#include <vector>
#include <iostream>

using namespace std;

template<typename T>
void print (ostream & os, vector <T>* vector);

template<typename T>
T distance (vector <T>* vector1, vector <T>* vector2);