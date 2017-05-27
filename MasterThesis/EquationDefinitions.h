#pragma once

const double COORD_LOWER_BOUND_X = 0;
const double COORD_UPPER_BOUND_X = 5;
const double COORD_LOWER_BOUND_Y = 0;
const double COORD_UPPER_BOUND_Y = 2;

const double delta_x = 1;
const double delta_y = 1;

const double EPSILON = 0.1;

namespace def{
	template <typename T> 
	T sigma (T x, T y){
		return -2 * (x * x + y * y);
	}

	template <typename T>
	double boundaryValue (T x, T y){
		return x * x * y * y;
	}
}