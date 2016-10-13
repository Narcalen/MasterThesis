#pragma once

const int COORD_LOWER_BOUND_X = 0;
const int COORD_UPPER_BOUND_X = 5;
const int COORD_LOWER_BOUND_Y = 0;
const int COORD_UPPER_BOUND_Y = 2;

const double delta_x = 1;
const double delta_y = 0.5;

const double EPSILON = 0.0001;

namespace def{
	template <typename T> 
	T sigma (T x, T y){
		return (T) 0;
	}

	template <typename T>
	double boundaryValue (T x, T y){
		return 2 * (x * x + y * y);
	}
}