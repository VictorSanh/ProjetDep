#pragma once

#include <iostream>
#include <vector> 
#include <../Eigen/Dense>

class Polynome {
	//Polynome Class
public:
	//Constructor : We make the assumption that the polynome we are looking for is degree 4.
	Polynome(double a0, double a1, double a2, double a3, double a4);

	//Vector of coefficients. coeff[i] is the coeff before x^i
	std::vector<double> coeff;

	//Evaluate the trained polynom in point x. Degree indicates the degree of derivate.
	double eval(double x, int degree);
};

