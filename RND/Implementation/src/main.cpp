#include <iostream>
#include <vector>
#include <math.h>   
#include <../Eigen/Dense>
//#include "Imagine/Graphics.h"

using namespace std;
//using namespace Imagine;

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif

const double delta = 0.00001;


class Polynome {
	//Polynome Class
public:
	//Constructor : We make the assumption that the polynome we are looking for is degree 4.
	Polynome(double a0, double a1, double a2, double a3, double a4);

	//Vector of coefficients. coeff[i] is the coeff before x^i
	vector<double> coeff;

	//Evaluate the trained polynom in point x. Degree indicates the degree of derivate.
	double eval(double x, int degree);
};

Polynome::Polynome(double a0, double a1, double a2, double a3, double a4) {
	coeff.push_back(a0);
	coeff.push_back(a1);
	coeff.push_back(a2);
	coeff.push_back(a3);
	coeff.push_back(a4);
}

void graphicDisplay(Polynome interp) {
	//This function aims to graphically display a polynome
	
	//Window W = openWindow(512, 512, "2D basics");
	/*for (int i = 0; i < 20; i += 2)
		drawPoint(i + 100, i + 200, BLACK);
	drawArrow(100, 470, 200, 450, BLUE);                    // An arrow from (100,450) to (200,450)
	drawArrow(300, 470, 200, 450, RED, 30, 10, 1);             // tip=30x10 pixels, style=1
	drawArrow(200, 450, 250, 400, BLACK, 20, 20, 2);
	click();*/
	//closeWindow(W);
}

double Polynome::eval(double x, int degree) {
	if (coeff.size() != 5) {
		cout << "The function has not been trained. Returning 0." << endl;
		return 0;
	}
	else {
		if (degree == 0)
			return coeff[0] + coeff[1] * x + coeff[2] * x * x + coeff[3] * x * x * x + coeff[4] * x * x * x * x;
		else if (degree == 1)
			return coeff[1] + 2 * coeff[2] * x + 3 * coeff[3] * x * x + 4 * coeff[4] * x * x * x;
		else if (degree == 2)
			return 2 * coeff[2] + 6 * coeff[3] * x + 12 * coeff[4] * x * x;
		else if (degree == 3)
			return  6 * coeff[3] + 24 * coeff[4] * x;
		else {
			cout << "Are you sure you want such a high degree evaluation ? Returning 0." << endl;
			return 0;
		}
	}
}

double f(double x) {
	//Test function.
	return 4*x*x*x*x + 3*x*x*x + 2*x*x + x + 1;
}

double phi(double x)
{
	//Cumulative normal distribution function
	double L, K, w;
	/* constants */
	double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
	double const a4 = -1.821255978, a5 = 1.330274429;

	L = fabs(x);
	K = 1.0 / (1.0 + 0.2316419 * L);
	w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K, 3) + a4 * pow(K, 4) + a5 * pow(K, 5));

	if (x < 0) {
		w = 1.0 - w;
	}
	return w;
}

double blackSholesValuation(double St, double X, double T, Polynome sigma) {
	double result = 0;
	double s = sigma.eval(St / X, 0);
	result += St*phi((log(St / X) + s*s / 2 * T)/(s*sqrt(T)));
	result -= X*phi((log(St / X) - s*s / 2 * T) / (s*sqrt(T)));
	return result;
}

double riskNeutralDistribEval(double St, double X, double T, Polynome sigma) {
	double result = 0;
	blackSholesValuation(St, X + delta, T, sigma) + blackSholesValuation(St, X - delta, T, sigma) - 2 * blackSholesValuation(St, X, T, sigma);
	return result / (delta*delta);
}

int main()
{
	Eigen::MatrixXd A(7, 5);
	A.setZero();
	Eigen::VectorXd b(7);
	for (int x = 0; x < 7; x++) {
		b(x) = f(x);
		A(x, 0) = 1;
		A(x, 1) = x;
		A(x, 2) = x*x;
		A(x, 3) = x*x*x;
		A(x, 4) = x*x*x*x;
	}

	cout << "A" << endl << A << endl;
	cout << "b" << endl << b << endl;

	Eigen::VectorXd coeff(5);
	coeff = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	cout << "The least-squares solution is:\n"
		<< coeff << endl;

	Polynome interpolation(coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));
	cout << "Test : Evaluation in x = 0 : " << interpolation.eval(1,0) << endl;
	cout << "Test : blackScholesEvaluation : " << blackSholesValuation(1, 2, 3, interpolation) << endl;
	//graphicDisplay(interpolation);

	return 0;
}
