#include "tools.h"

double f(double x) {
	//Test function.
	return 4 * x*x*x*x + 3 * x*x*x + 2 * x*x + x + 1;
}

double phi(double x)
{
	//Cumulative normal distribution function
	double L, K, w;
	//constants
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

double blackSholesValuation(double St, double X, double T, double r, Polynome sigma) {
	//Price given by the Black Scholes formula for call option.

	double result = 0;
	double m = X - St; //moneyness
	double s = sigma.eval(m, 0); //sigma is function of r and moneyness expressed in %.

	result += St*phi((log(St/X) + s*s / 2 * T) / (s*sqrt(T)));
	result -= X*exp(-r*T)*phi((log(St/X) - s*s / 2 * T) / (s*sqrt(T)));
	return result;
}

double riskNeutralDistribEval(double St, double X, double T, double r, Polynome sigma) {
	//Evaluation of the risk-neutral distribution.
	double result = 0;

	result += blackSholesValuation(St, X + delta, T, r, sigma);
	result += +blackSholesValuation(St, X - delta, T, r, sigma);
	result -= 2 * blackSholesValuation(St, X, T, r, sigma);
	return exp(r*T)*result / (delta*delta);
}

double riskNeutralCumulDistribEval(double St, double X, double T, double r, Polynome sigma) {
	//Evaluation of the cumulative distribution
	double result = 1;

	result += exp(r*T) / delta*blackSholesValuation(St, X + delta, T, r, sigma);
	result -= exp(r*T) / delta*blackSholesValuation(St, X - delta, T, r, sigma);
	return result;
}

void fillingMatrix(Eigen::MatrixXd& A, double x) {
	//Adding a row to a Eigen::MatrixXd
	assert(A.cols() == 5);
	A.conservativeResize(A.rows()+1, A.cols());
	for (int i = 0; i < 5; i++)
		A(A.rows()-1, i) = pow(x, i);
}