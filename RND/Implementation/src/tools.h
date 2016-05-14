#pragma once

#include <math.h> 
#include "polynome.h"

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif

const double delta = 0.00001;

//Test function
double f(double x);

//Cumulative normal distribution function
double phi(double x);

//Price given by the Black Scholes formula for call option.
double blackSholesValuation(double St, double X, double T, double r, Polynome sigma);

//Evaluation of the risk-neutral distribution.
double riskNeutralDistribEval(double St, double X, double T, double r, Polynome sigma);

//Evaluation of the cumulative distribution
double riskNeutralCumulDistribEval(double St, double X, double T, double r, Polynome sigma);

//Adding a row to a Eigen::MatrixXd
void fillingMatrix(Eigen::MatrixXd& A, double x);