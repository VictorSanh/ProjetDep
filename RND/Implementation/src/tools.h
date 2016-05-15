#pragma once

#include <math.h> 
#include <fstream>
#include "polynome.h"

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif

const double delta = 0.0001;

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



///////////////////////
//For Graphical Purpose
///////////////////////

//Writing the interpolated smile volatility to a csv file.
void writingSmileToCSV(std::string nom_fichier, std::vector <double> abscisse, Polynome vol);

//Writing the cumulated distribution to a csv file.
void writingCumulDistrivToCSV(std::string nom_fichier, std::vector <double> abscisse, Polynome vol, double St, double T);

//Writing all the information in CSVs files for graphical purpose (using plot from python).
void writingToCsv(double St, double T, Polynome vol);