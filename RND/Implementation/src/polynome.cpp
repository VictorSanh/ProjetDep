#include "polynome.h"

Polynome::Polynome(double a0, double a1, double a2, double a3, double a4) {
	coeff.push_back(a0);
	coeff.push_back(a1);
	coeff.push_back(a2);
	coeff.push_back(a3);
	coeff.push_back(a4);
}

double Polynome::eval(double x, int degree) {
	if (coeff.size() != 5) {
		std::cout << "The function has not been trained. Returning 0." << std::endl;
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
			std::cout << "Are you sure you want such a high degree evaluation ? Returning 0." << std::endl;
			return 0;
		}
	}
}
