#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;

class Polynome {
public:
	Polynome(float a0, float a1, float a2, float a3, float a4);
	vector<float> coeff;
	float eval(float x, int degree);
};

Polynome::Polynome(float a0, float a1, float a2, float a3, float a4) {
	coeff.push_back(a0);
	coeff.push_back(a1);
	coeff.push_back(a2);
	coeff.push_back(a3);
	coeff.push_back(a4);
}


float Polynome::eval(float x, int degree) {
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

float f(float x) {
	return 4*x*x*x*x + 3*x*x*x + 2*x*x + x + 1;
}

int main()
{
	Eigen::MatrixXf A(7, 5);
	A.setZero();
	Eigen::VectorXf b(7);
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

	Eigen::VectorXf coeff(5);
	coeff = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	cout << "The least-squares solution is:\n"
		<< coeff << endl;

	Polynome interpolation(coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));
	cout << "Test : Evaluation in x = 0 : " << interpolation.eval(1,0) << endl;

	return 0;
}
