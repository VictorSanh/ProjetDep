#include <iostream>
#include <Eigen/Dense>

using namespace std;


float f(float x) {
	return x*x*x*x + x*x*x + x*x + x + 1;
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
	cout << "The least-squares solution is:\n"
		<< A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b) << endl;
	return 0;
}
