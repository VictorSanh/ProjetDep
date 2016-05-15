#include "tools.h"
#include "polynome.h"

using namespace std;

int main()
{
	Eigen::MatrixXd A(0,5);
	A.setZero();
	Eigen::VectorXd b(9);
	
	fillingMatrix(A, -2);
	//b(0) = 36.266;
	b(0) = 32.5790;
	fillingMatrix(A, -1);
	//b(1) = 9.847;
	b(1) = 28.934;
	fillingMatrix(A, -0.50);
	//b(2) = 2.614;
	b(2) = 27.8261;
	fillingMatrix(A, 0);
	//b(3) = 1.753;
	b(3) = 27.0250;
	fillingMatrix(A, 0.50);
	//b(4) = -0.324;
	b(4) = 26.4866;
	fillingMatrix(A, 1);
	//b(5) = 0.233;
	b(5) = 26.1151;
	fillingMatrix(A, 2);
	//b(6) = 1.879;
	b(6) = 25.7388;
	fillingMatrix(A, -0.25);
	b(7) = 27.3975;
	fillingMatrix(A, 0.25);
	b(8) = 26.7361;

	cout << "A" << endl << A << endl;
	cout << "b" << endl << b << endl;

	Eigen::VectorXd coeff(5);
	coeff = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	cout << "The least-squares solution is:\n"
		<< coeff << endl <<endl;

	Polynome interpolation(coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));

	
	double St = 4.09;
	/*cout << "-2 : " << riskNeutralCumulDistribEval(St, St+4, 1, 0, interpolation) << endl;
	cout << "-2 : " << riskNeutralCumulDistribEval(St, St - 3, 1, 0, interpolation) << endl;
	cout << "-2 : " << riskNeutralCumulDistribEval(St, St - 4, 1, 0, interpolation) << endl;
	cout << "-2 : " << riskNeutralCumulDistribEval(St, St - 4, 1, 0, interpolation) << endl;
	cout << "-2 : " << riskNeutralCumulDistribEval(St, St - 4, 1, 0, interpolation) << endl;*/

	//cout << "BLACK SCHOLES VALUATION" << endl;
	//cout << blackSholesValuation(St, 2.09, 2, 0, interpolation) << endl;
	//cout << blackSholesValuation(St, 3.09, 2, 0, interpolation) << endl;
	//cout << blackSholesValuation(St, 1.09, 2, 0, interpolation) << endl;
	//cout << blackSholesValuation(St, 0.09, 2, 0, interpolation) << endl;

	writingToCsv(St, 2, interpolation);
	return 0;
}
