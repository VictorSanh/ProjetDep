#include "tools.h"
#include "polynome.h"
//#include "Imagine/Graphics.h"

using namespace std;
//using namespace Imagine;




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




int main()
{
	Eigen::MatrixXd A(0,5);
	A.setZero();
	Eigen::VectorXd b(7);
	
	fillingMatrix(A, -2);
	b(0) = 36.266;
	fillingMatrix(A, -1);
	b(1) = 9.847;
	fillingMatrix(A, -0.50);
	b(2) = 2.614;
	fillingMatrix(A, 0);
	b(3) = 1.753;
	fillingMatrix(A, 0.50);
	b(4) = -0.324;
	fillingMatrix(A, 1);
	b(5) = 0.233;
	fillingMatrix(A, 2);
	b(6) = 1.879;

	cout << "A" << endl << A << endl;
	cout << "b" << endl << b << endl;

	Eigen::VectorXd coeff(5);
	coeff = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	cout << "The least-squares solution is:\n"
		<< coeff << endl;

	Polynome interpolation(coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));
	


	//graphicDisplay(interpolation);

	return 0;
}
