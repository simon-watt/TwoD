#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

VectorXd matrixToCol(MatrixXd A)
{
	int n=A.rows(),m=A.cols();
	VectorXd col(n*m);

	for (int i=0;i<m;i++)
		for (int j=0;j<n;j++)
			col(i*n+j)=A(j,i);

	return col;
}

int main()
{
	MatrixXd A(3,4);

	A << 1,2,3,4,5,6,7,8,9,10,11,12;

	cout << A << endl;

	VectorXd c;

	c=matrixToCol(A);
	cout << c << endl;
}

