#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

VectorXd matrixToRow(MatrixXd A)
{
	int n=A.rows(),m=A.cols();
	VectorXd row(n*m);

	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			row(i*m+j)=A(i,j);

	return row;
}

int main()
{
	MatrixXd A(3,4);

	A << 1,2,3,4,5,6,7,8,9,10,11,12;

	cout << A << endl;

	VectorXd r;

	r=matrixToRow(A);
	cout << r << endl;
}

