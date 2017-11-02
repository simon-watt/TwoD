#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd colToMatrix(VectorXd col,int n,int m)
{
	MatrixXd A(n,m);

	for (int i=0;i<m;i++)
		for (int j=0;j<n;j++)
			A(j,i)=col(i*n+j);

	return A;
}

int main()
{
	int n=3,m=4;
	VectorXd c(n*m);
	MatrixXd A;

	c << 1,5,9,2,6,10,3,7,11,4,8,12;

	cout << c << endl;

	A=colToMatrix(c,n,m);
	cout << A << endl;
}

