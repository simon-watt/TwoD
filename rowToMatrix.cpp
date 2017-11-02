#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd rowToMatrix(VectorXd row,int n,int m)
{
	MatrixXd A(n,m);

	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			A(i,j)=row(i*m+j);

	return A;
}

int main()
{
	int n=3,m=4;
	VectorXd r(n*m);
	MatrixXd A;

	r << 1,2,3,4,5,6,7,8,9,10,11,12;

	cout << r << endl;

	A=rowToMatrix(r,n,m);
	cout << A << endl;
}

