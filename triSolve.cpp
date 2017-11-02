#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

VectorXd triSolve(MatrixXd A,VectorXd b)
{
	int n=A.rows();
	VectorXd sol(n),gam(n);
	double beta;

//	A(i,0) = sub-diagonal (a)
//	A(i,1) = diagonal (b)
//	A(i,2) = super-diagonal (c)

	beta=A(0,1);
	sol(0)=b(0)/beta;
	for (int j=1;j<n;j++)
	{
		gam(j)=A(j-1,2)/beta;
		beta=A(j,1)-A(j,0)*gam(j);
		sol(j)=(b(j)-A(j,0)*sol(j-1))/beta;
	}

	for (int j=n-2;j>=0;j--)
		sol(j)-=gam(j+1)*sol(j+1);

	return sol;
}

MatrixXd triToFull(MatrixXd T)
{
	int n=T.rows();
	MatrixXd A(n,n);

	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			A(i,j)=0;

	for (int i=0;i<n;i++)
	{
		A(i,i)=T(i,1);
		if (i>0)
			A(i,i-1)=T(i,0);
		if (i<n-1)
			A(i,i+1)=T(i,2);
	}

	return A;
}

int main()
{
	int n=3;
	MatrixXd A(n,3);
	VectorXd b(n);

	A << 1,1,2,1,3,3,2,4,1;

	cout << A << endl;

	cout << triToFull(A) << endl;
	
	b << 1,2,3;

	cout << b << endl;

	VectorXd sol=triSolve(A,b);
	cout << sol << endl;

	cout << "check" << endl << triToFull(A)*sol << endl;
}
	

	
