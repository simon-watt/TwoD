#include<iostream>
#include<eigen3/Eigen/Dense>
#include<cmath>

using namespace std;
using namespace Eigen;

double initF(double x,double y)
{
//	Gaussian with mean zero and std dev one

	double pi=4.0*atan(1.0);
	double f1=1.0/sqrt(2.0*pi)*exp(-pow(x,2)/2.0);
	double f2=1.0/sqrt(2.0*pi)*exp(-pow(y,2)/2.0);

	return f1*f2;
}

void init(MatrixXd &U,VectorXd x,VectorXd y)
{
	int N=x.size(),M=y.size();
	
	U.resize(N,M);

	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			U(i,j)=initF(x(i),y(j));
}

int main()
{
	int N=101,M=201;
	double x0=-10.0,x1=10.0,y0=-10.0,y1=10.0;
	double dx=(x1-x0)/(N-1),dy=(y1-y0)/(M-1);
	MatrixXd U;
	VectorXd x(N),y(M);

	for (int i=0;i<N;i++)
	{
		x(i)=i*dx+x0;
	}

	for (int j=0;j<M;j++)
	{
		y(j)=j*dy+y0;
	}

	init(U,x,y);

	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			cout << x(i) << " " << y(j) << " " << U(i,j) << endl;
}	
