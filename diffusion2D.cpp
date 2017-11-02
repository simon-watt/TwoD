#include<iostream>
#include<eigen3/Eigen/Dense>
#include<cmath>
#include<fstream>

using namespace std;
using namespace Eigen;

double initF(double x,double y)
{
//	Gaussian with mean zero and std dev one

	double pi=4.0*atan(1.0);
	double f1=1.0/sqrt(4.0*pi)*exp(-pow(x,2)/4.0);
	double f2=1.0/sqrt(4.0*pi)*exp(-pow(y,2)/4.0);

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

VectorXd initDim(double x0,double x1,int n)
{
	VectorXd x(n);
	double dx=1.0*(x1-x0)/(n-1);
	for (int i=0;i<n;i++)
		x(i)=dx*i+x0;

	return x;
}

MatrixXd colToMatrix(VectorXd col,int n,int m)
{
        MatrixXd A(n,m);

        for (int i=0;i<m;i++)
                for (int j=0;j<n;j++)
                        A(j,i)=col(i*n+j);

        return A;
}

MatrixXd rowToMatrix(VectorXd row,int n,int m)
{
        MatrixXd A(n,m);

        for (int i=0;i<n;i++)
                for (int j=0;j<m;j++)
                        A(i,j)=row(i*m+j);

        return A;
}

VectorXd triSolve(MatrixXd A,VectorXd b)
{                                       
        int n=A.rows();                 
        VectorXd sol(n),gam(n);         
        double beta;                    

//      A(i,0) = sub-diagonal (a)
//      A(i,1) = diagonal (b)    
//      A(i,2) = super-diagonal (c)

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

void integrate(MatrixXd &U,double dt,double dx,double dy)
{
	int n=U.rows(),m=U.cols();
	VectorXd rhs(n*m);
	MatrixXd T(n*m,3); // Tu = r
	MatrixXd Uh(n,m);
	VectorXd sol(n*m);

	for (int j=0;j<m;j++)
	{
		for (int i=0;i<n;i++)
		{
			int index=j*n+i;
			if (j==0 || j==m-1)
			{
				rhs(index)=0.0;
				T(index,0)=0.0;
				T(index,1)=1.0;
				T(index,2)=0.0;
			}
			else
			{
				double Uyy=(U(i,j+1)-2.0*U(i,j)+U(i,j-1))/pow(dy,2);
				rhs(index)=2.0/dt*U(i,j)+Uyy;
				T(index,0)=-1.0/pow(dx,2);
				T(index,1)=2.0/dt+2.0/pow(dx,2);
				T(index,2)=-1.0/pow(dx,2);
			}
		}
	}
	sol=triSolve(T,rhs);
	Uh=colToMatrix(sol,n,m);

	for (int i=0;i<n;i++)
	{
		for (int j=0;j<m;j++)
		{
			int index=i*m+j;
			if (i==0 || i==n-1)
			{
				rhs(index)=0.0;
				T(index,0)=0.0;
				T(index,1)=1.0;
				T(index,2)=0.0;
			}
			else
			{
				double Uxx=(Uh(i+1,j)-2.0*Uh(i,j)+Uh(i-1,j))/pow(dx,2);
				rhs(index)=2.0/dt*Uh(i,j)+Uxx;
				T(index,0)=-1.0/pow(dy,2);
				T(index,1)=2.0/dt+2.0/pow(dy,2);
				T(index,2)=-1.0/pow(dy,2);
			}
		}
	}
	sol=triSolve(T,rhs);
	U=rowToMatrix(sol,n,m);
}

int main()
{
	int N=201,M=201;
	double x0=-10.0,x1=10.0,y0=-10.0,y1=10.0;
	double dx=(x1-x0)/(N-1),dy=(y1-y0)/(M-1);
	MatrixXd U;
	VectorXd x,y;

	x=initDim(x0,x1,N);
	y=initDim(y0,y1,M);

	init(U,x,y);

	double dt=1e-2,t=0;

	ofstream out2("peak.dat");

	while (t<1)
	{
		integrate(U,dt,dx,dy);
		t+=dt;
		cout << "t = " << t << endl;
		int ii=(N-1)/2,jj=(M-1)/2;
		out2 << t << " " << x(ii) << " " << y(jj) << " " << U(ii,jj) << endl;
	}
	out2.close();

	ofstream out("data.dat");
	for (int i=0;i<N;i++)
		for (int j=0;j<M;j++)
			out << x(i) << " " << y(j) << " " << U(i,j) << endl;
	out.close();
}	
