
void RearrangeMatrix(int n,double **A, double *b) {
//    int n = A.size(); 
    double **Ab;
    Ab = new double *[n];
    for(int i=0; i<n; i++) Ab[i] = new double[n+1];
    
    for(int i=0;i<n;i++) {
    	for(int j=0;j<n;j++) {
    		Ab[i][j] = A[i][j];
    	}
    	Ab[i][n] = b[i];
    }
    		
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(Ab[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(Ab[k][i]) > maxEl) {
                maxEl = abs(Ab[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = Ab[maxRow][k];
            Ab[maxRow][k] = Ab[i][k];
            Ab[i][k] = tmp;
        }
	}
    for(int i=0;i<n;i++) {
    	for(int j=0;j<n;j++) {
    		A[i][j] = Ab[i][j];
    	}
    	b[i] = Ab[i][n];
    }

}

void Jacobi(int N,double **A, double *b, double *u,  double tol)
{
	double *uold;
	uold = new double[N];
	for(int i = 0; i<N;i++)
	    uold[i] = u[i];

  int iter = 0;
  double errMax = 1.0;
	while(errMax > tol) {
		for (int i = 0; i < N; i++) {		
			double sum = 0.;		
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += A[i][j] * uold[j];
				}
			}
			u[i] = (b[i] - sum) / (A[i][i]);
		}	
		errMax =0.0;
		for(int i=0;i<N;i++) {
			double error = fabs((uold[i] - u[i]))/fabs(uold[i]);
			if(error > errMax) errMax = error;
		}		
		for(int i=0; i<N;i++) {
		    uold[i] = u[i];
		}
		iter++;
	}
}

void GaussSe(int N,double **A, double *b, double *u,  double tol)
{
    double *uold;
	uold = new double[N];
	for(int i = 0; i<N;i++) uold[i] = u[i];

  RearrangeMatrix(N, A, b);

  int iter = 0;
  double errMax = 1.0;
	while(errMax > tol) {
		for (int i = 0; i < N; i++) {		
			double sum = 0.;		
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += A[i][j] * u[j];
				}
			}
//			std::cout << "A " << i << ", " << A[i][i] << std::endl;
			u[i] = (b[i] - sum) / (A[i][i]);
//			std::cout << "b sum " << i << ", " << b[i]  << ", " << sum << std::endl;
//			std::cout << "u " << u[i] << std::endl;
		}	
		errMax =0.0;
		for(int i=0;i<N;i++) {
			double error = fabs(uold[i] - u[i])/fabs(uold[i]);
			if(error > errMax) errMax = error;
		}		
		for(int i=0; i<N;i++) {
		    uold[i] = u[i];
		}
		iter++;
	}
}

void GaussSe2(int N,double **A, double *b, double *x,  double tol)
{
    int iter=0;
    double errMax = 1.0;
	double *c, *xold;
    c = new double[N];
    xold = new double[N];
	for(int i=0;i<N;i++) xold[i] = x[i];
    RearrangeMatrix(N, A, b);

    while(errMax > tol) {
		iter++;
		for(int i=0;i<N;i++) {
	    	c[i] = b[i] /A[i][i];
	    	for(int j=0;j<N;j++) {
				if(j == i) continue;
				c[i] -= (A[i][j]/A[i][i])*x[j];
		    	x[i] = c[i];
			}
		}
		errMax = 0.0;
		for(int i=0;i<N;i++) {
		    double err = fabs(xold[i] - x[i]);
	    	if(err > errMax) errMax = err;
		}
		for(int i=0;i<N;i++)
	    	xold[i] = x[i];
//			std::cout << "iter = " << iter << " : error = " << errMax << std::endl;
	}
    delete[] c;
    delete[] xold;

}

void SOR(int N,double **A, double *b, double *u,  double tol)
{
	double *uold, omega = 1.4;
	uold = new double[N];
	for(int i = 0; i<N;i++) uold[i] = u[i];
    RearrangeMatrix(N, A, b);
  int iter = 0;
  double errMax = 1.0;
	while(errMax > tol) {
		for (int i = 0; i < N; i++) {		
			double sum = 0.;		
			for (int j = 0; j < N; j++) {
				if (i != j) {
					sum += A[i][j] * u[j];
				}
			}
			u[i] = (b[i] - sum) / (A[i][i]);
			u[i] = (1.-omega)*uold[i] + omega * u[i];
		}	
		errMax =0.0;
		for(int i=0;i<N;i++) {
			double error = fabs(uold[i] - u[i])/fabs(uold[i]);
			if(error > errMax) errMax = error;
		}		
		for(int i=0; i<N;i++) {
		    uold[i] = u[i];
		}
		iter++;
	}
}

// DIRECT SOLVER : Gauss Elimination, LU Decomposition, TDMA

 void GaussElim(int n, double **A, double *b, double *x) {

//    int n = A.size();
    double **Ab;
    Ab = new double *[n];
    for(int i=0; i<n; i++) Ab[i] = new double[n+1];
    
    for(int i=0;i<n;i++) {
    	for(int j=0;j<n;j++) {
    		Ab[i][j] = A[i][j];
    	}
    	Ab[i][n] = b[i];
    }
    		
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(Ab[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(Ab[k][i]) > maxEl) {
                maxEl = abs(Ab[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = Ab[maxRow][k];
            Ab[maxRow][k] = Ab[i][k];
            Ab[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -Ab[k][i]/Ab[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    Ab[k][j] = 0;
                } else {
                    Ab[k][j] += c * Ab[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
//    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = Ab[i][n]/Ab[i][i];
        for (int k=i-1;k>=0; k--) {
            Ab[k][n] -= Ab[k][i] * x[i];
        }
    }
//    return x;
}

//-------------------------------------


void LUdcmp(int N, double **a,int *indx)
{   
    double TINY = 1.0e-20;
    double *vv;
    vv = new double[N];  // vv stores the implicit scaling of each row
	double d=1.0;
    
	for (int i=0;i < N;i++) {
		double big=0.0;
		for (int j=0;j<N;j++)
	    	if (fabs(a[i][j]) > big) big=fabs(a[i][j]);
		if (big == 0.0) std::cout << "Singular matrix in routine ludcmp"<< std::endl;
		vv[i]=1.0/big;  // save the scaling
    }
	int imax=0;
    for (int j=0;j<N;j++) {
		for (int i=0;i<j;i++) {
	    	double sum=a[i][j];
	    	for (int k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	    	a[i][j]=sum;
		}
		double big=0.0;
		for (int i=j;i<N;i++) {
	    	double sum=a[i][j];
	    	for (int k=0;k<j;k++) sum -= a[i][k]*a[k][j];
	    	a[i][j]=sum;
	    	double dum = vv[i]*fabs(sum);
	    	if (dum >= big) {
				big=dum;
				imax=i;
	    	}
		}
		if (j != imax) {
	    	for (int k=0;k<N;k++) {
				double dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
	    	}
	    	d = -d;
	    	vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != N-1) {
	    	double dum=1.0/(a[j][j]);
	    	for (int i=j+1;i<N;i++) a[i][j] *= dum;
		}
    }
    delete[] vv;
}

void LUbksb(int N, double  **a, int *indx, double *b)
{
    int ii = -1;
   
    for (int i=0;i<N;i++) {
		int ip=indx[i];
		double sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
	    	for (int j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.) ii=i;
		b[i]=sum;
    }
    for (int i=N-1;i>=0;i--) {
		double sum=b[i];
		for (int j=i+1;j<N;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void LU(int N ,double **A, double *b, double *x)
{
	int *indx;
    indx = new int[N];
    LUdcmp(N, A, indx);
    LUbksb(N, A, indx,b);
    for(int i=0;i<N;i++) x[i] = b[i];	
    delete[] indx;
 //   delete[] xold;

}


void TDMA (int N,double **A, double *b, double *u) {
	int i;
	double *P, *Q, denom;
    P = new double[N];
	Q = new double[N];

	P[0] = -A[0][1]/A[0][0];
	Q[0] = b[0] / A[0][0];
	for(int i=1;i<N-1;i++) {
		denom = (A[i][i] + A[i][i-1] *P[i-1]);
		P[i] = -A[i][i+1]/denom;
		Q[i] = (b[i] - A[i][i-1]*Q[i-1])/denom;
	}
	u[N-1]= (b[N-1] - A[N-1][N-2] * Q[N-2])/(A[N-1][N-1] + A[N-1][N-2]*P[N-2]);

	for(int i=N-2;i>=0;i--) {
		u[i] = P[i] * u[i+1] + Q[i];
	}
}
