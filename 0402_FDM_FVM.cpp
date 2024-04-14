#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "LAESolver.h"

using namespace std;



int main(){
    double *x, *u, **A, *B;
    double dX(0.0), Length(1.0), PI(3.141592);
    int N; //Num of slice

    ofstream fout;
    fout.open("FDM_FVM(0.05).txt");

    //Modeling
    dX=0.05;
    N=(Length/dX) + 1;
    x = new double[N];
    x[0]=0.0;
    for(int i=0; i<N-1; i++){
        x[i+1]=x[i]+dX;
        // cout << x[i+1] << endl;
    }

    //Matrix, Allocation
    u = new double[N];
    B = new double [N];
    A = new double* [N];
    for (int i=0; i<N; i++){
        A[i] = new double[N];
    }

    //Set A and B matrix
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            A[i][j] = 0.0;
        }
        B[i]=0.0;
        u[i]=0.0;
    }

    for (int i=1; i<N-1; i++){
        A[i][i] = -2.0;
        A[i][i+1] = 1.0;
        A[i][i-1] = 1.0;
        B[i] = -sin(2*PI*x[i])*dX*dX;
    }

    //Boundary Condition
    A[0][0] = 1.0; B[0] = 0.0;
    A[N-1][N-1] = 1.0; B[N-1] = 0.0;
    
    //Solver for U
    GaussElim(N,A,B,u);
    double Uexact(0); 

    for (int i=0; i<N; i++){
        Uexact = sin(2*PI*x[i])/(4*PI*PI); // theory
        cout << setw(15) << x[i] << setw(15) << u[i] << setw(15) << Uexact << endl;
        fout << setw(15) << x[i] << setw(15) << u[i] << setw(15) << Uexact << endl;
    }

    return 0;
}