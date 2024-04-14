#include <iostream>
#include <cmath>
#include "LAESolver.h"
//#include "Lasol2.h"
using namespace std;

int main() {
   int N(3);
   double **A, *x, *b;

   x = new double[N];
   b = new double[N];
   A = new double *[N];
   for(int i = 0 ;i <N; ++i) A[i] = new double [N];
    
   for(int i = 0;i<N; i++) x[i] = 0.0;

/* EX1
   A[0][0] = 1.; A[0][1] = 1.; A[0][2] = -1.;
   A[1][0] = 2.0; A[1][1] = -1.; A[1][2] = 3.; 
   A[2][0] = 1; A[2][1] = 2; A[2][2] = 1.0;

   b[0] = 0; b[1] = 9 ; b[2] = 8;
*/
/* EX2
   A[0][0] = 2.; A[0][1] = 1.0; A[0][2] = -1.;
   A[1][0] = 1.0; A[1][1] = -1.; A[1][2] = 3.; 
   A[2][0] = 3.; A[2][1] = 2.; A[2][2] = 1.0;

   b[0] = 8.0; b[1] = -4 ; b[2] = 11.0;
*/

// EX3
   A[0][0] = 3.0; A[0][1] = 2.0; A[0][2] = -1;
   A[1][0] = -1.0; A[1][1] = 3; A[1][2] = 2; 
   A[2][0] = 1; A[2][1] = -1; A[2][2] = -1;

   b[0] = 10.0; b[1] = 5 ; b[2] = -1.0;


   GaussElim(N, A, b, x);

    cout << " X = ( " ;
   for(int i = 0; i<N ; ++i) {
       cout << " " << x[i] << ", " ;
   }
   cout << ") " <<endl;
         






    return 0;
}
