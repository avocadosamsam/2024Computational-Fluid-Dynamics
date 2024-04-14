#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double f1(double x, double y){
    double del_y{};
    del_y = (2*sin(3*x)-(pow(x,2)*pow(y,2)))*exp(-y);

    return del_y;
}

int main(){
    double t{}, tEnd{1.0}, step{0.01};
    double *Y;
    double k1{}, k2{}, k3{}, k4{};
    int N;

    ofstream fout;
    fout.open("problem_14.txt");

    N = int(tEnd/step);
    Y = new double[N];

    Y[0] = 5;

    for(int i=0; i<N; i++){
        k1 = step*f1(t, Y[i]);
        //cout << setw(10) << k1 << endl;
        k2 = step*f1(t+(step/2), Y[i]+(k1/2));
        k3 = step*f1(t+(step/2), Y[i]+(k2/2));
        k4 = step*f1(t+step, Y[i]+k3);
        Y[i+1] = Y[i] + (k1+2*k2+2*k3+k4)/6;

        t+=step;

        cout << setw(10) << t << setw(15) << Y[i+1] << endl;
        fout << setw(10) << t << setw(15) << Y[i+1] << endl;

    }

    return 0;
e