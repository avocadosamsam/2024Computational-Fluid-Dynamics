#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double f1(double t, double Vr){
    double me(0.1), ve(100), g(9.81);
    double mr, m0(1.0), delVr;
    mr=m0-me*t;
    delVr=(me*ve)/mr-g;

    return delVr;
}

int main(){

    return 0;
}