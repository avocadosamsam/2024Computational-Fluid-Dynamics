#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double f1(double t, double Vr){
    double me(0.1), ve(100), g(9.81);
    double mr,m0(1.0), delVr;
    mr= m0- me*t;
    delVr=(me*ve)/mr-g;
    
    return delVr;
}

int main(){
    double t(0), tEnd(1.0), dTime;
    double *Vr;
    double *zr;
    double k1,k2,k3,k4;
    int N;
    ofstream fout;
    fout.open("rocket.out");

    dTime=0.01;
    N = int(tEnd/dTime);
    
    Vr = new double[N]; //동적할당
    zr= new double[N];
    Vr[0] = 0.0; //initial condition

    for(int i=0;i<N;i++){
        //Vr[i+1] = Vr[i] + dTime * f1(t,Vr[i]);

        //룬지쿠타
        // k1=dTime*f1(t,Vr[i]);
        // k2=dTime*f1(t+dTime,Vr[i]+k1);
        // Vr[i+1]=Vr[i]+0.5*(k1+k2);

        // 룬자쿠타 4차
        k1=dTime*f1(t,Vr[i]);
        k2=dTime*f1(t+(dTime/2),Vr[i]+(k1/2));
        k3=dTime*f1(t+(dTime/2),Vr[i]+(k2/2));
        k4=dTime*f1(t+dTime,Vr[i]+k3);
        Vr[i+1]=Vr[i]+(1/6)*(k1+2*k2+2*k3+k4);

        //시간에 따른 로켓의 위치변화
        zr[i+1]=Vr[i]+dTime*(Vr[i+1]+Vr[i])/2;

        t += dTime;
        cout << setw(15) << t << setw(15) << Vr[i+1] << setw(15) << zr[i+1] << endl;
        fout << setw(15) << t << setw(1) << Vr[i+1] << setw(15) << zr[i+1] << endl;

    }
    

    return 0;
}