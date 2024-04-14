// double Dfunc(double t, double y){
//     double dydt;
//     dydt='given function of (t,y)';
//     return dydt;
// }
// y[0]='given value';
// for (int i=0; i<n;i++){
//     k1=delt*Dfunc(t,y[i]);
//     k2=delt*Dfunc(t+delt, y[i]+k1);
//     y[i+1]=y[i]+0.5*(k1+k2);
// }


#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double f1(double t, double Vr){
    double me{0.1}, ve{100},g{9.81}; //me:엔진에서 발사되는 질량유량, ve:분사속도
    double mr{}, m0{1.0}, delVr{};
    mr=m0-me*t;
    delVr=(me*ve)/mr-g;

    return delVr;
}

int main(){
    double t{}, tEnd{1.0}, step_time{0.01};
    double *Vr_Eu, *Vr_RK2, *Vr_RK4;
    double *zr_1, *zr_2, *zr_3;
    double k1,k2,k3,k4; //RK4
    //double k1_2, k2_2; //RK2
    
    int N;

    ofstream fout;
    fout.open("rocket_Vr_Zr.txt");

    N = int(tEnd/step_time);

    //Vr_Eu = new double[N];
    //Vr_RK2 = new double[N];
    Vr_RK4 = new double[N];
    zr_1 = new double[N];
    zr_2 = new double[N];
    zr_3 = new double[N];
    

    //zr = new double[N];
    //Vr_Eu[0] = 0.0;
    //Vr_RK2[0] = 0.0;
    Vr_RK4[0] = 0.0;
    zr_1[0] = 0.0;
    zr_2[0] = 0.0;
    zr_3[0] = 0.0;

    for (int i=0; i<N; i++){
        k1 = step_time*f1(t,Vr_RK4[i]);
        k2 = step_time*f1(t+(step_time/2), Vr_RK4[i]+(k1/2));
        k3 = step_time*f1(t+(step_time/2), Vr_RK4[i]+(k2/2));
        k4 = step_time*f1(t+step_time, Vr_RK4[i]+k3);
        Vr_RK4[i+1] = Vr_RK4[i] + (k1+2*k2+2*k3+k4)/6;

        //Zr : 로켓의 위치
        zr_1[i+1] = zr_1[i] + step_time*Vr_RK4[i+1]; //방법 1
        zr_2[i+1] = zr_2[i] + (step_time/2)*(Vr_RK4[i]+Vr_RK4[i+1]); //방법2
        t += step_time;
    }

    t = 0;
    zr_3[1] = 0.0;
    for (int i=1; i<N+1; i++){
        //Euler
        //Vr_Eu[i+1] = Vr_Eu[i] + step_time*f1(t, Vr_Eu[i]);

        //RK2
        //k1_2 = step_time*f1(t, Vr_RK2[i]);
        //k2_2 = step_time*f1(t+step_time, Vr_RK2[i]+k1_2);
        //Vr_RK2[i+1] = Vr_RK2[i] + (k1_2 + k2_2)/2;

        //RK4, Vr : 로켓의 속도
        // k1 = step_time*f1(t,Vr_RK4[i]);
        // k2 = step_time*f1(t+(step_time/2), Vr_RK4[i]+(k1/2));
        // k3 = step_time*f1(t+(step_time/2), Vr_RK4[i]+(k2/2));
        // k4 = step_time*f1(t+step_time, Vr_RK4[i]+k3);
        // Vr_RK4[i+1] = Vr_RK4[i] + (k1+2*k2+2*k3+k4)/6;

        //Zr : 로켓의 위치
        zr_3[i+1] = 2*zr_3[i]-zr_3[i-1]+(step_time*step_time)*f1(t,Vr_RK4[i]);

        t += step_time;
        // Only Vr(Euler, RK2, RK4)
        //cout << setw(10) << t << setw(15) << Vr_Eu[i+1] << setw(15) << Vr_RK2[i+1] << setw(15) << Vr_RK4[i+1] << endl;
        //fout << setw(10) << t << setw(15) << Vr_Eu[i+1] << setw(15) << Vr_RK2[i+1] << setw(15) << Vr_RK4[i+1] << endl;

        // Vr(RK4) & Zr(1,2,3)
    }

    for(int i=0;i<N;i++){
        cout << setw(15) << Vr_RK4[i+1] << setw(15) << zr_1[i+1] << setw(15) << zr_2[i+1] << setw(15) << zr_3[i+1] << endl;
        fout << setw(15) << Vr_RK4[i+1] << setw(15) << zr_1[i+1] << setw(15) << zr_2[i+1] << setw(15) << zr_3[i+1] << endl;
    }

    return 0;
}