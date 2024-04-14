#include <iostream>
#include <cmath>

using namespace std;

// int func_distance(double side1, double side2){
//     return sqrt(side1*side1+side2*side2);
// }

// int main(){
//     double x1=1, y1=5, x2=4, y2=7;
//     double side1{}, side2{}, distance{};
//     side1= x2-x1;
//     side2= y2-y1;
//     distance = func_distance(side1, side2);
//     cout << "The distance between the two point is "<< distance << endl;

//     return 0;
// }

// void swap(double &x, double &y){
//     double tmp{};
//     tmp = x;
//     x = y;
//     y = tmp;
// }

// int main(){
//     double a{30}, b{45};
//     cout << "a= "<< a << "b= " << b << endl;

//     swap(a,b);
//     cout << "a= "<< a << "b= " << b << endl;
//     return 0;
// }


double FrictionFactor(double Re, double D, double epsilon){
    double f{}, denomiator{};

    if (Re <= 2400)
        f = 64/Re;
    else 
        denomiator = -1.8*log10((6.9/Re)+pow((epsilon/D)/3.7,1.11));
        f = pow((1/denomiator),2);

    return f;
}

int main(){
    double f{}, D{0.01}, epsilon{1.e-6};
    double Re{2000}, deltaRe{};

    deltaRe = (5000-2000)/100;
    for(int i=0; i<100; ++i){
        f= FrictionFactor(Re, D, epsilon);
        cout << "Re= "<< Re << " " << "f= " << f << endl;
        Re += deltaRe;
    }

    return 0;
}

