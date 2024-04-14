#include <iostream>

using namespace std;

int main(){
    float *xp, *yp, *zp;
    float x{1.0}, y{2.0}, z(3.0);

    xp=&x;
    yp=&y;
    zp=&z;

    cout << "x = " << x  << ", *xp = " << *xp << endl;
    cout << "x = " << x  << ", *xp = " << *xp << endl;
    cout << "x = " << x  << ", *xp = " << *xp << endl;


    return 0;
}