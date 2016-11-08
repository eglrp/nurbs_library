#include <iostream>
#include "NurbsCurve.h"
using namespace std;
int main() {
    NurbsCurve<double, 2> a;
    cout<<endl;
    cout<<a.getDegree();
    return 0;
}