#include <iostream>
#include "NurbsCurve.h"
using namespace std;
int main() {
    NurbsCurve<double, 2> a;
    cout<<a.getDegree();
    cout<<endl;
    return 0;
}