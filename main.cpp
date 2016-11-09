#include <iostream>
#include "NurbsCurve.h"
using namespace std;
int main() {
    Vector2d point1(0, 1);
    Vector2d point2(2, 4);
    Vector2d point3(5, 5);
    Vector2d point4(3, 2);

    VectorXd knotVector(7);
    knotVector << 0, 0, 0, .5, 1, 1, 1;
    Matrix<Matrix<double, 2, 1>, Dynamic, 1> points(4);
    points << point1, point2, point3, point4;
    NurbsCurve<double, 2> a(points, knotVector, 2);
    cout << a.getDOF();
    return 0;
}