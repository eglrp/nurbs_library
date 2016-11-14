#include <iostream>
#include "NurbsCurve.h"
#include <map>
#include <vector>
#include <cmath>

#include "gnuplot-iostream.h"
using namespace std;
int main() {
    Vector2d point1(0, 0);
    Vector2d point2(2, 2);
    Vector2d point3(5, 5);
    Vector2d point4(3, 3);
    Vector2d point5(-3, 3);
    Vector2d point6(7, 7);
    Vector2d point7(9, 9);
    Vector2d point8(10, 10);

    VectorXd knotVector(11);
    knotVector << 0,0,0,1,2,4,4,5,6,6,6;
    Matrix<Matrix<double, 2, 1>, Dynamic, 1> points(8);
    points << point1, point2, point3, point4, point5, point6, point7, point8;
    NurbsCurve<double, 2> a(points, knotVector, 2);

    cout << a(4) << endl;
    NurbsCurve<double, 2> b(a);
    a.degreeElevate(3);
    cout << a(4) << endl;
    a.knotInsertion(3.5, 3);
    cout << a(4) << endl;
    cout << a.getKnotVector() << endl;
    VectorXd X(5);
    X << 2.4, 2.6, 3.3, 4.4, 4.4;
    a.refineKnotVectorCurve(X);
    cout << a.getKnotVector() << endl;
    cout << a(4) << endl;
    return 0;
}