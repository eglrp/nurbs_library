//
// Created by miaodi on 11/8/16.
//

#ifndef UNTITLED_NURBSCURVE_H
#define UNTITLED_NURBSCURVE_H


#include <eigen3/Eigen/Dense>

using namespace Eigen;

template<class T, int N>
class NurbsCurve {
public:
    //constructor
    NurbsCurve();

    NurbsCurve(const Matrix<Matrix<T, N, 1>, Dynamic, 1> &P1, const Matrix<T, Dynamic, 1> &U1, int deg = 3);

    const Matrix<Matrix<T, N, 1>, Dynamic, 1> &getControlPoints() const { return P; }

    const Matrix<T, Dynamic, 1> &getKnotVector() const { return U; }

    const Matrix<T, N, 1> getControlPoint(int i) { return P(i); }

    int getDegree() const { return deg_; }

    int getDOF() const { return dof; }

    void resetdof() { dof = P.rows(); }

protected:
    Matrix<Matrix<T, N, 1>, Dynamic, 1> P;  //the vector of control points
    Matrix<T, Dynamic, 1> U;  //knot vector
    T deg_; //degree of spline curve
private:
    int dof; //degree of freedom
};

template<class T, int N>
NurbsCurve<T, N>::NurbsCurve(): P(1), U(1), deg_(0) {
    resetdof();
}

template<class T, int N>
NurbsCurve<T, N>::NurbsCurve(const Matrix<Matrix<T, N, 1>, Dynamic, 1> &P1, const Matrix<T, Dynamic, 1> &U1, int degree) {
    int nSize = P1.rows();
    int mSize = U1.rows();
    deg_ = degree;
    // Exception is needed
    P.resize(nSize);
    U.resize(mSize);
    P = P1;
    U = U1;

}

#endif //UNTITLED_NURBSCURVE_H
