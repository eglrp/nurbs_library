//
// Created by miaodi on 11/8/16.
//

#ifndef NURBSCURVE_H
#define NURBSCURVE_H


#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <iostream>

using namespace Eigen;

template<class T, int N>
class NurbsCurve {
public:
    //constructor
    NurbsCurve();

    NurbsCurve(const NurbsCurve<T, N> &nurb);

    NurbsCurve(const Matrix<Matrix<T, N, 1>, Dynamic, 1> &P1, const Matrix<T, Dynamic, 1> &U1, int deg = 3);

    inline const Matrix<Matrix<T, N, 1>, Dynamic, 1> &getControlPoints() const { return P; }

    inline const Matrix<T, Dynamic, 1> &getKnotVector() const { return U; }

    inline const Matrix<T, N, 1> getControlPoint(int i) const { return P(i); }

    inline int getDegree() const { return deg_; }

    inline int getDOF() const { return dof; }

    const Matrix<T, N, 1> operator()(T u) const;

    inline void resetdof() { dof = static_cast<int>(P.rows()); }

    void resize(int n, int Deg);

    int findSpan(T u) const;

    void dersBasisFuns(int n, T u, int span, Matrix<T, Dynamic, Dynamic> &ders) const;

    void degreeElevate(int t);

    void knotInsertion(T u, int r);

    void refineKnotVectorCurve(Matrix<T, Dynamic, 1> X);

protected:
    Matrix<Matrix<T, N, 1>, Dynamic, 1> P;  //the vector of control points
    Matrix<T, Dynamic, 1> U;  //knot vector
    int deg_; //degree of spline curve
private:
    int dof; //degree of freedom
};


#endif //NURBSCURVE_H
