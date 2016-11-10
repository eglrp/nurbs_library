//
// Created by miaodi on 11/8/16.
//

#ifndef UNTITLED_NURBSCURVE_H
#define UNTITLED_NURBSCURVE_H


#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <iostream>
using namespace Eigen;

template<class T, int N>
class NurbsCurve {
public:
    //constructor
    NurbsCurve();

    NurbsCurve(const NurbsCurve<T,N>& nurb);

    NurbsCurve(const Matrix<Matrix<T, N, 1>, Dynamic, 1> &P1, const Matrix<T, Dynamic, 1> &U1, int deg = 3);

    const Matrix<Matrix<T, N, 1>, Dynamic, 1> &getControlPoints() const { return P; }

    const Matrix<T, Dynamic, 1> &getKnotVector() const { return U; }

    const Matrix<T, N, 1> getControlPoint(int i) const { return P(i); }

    int getDegree() const { return deg_; }

    int getDOF() const { return dof; }

    const Matrix<T, N, 1> operator()(T u) const;

    void resetdof() { dof = P.rows(); }

    void resize(int n, int Deg);

    int findSpan(T u) const;

    void dersBasisFuns(int n, T u, int span, Matrix<T, Dynamic, Dynamic> &ders) const;

    void degreeElevate(int t);

protected:
    Matrix<Matrix<T, N, 1>, Dynamic, 1> P;  //the vector of control points
    Matrix<T, Dynamic, 1> U;  //knot vector
    int deg_; //degree of spline curve
private:
    int dof; //degree of freedom
};

template<class T, int N>
NurbsCurve<T, N>::NurbsCurve(): P(1), U(1), deg_(0) {
    resetdof();
}

template<class T, int N>
NurbsCurve<T,N>::NurbsCurve(const NurbsCurve &nurb): P(nurb.P),U(nurb.U),deg_(nurb.deg_){
    resetdof();
}
template<class T, int N>
NurbsCurve<T, N>::NurbsCurve(const Matrix<Matrix<T, N, 1>, Dynamic, 1> &P1, const Matrix<T, Dynamic, 1> &U1,
                             int degree): P(P1), U(U1), deg_(degree) {
    resetdof();
}

template<class T, int N>
int NurbsCurve<T, N>::findSpan(T u) const {
    if (u >= U[dof])
        return dof - 1;
    if (u <= U[deg_])
        return deg_;

    int low = 0;
    int high = dof + 1;
    int mid = (low + high) / 2;

    while (u < U[mid] || u >= U[mid + 1]) {
        if (u < U[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

template<class T, int N>
void NurbsCurve<T, N>::dersBasisFuns(int n, T u, int span, Matrix<T, Dynamic, Dynamic> &ders) const {
    T *left = new T[2 * (deg_ + 1)];
    T *right = &left[deg_ + 1];

    Matrix<T, Dynamic, Dynamic> ndu(deg_ + 1, deg_ + 1);
    T saved, temp;
    int j, r;

    ders.resize(n + 1, deg_ + 1);

    ndu(0, 0) = 1.0;
    for (j = 1; j <= deg_; j++) {
        left[j] = u - U[span + 1 - j];
        right[j] = U[span + j] - u;
        saved = 0.0;

        for (r = 0; r < j; r++) {
            // Lower triangle
            ndu(j, r) = right[r + 1] + left[j - r];
            temp = ndu(r, j - 1) / ndu(j, r);
            // Upper triangle
            ndu(r, j) = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        ndu(j, j) = saved;
    }

    for (j = deg_; j >= 0; --j)
        ders(0, j) = ndu(j, deg_);

    // Compute the derivatives
    Matrix<T, Dynamic, Dynamic> a(deg_ + 1, deg_ + 1);
    for (r = 0; r <= deg_; r++) {
        int s1, s2;
        s1 = 0;
        s2 = 1; // alternate rows in array a
        a(0, 0) = 1.0;
        // Compute the kth derivative
        for (int k = 1; k <= n; k++) {
            T d;
            int rk, pk, j1, j2;
            d = 0.0;
            rk = r - k;
            pk = deg_ - k;

            if (r >= k) {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }

            if (rk >= -1) {
                j1 = 1;
            } else {
                j1 = -rk;
            }

            if (r - 1 <= pk) {
                j2 = k - 1;
            } else {
                j2 = deg_ - r;
            }

            for (j = j1; j <= j2; j++) {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if (r <= pk) {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }
            ders(k, r) = d;
            j = s1;
            s1 = s2;
            s2 = j; // Switch rows
        }
    }

    // Multiply through by the correct factors
    r = deg_;
    for (int k = 1; k <= n; k++) {
        for (j = deg_; j >= 0; --j)
            ders(k, j) *= r;
        r *= deg_ - k;
    }
    delete[] left;
}

template<class T, int N>
const Matrix<T, N, 1> NurbsCurve<T, N>::operator()(T u) const {
    Matrix<T, Dynamic, Dynamic> Nb;
    int span = findSpan(u);

    dersBasisFuns(0, u, span, Nb);

    Matrix<T, N, 1> p=Matrix<T, N, 1>::Zero();
    for (int i = deg_; i >= 0; --i) {
        p += Nb(0, i) * P(span - deg_ + i);
    }
    return p;
}

/*!
  \brief Resizes a NURBS curve

  Resizes a NURBS curve. The old values are lost and new ones
  have to be created.

  \param n  the new number of control points for the curve
  \param Deg  the new degree for the curve

  \author Philippe Lavoie
  \date 24 January 1997
*/
template<class T, int N>
void NurbsCurve<T, N>::resize(int n, int Deg) {
    deg_ = Deg;
    P.resize(n);
    U.resize(n + deg_ + 1);
}

template<class T>
void binomialCoef(Matrix<T, Dynamic, Dynamic> &Bin) {
    int n, k;
    // Setup the first line
    Bin(0, 0) = 1.0;
    for (k = static_cast<int>(Bin.cols()) - 1; k > 0; --k)
        Bin(0, k) = 0.0;
    // Setup the other lines
    for (n = 0; n < static_cast<int>(Bin.rows()) - 1; n++) {
        Bin(n + 1, 0) = 1.0;
        for (k = 1; k < static_cast<int>(Bin.cols()); k++)
            if (n + 1 < k)
                Bin(n, k) = 0.0;
            else
                Bin(n + 1, k) = Bin(n, k) + Bin(n, k - 1);
    }
}

template<class T, int N>
void NurbsCurve<T, N>::degreeElevate(int t) {
    if (t <= 0) {
        return;
    }

    NurbsCurve<T, N> c(*this);
    int i, j, k;

    int n = dof-1;
    int p = c.deg_;
    int m = n + p + 1;
    int ph = p + t;
    int ph2 = ph / 2;
    Matrix<T, Dynamic, Dynamic> bezalfs=Matrix<T,Dynamic,Dynamic>::Zero(p + t + 1, p + 1); // coefficients for degree elevating the Bezier segment
    Matrix<Matrix<T, N, 1>, Dynamic, 1> bpts(p + 1); // pth-degree Bezier control points of the current segment
    Matrix<Matrix<T, N, 1>, Dynamic, 1> ebpts(p + t + 1); // (p+t)th-degree Bezier control points of the  current segment
    Matrix<Matrix<T, N, 1>, Dynamic, 1> Nextbpts(p - 1); // leftmost control points of the next Bezier segment
    Matrix<T, Dynamic, 1> alphas=Matrix<T, Dynamic, 1>::Zero(p - 1); // knot instertion alphas.

    // Compute the binomial coefficients
    Matrix<T, Dynamic, Dynamic> Bin=Matrix<T,Dynamic,Dynamic>::Zero(ph + 1, ph2 + 1);
    binomialCoef(Bin);

    // Compute Bezier degree elevation coefficients
    T inv, mpi;
    bezalfs(0, 0) = bezalfs(ph, p) = 1.0;
    for (i = 1; i <= ph2; i++) {
        inv = 1.0 / Bin(ph, i);
        mpi = std::min(p, i);
        for (j = std::max(0, i - t); j <= mpi; j++) {
            bezalfs(i, j) = inv * Bin(p, j) * Bin(t, i - j);
        }
    }

    for (i = ph2 + 1; i < ph; i++) {
        mpi = std::min(p, i);
        for (j = std::max(0, i - t); j <= mpi; j++)
            bezalfs(i, j) = bezalfs(ph - i, p - j);
    }

    resize(static_cast<int>(c.getControlPoints().cols() + c.getControlPoints().cols() * t), ph); // Allocate more control points than necessary

    int mh = ph;
    int kind = ph + 1;
    T ua = c.U(0);
    T ub = 0.0;
    int r = -1;
    int oldr;
    int a = p;
    int b = p + 1;
    int cind = 1;
    int rbz, lbz = 1;
    int mul, save, s;
    T alf;
    int first, last, kj;
    T den, bet, gam, numer;

    P(0) = c.P(0);
    for (i = 0; i <= ph; i++) {
        U(i) = ua;
    }

    // Initialize the first Bezier segment
    for (i = 0; i <= p; i++)
        bpts(i) = c.P(i);
    while (b < m) { // Big loop thru knot vector
        i = b;
        while (b < m && c.U(b) == c.U(b + 1)) // for some odd reasons... == doesn't work
            b++;
        mul = b - i + 1;
        mh += mul + t;
        ub = c.U(b);
        oldr = r;
        r = p - mul;
        if (oldr > 0)
            lbz = (oldr + 2) / 2;
        else
            lbz = 1;
        if (r > 0)
            rbz = ph - (r + 1) / 2;
        else
            rbz = ph;
        if (r > 0) { // Insert knot to get Bezier segment
            numer = ub - ua;
            for (k = p; k > mul; k--) {
                alphas(k - mul - 1) = numer / (c.U(a + k) - ua);
            }
            for (j = 1; j <= r; j++) {
                save = r - j;
                s = mul + j;
                for (k = p; k >= s; k--) {
                    bpts(k) = alphas(k - s) * bpts(k) + (1.0 - alphas(k - s)) * bpts(k - 1);
                }
                Nextbpts(save) = bpts(p);
            }
        }

        for (i = lbz; i <= ph; i++) { // Degree elevate Bezier,  only the points lbz,...,ph are used
            ebpts(i) = Matrix<T,Dynamic,1>::Zero(N);
            mpi = std::min(p, i);
            for (j = std::max(0, i - t); j <= mpi; j++)
                ebpts(i) += bezalfs(i, j) * bpts(j);
        }

        if (oldr > 1) { // Must remove knot u=c.U[a] oldr times
            // if(oldr>2) // Alphas on the right do not change
            //	alfj = (ua-U[kind-1])/(ub-U[kind-1]) ;
            first = kind - 2;
            last = kind;
            den = ub - ua;
            bet = (ub - U(kind - 1)) / den;
            for (int tr = 1; tr < oldr; tr++) { // Knot removal loop
                i = first;
                j = last;
                kj = j - kind + 1;
                while (j - i > tr) { // Loop and compute the new control points for one removal step
                    if (i < cind) {
                        alf = (ub - U(i)) / (ua - U(i));
                        P(i) = alf * P(i) + (1.0 - alf) * P(i - 1);
                    }
                    if (j >= lbz) {
                        if (j - tr <= kind - ph + oldr) {
                            gam = (ub - U(j - tr)) / den;
                            ebpts(kj) = gam * ebpts(kj) + (1.0 - gam) * ebpts(kj + 1);
                        } else {
                            ebpts(kj) = bet * ebpts(kj) + (1.0 - bet) * ebpts(kj + 1);
                        }
                    }
                    ++i;
                    --j;
                    --kj;
                }
                --first;
                ++last;
            }
        }

        if (a != p) // load the knot u=c.U[a]
            for (i = 0; i < ph - oldr; i++) {
                U(kind++) = ua;
            }
        for (j = lbz; j <= rbz; j++) { // load control points onto the curve
            P(cind++) = ebpts(j);
        }

        if (b < m) { // Set up for next pass thru loop
            for (j = 0; j < r; j++)
                bpts(j) = Nextbpts(j);
            for (j = r; j <= p; j++)
                bpts(j) = c.P(b - p + j);
            a = b;
            b++;
            ua = ub;
        } else {
            for (i = 0; i <= ph; i++)
                U(kind + i) = ub;
        }
    }
    resize(mh - ph, ph); // Resize to the proper number of control points

}




#endif //UNTITLED_NURBSCURVE_H
