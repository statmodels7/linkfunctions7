#ifndef LF7_LINK_POINTS_H
#define LF7_LINK_POINTS_H

#include <Rcpp.h>
#include <cfloat>
#include <cmath>

// The inverse links' derivatives at ONE point, written once: the vector
// kernels of link_kernels.cpp loop over them, and the scalar C entries of
// lf7_ccallable.cpp call them, so the two routes run the same expressions.

namespace lf7 {

// the logistic derivative polynomials in p: shared by the logit (as is),
// the doubly bounded link (scaled by the width) and the softplus (an
// antiderivative of the logistic, so shifted one place down)
inline double logistic_poly(double p, int k) {
    double pq = p * (1.0 - p);
    switch (k) {
    case 1: return pq;
    case 2: return pq * (1.0 - 2.0 * p);
    case 3: return pq * (1.0 + p * (-6.0 + 6.0 * p));
    case 4: return pq * (1.0 + p * (-14.0 + p * (36.0 - 24.0 * p)));
    // P_{k+1} = (1 - 2p) P_k + pq P_k', so the fifth follows from the fourth
    // and the orders above it are generated rather than transcribed.
    case 5: return pq * (1.0 + p * (-30.0 + p * (150.0 +
                                    p * (-240.0 + 120.0 * p))));
    default: return NA_REAL;
    }
}

// The same polynomials written in p, q = 1 - p and pq, with q computed from
// eta and not as 1 - p: at eta = 37 the subtraction returns q = 0 and every
// order 0, where pq is 8.5e-17. P_3 = pq(1 - 6pq), P_4 = pq(q - p)(1 - 12pq)
// and P_5 = pq(1 - 30pq + 120(pq)^2) expand to logistic_poly()'s forms.
inline double logit_inv(double eta, int k) {
    double p = 1.0 / (1.0 + std::exp(-eta));
    double q = 1.0 / (1.0 + std::exp(eta));
    double pq = p * q;
    switch (k) {
    case 1: return pq;
    case 2: return pq * (q - p);
    case 3: return pq * (1.0 - 6.0 * pq);
    case 4: return pq * (q - p) * (1.0 - 12.0 * pq);
    case 5: return pq * (1.0 + pq * (-30.0 + 120.0 * pq));
    default: return NA_REAL;
    }
}

inline double probit_inv(double e, int k) {
    double phi = R::dnorm4(e, 0.0, 1.0, 0);
    switch (k) {
    case 1: return phi;
    case 2: return -e * phi;
    case 3: return (e * e - 1.0) * phi;
    case 4: return (3.0 * e - e * e * e) * phi;
    case 5: {
        double e2 = e * e;
        return (e2 * e2 - 6.0 * e2 + 3.0) * phi;
    }
    default: return NA_REAL;
    }
}

inline double cloglog_inv(double eta, int k) {
    double w = std::exp(eta);
    double E = std::exp(eta - w);
    switch (k) {
    case 1: return E;
    case 2: return E * (1.0 - w);
    case 3: return E * (1.0 + w * (-3.0 + w));
    case 4: return E * (1.0 + w * (-7.0 + w * (6.0 - w)));
    case 5: return E * (1.0 + w * (-15.0 + w * (25.0 + w * (-10.0 + w))));
    default: return NA_REAL;
    }
}

inline double loglog_inv(double eta, int k) {
    double z = std::exp(-eta);
    double E = std::exp(-z);
    switch (k) {
    case 1: return E * z;
    case 2: return E * z * (z - 1.0);
    case 3: return E * z * (1.0 + z * (-3.0 + z));
    case 4: return E * z * (-1.0 + z * (7.0 + z * (-6.0 + z)));
    case 5: return E * z * (1.0 + z * (-15.0 + z * (25.0 + z * (-10.0 + z))));
    default: return NA_REAL;
    }
}

inline double cauchit_inv(double e, int k) {
    double u = 1.0 + e * e;
    switch (k) {
    case 1: return 1.0 / (M_PI * u);
    case 2: return -2.0 * e / (M_PI * u * u);
    case 3: return 2.0 * (3.0 * e * e - 1.0) / (M_PI * u * u * u);
    case 4: return 24.0 * e * (1.0 - e * e) / (M_PI * u * u * u * u);
    case 5: {
        double u2 = u * u, e2 = e * e;
        return 24.0 * (1.0 - 10.0 * e2 + 5.0 * e2 * e2) / (M_PI * u2 * u2 * u);
    }
    default: return NA_REAL;
    }
}

// Written in t = tanh(eta) and u = sech^2(eta) = 1 - t^2, with u computed
// from eta as 4a / (1 + a)^2, a = exp(-2|eta|), and not as 1 - t^2: at
// eta = 20 the subtraction returns 0 where u is 1.7e-17. Each order is u
// times a polynomial in t and u, from dt = u and du = -2tu.
inline double rhobit_inv(double eta, int k) {
    double t = std::tanh(eta), t2 = t * t;
    double a = std::exp(-2.0 * std::fabs(eta));
    double u = 4.0 * a / ((1.0 + a) * (1.0 + a));
    switch (k) {
    case 1: return u;
    case 2: return -2.0 * t * u;
    case 3: return u * (4.0 * t2 - 2.0 * u);
    case 4: return 8.0 * t * u * (2.0 * u - t2);
    case 5: return u * (16.0 * u * u - 88.0 * t2 * u + 16.0 * t2 * t2);
    default: return NA_REAL;
    }
}

// exp_floored() of link_class.R: pmax(exp(eta), exp_floor) with
// exp_floor = (24 / double.xmax)^0.25, written as the R constant computes it
inline double exp_floored(double eta) {
    static const double floor_ = std::pow(24.0 / DBL_MAX, 0.25);
    double v = std::exp(eta);
    if (!ISNAN(v) && v < floor_) v = floor_;
    return v;
}

// R's `^`, which the R methods use: x * x at an exponent of two, R_pow()
// otherwise
inline double rpow(double x, double y) {
    return y == 2.0 ? x * x : R_pow(x, y);
}

} // namespace lf7

#endif
