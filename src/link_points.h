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

inline double logit_inv(double eta, int k) {
    double p = 1.0 / (1.0 + std::exp(-eta));
    return logistic_poly(p, k);
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

inline double rhobit_inv(double eta, int k) {
    double t = std::tanh(eta), t2 = t * t;
    switch (k) {
    case 1: return 1.0 - t2;
    case 2: return -2.0 * t * (1.0 - t2);
    case 3: return -2.0 + 8.0 * t2 - 6.0 * t2 * t2;
    case 4: return t * (16.0 + t2 * (-40.0 + 24.0 * t2));
    case 5: return 16.0 + t2 * (-136.0 + t2 * (240.0 - 120.0 * t2));
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
