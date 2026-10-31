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
    // P_6 = pq (1 - 2p) (1 - 60 pq + 360 (pq)^2), from d/deta of the fifth
    case 6: return pq * (1.0 - 2.0 * p) * (1.0 + pq * (-60.0 + 360.0 * pq));
    default: return NA_REAL;
    }
}

// The same polynomials written in p, q = 1 - p and pq, with q computed from
// eta and not as 1 - p: at eta = 37 the subtraction returns q = 0 and every
// order 0, where pq is 8.5e-17. P_3 = pq(1 - 6pq), P_4 = pq(q - p)(1 - 12pq)
// and P_5 = pq(1 - 30pq + 120(pq)^2) expand to logistic_poly()'s forms.
inline double logit_inv(double eta, int k) {
    // Beyond |eta| = 600 the value is exp(-|eta|) to the working precision
    // (pq = a / (1 + a)^2 with a < 1e-261), signed by the parity of the order
    // and the side. The direct form returns 0 from |eta| = 709.78, where
    // exp(|eta|) overflows, although the value is a number down to 745.
    if (std::fabs(eta) > 600.0) {
        double a = std::exp(-std::fabs(eta));
        switch (k) {
        case 1: case 3: case 5: return a;
        case 2: case 4: case 6: return eta > 0.0 ? -a : a;
        default: return NA_REAL;
        }
    }
    double p = 1.0 / (1.0 + std::exp(-eta));
    double q = 1.0 / (1.0 + std::exp(eta));
    double pq = p * q;
    switch (k) {
    case 1: return pq;
    case 2: return pq * (q - p);
    case 3: return pq * (1.0 - 6.0 * pq);
    case 4: return pq * (q - p) * (1.0 - 12.0 * pq);
    case 5: return pq * (1.0 + pq * (-30.0 + 120.0 * pq));
    case 6: return pq * (q - p) * (1.0 + pq * (-60.0 + 360.0 * pq));
    default: return NA_REAL;
    }
}

// s * logit_inv(eta, k) where the product, not the factor, is the quantity
// that must not underflow: beyond |eta| = 600 it is formed as
// exp(log|s| - |eta|), so that a large s (the width of a doubly bounded link,
// a power of the softplus scale) keeps the representable values.
inline double logit_inv_scaled(double eta, int k, double s) {
    if (std::fabs(eta) > 600.0) {
        if (ISNAN(eta) || ISNAN(s) || s == 0.0) return eta * s;
        if (k < 1 || k > 6) return NA_REAL;
        // the sign is the parity of the order and the side of eta; the value
        // itself may already be the subnormal 0 that this branch exists for
        double sg = (k % 2 == 1) ? 1.0 : (eta > 0.0 ? -1.0 : 1.0);
        if (s < 0.0) sg = -sg;
        return sg * std::exp(std::log(std::fabs(s)) - std::fabs(eta));
    }
    return s * logit_inv(eta, k);
}

// phi(e) P_k(e) with P_k the Hermite polynomials: -e phi, (e^2 - 1) phi, ...
// Where phi underflows (|e| > 37.1) the product is formed through
// exp_times(), so that the values that are still numbers up to |e| = 38.6
// are not flushed to zero or computed from a subnormal phi; beyond 1e30
// phi is 0 and e^5 would overflow, and the value is 0.
inline double probit_inv(double e, int k) {
    double phi = R::dnorm4(e, 0.0, 1.0, 0);
    double e2 = e * e, P;
    switch (k) {
    case 1: P = 1.0; break;
    case 2: P = -e; break;
    case 3: P = e2 - 1.0; break;
    case 4: P = 3.0 * e - e * e2; break;
    case 5: P = e2 * e2 - 6.0 * e2 + 3.0; break;
    case 6: P = -e * (e2 * e2 - 10.0 * e2 + 15.0); break;
    default: return NA_REAL;
    }
    if (std::fabs(e) > 1e30) return 0.0;
    if (phi > 1e-290 || ISNAN(phi) || ISNAN(P) || P == 0.0) return phi * P;
    return std::copysign(std::exp(-0.5 * e2 - 0.91893853320467274178 + std::log(std::fabs(P))), P);
}

// exp(a) * P without a premature underflow: where exp(a) is below 1e-290 the
// product is formed as exp(a + log|P|), so that a result still representable
// (a subnormal, say) is not flushed to zero by the first factor. P is a
// finite polynomial value there; a NaN in either factor is returned as is.
inline double exp_times(double a, double P) {
    double E = std::exp(a);
    if (E > 1e-290 || ISNAN(E) || ISNAN(P) || P == 0.0) return E * P;
    return std::copysign(std::exp(a + std::log(std::fabs(P))), P);
}

// h(eta) = 1 - exp(-exp(eta)). With w = exp(eta), the derivative of order k is
// exp(eta - w) Q_k(w), Q_1 = 1, Q_{k+1} = (1 - w) Q_k + w Q_k'. For w > 1e30
// exp(-w) is 0 and the polynomial would overflow (w^5 from w = 1e61.6, so
// 0 * Inf = NaN, and Inf for w = Inf at eta > 709.78): the value is 0 there.
inline double cloglog_inv(double eta, int k) {
    double w = std::exp(eta);
    double P;
    switch (k) {
    case 1: P = 1.0; break;
    case 2: P = 1.0 - w; break;
    case 3: P = 1.0 + w * (-3.0 + w); break;
    case 4: P = 1.0 + w * (-7.0 + w * (6.0 - w)); break;
    case 5: P = 1.0 + w * (-15.0 + w * (25.0 + w * (-10.0 + w))); break;
    case 6: P = 1.0 + w * (-31.0 + w * (90.0 + w * (-65.0 + w * (15.0 - w)))); break;
    default: return NA_REAL;
    }
    if (w > 1e30) return 0.0;
    return exp_times(eta - w, P);
}

// h(eta) = exp(-exp(-eta)). With z = exp(-eta), the derivative of order k is
// (-1)^{k-1} exp(-z) z Q_k(z), the polynomials of the cloglog; exp(-z) z is
// formed as exp(-eta - z). Zero for z > 1e30, as in the cloglog.
inline double loglog_inv(double eta, int k) {
    double z = std::exp(-eta);
    double P;
    switch (k) {
    case 1: P = 1.0; break;
    case 2: P = z - 1.0; break;
    case 3: P = 1.0 + z * (-3.0 + z); break;
    case 4: P = -1.0 + z * (7.0 + z * (-6.0 + z)); break;
    case 5: P = 1.0 + z * (-15.0 + z * (25.0 + z * (-10.0 + z))); break;
    case 6: P = -(1.0 + z * (-31.0 + z * (90.0 + z * (-65.0 + z * (15.0 - z))))); break;
    default: return NA_REAL;
    }
    if (z > 1e30) return 0.0;
    return exp_times(-eta - z, P);
}

// h(eta) = 1/2 + atan(eta)/pi. With a = 1/(1 + e^2) and b = e^2 a = 1 - a,
// every order is a power of a times a polynomial in (a, b) and the factor
// e a, which stay finite and exact for any |e| (u = 1 + e^2 overflows at
// 1.3e154, and u^2, u^3, ... earlier, while the value is still a number).
// b is formed as e^2 a for |e| <= 1 and as 1 - a beyond, so that neither
// form cancels; a is (1/e)/e once e^2 overflows.
inline double cauchit_inv(double e, int k) {
    double a = (std::fabs(e) > 1e150) ? (1.0 / e) / e : 1.0 / (1.0 + e * e);
    double b = (std::fabs(e) <= 1.0) ? e * e * a : 1.0 - a;
    double ea = e * a;
    switch (k) {
    case 1: return a / M_PI;
    case 2: return -2.0 * ea * a / M_PI;
    case 3: return 2.0 * (3.0 * b - a) * a * a / M_PI;
    case 4: return 24.0 * ea * (a - b) * a * a / M_PI;
    case 5: return 24.0 * (a * a - 10.0 * a * b + 5.0 * b * b) * a * a * a / M_PI;
    case 6: return -240.0 * ea * a * a * a * (3.0 * b * b - 10.0 * a * b + 3.0 * a * a) / M_PI;
    default: return NA_REAL;
    }
}

// Written in t = tanh(eta) and u = sech^2(eta) = 1 - t^2, with u computed
// from eta as 4a / (1 + a)^2, a = exp(-2|eta|), and not as 1 - t^2: at
// eta = 20 the subtraction returns 0 where u is 1.7e-17. Each order is u
// times a polynomial in t and u, from dt = u and du = -2tu.
inline double rhobit_inv(double eta, int k) {
    // Beyond |eta| = 300 the value is 4 (-2 t)^(k-1) exp(-2 |eta|), t = +-1, to
    // the working precision (the corrections are O(exp(-600))), formed in
    // one exponential so that the result keeps its precision down to the
    // smallest subnormal; through u = 4a / (1 + a)^2 the subnormal a would
    // be multiplied by the coefficient, up to 128, with its rounding error.
    if (std::fabs(eta) > 300.0) {
        if (k < 1 || k > 6) return NA_REAL;
        double s = eta > 0.0 ? 1.0 : -1.0, c = 4.0;
        for (int j = 1; j < k; ++j) c *= -2.0 * s;
        return std::copysign(std::exp(std::log(std::fabs(c)) - 2.0 * std::fabs(eta)), c);
    }
    double t = std::tanh(eta), t2 = t * t;
    double a = std::exp(-2.0 * std::fabs(eta));
    double u = 4.0 * a / ((1.0 + a) * (1.0 + a));
    switch (k) {
    case 1: return u;
    case 2: return -2.0 * t * u;
    case 3: return u * (4.0 * t2 - 2.0 * u);
    case 4: return 8.0 * t * u * (2.0 * u - t2);
    case 5: return u * (16.0 * u * u - 88.0 * t2 * u + 16.0 * t2 * t2);
    case 6: return -16.0 * t * u * (17.0 * u * u - 26.0 * t2 * u + 2.0 * t2 * t2);
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
