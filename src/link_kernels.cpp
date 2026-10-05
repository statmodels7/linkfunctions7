#include <Rcpp.h>
#include "link_points.h"
using namespace Rcpp;

// Compiled kernels for the transcendental links' derivatives. The formulas
// are the ones the R methods carried, one scalar loop per element: the
// transcendental call is the same in both languages, but an order-four
// polynomial evaluated as vectorized R allocates a dozen temporaries, and
// removing them measured 1.75x at n = 1e4 and 2.0x at 1e6 on the logit's
// fourth inverse derivative. Every kernel takes the order as an argument, so
// each S7 method body stays one line. The independent reference is
// check_link(), which validates every order against numDeriv and the
// inverse function theorem, and the extremes suite, which walks the tails.

// [[Rcpp::export]]
NumericVector lk_logit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::logit_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_logistic_poly_cpp(NumericVector p, int k) {
    R_xlen_t n = p.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::logistic_poly(p[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_logit_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double t = theta[i], q = 1.0 - t;
        switch (k) {
        case 1: out[i] = 1.0 / (t * q); break;
        case 2: out[i] = (2.0 * t - 1.0) / ((t * q) * (t * q)); break;
        case 3: out[i] = 2.0 / (t * t * t) + 2.0 / (q * q * q); break;
        case 4: {
            double t4 = t * t * t * t, q4 = q * q * q * q;
            out[i] = -6.0 / t4 + 6.0 / q4;
            break;
        }
        case 5: {
            double t5 = t * t * t * t * t, q5 = q * q * q * q * q;
            out[i] = 24.0 / t5 + 24.0 / q5;
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_probit_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double e = R::qnorm5(theta[i], 0.0, 1.0, 1, 0);
        double phi = R::dnorm4(e, 0.0, 1.0, 0);
        switch (k) {
        case 1: out[i] = 1.0 / phi; break;
        case 2: out[i] = e / (phi * phi); break;
        case 3: out[i] = (1.0 + 2.0 * e * e) / (phi * phi * phi); break;
        case 4: out[i] = (7.0 * e + 6.0 * e * e * e) / (phi * phi * phi * phi);
            break;
        case 5: {
            double e2 = e * e, p2 = phi * phi;
            out[i] = (7.0 + 46.0 * e2 + 24.0 * e2 * e2) / (p2 * p2 * phi);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_probit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::probit_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cloglog_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double v = 1.0 - theta[i], L = std::log1p(-theta[i]);
        switch (k) {
        case 1: out[i] = -1.0 / (v * L); break;
        case 2: out[i] = -(L + 1.0) / (v * v * L * L); break;
        case 3: out[i] = -(2.0 * L * L + 3.0 * L + 2.0) / (v * v * v * L * L * L); break;
        case 4: {
            double L2 = L * L;
            out[i] = -(6.0 * L * L2 + 11.0 * L2 + 12.0 * L + 6.0) /
                (v * v * v * v * L2 * L2);
            break;
        }
        case 5: {
            double L2 = L * L, v2 = v * v;
            out[i] = -(24.0 * L2 * L2 + 50.0 * L2 * L + 70.0 * L2 +
                       60.0 * L + 24.0) / (v2 * v2 * v * L2 * L2 * L);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cloglog_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::cloglog_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_loglog_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double t = theta[i], l = std::log(t);
        switch (k) {
        case 1: out[i] = -1.0 / (t * l); break;
        case 2: out[i] = (1.0 + l) / (t * t * l * l); break;
        case 3: out[i] = -(2.0 + 3.0 * l + 2.0 * l * l) / (t * t * t * l * l * l); break;
        case 4: {
            double l2 = l * l;
            out[i] = (6.0 + 12.0 * l + 11.0 * l2 + 6.0 * l * l2) /
                (t * t * t * t * l2 * l2);
            break;
        }
        case 5: {
            double l2 = l * l, t2 = t * t;
            out[i] = -(24.0 + 60.0 * l + 70.0 * l2 + 50.0 * l2 * l +
                       24.0 * l2 * l2) / (t2 * t2 * t * l2 * l2 * l);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_loglog_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::loglog_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cauchit_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double e = R::qcauchy(theta[i], 0.0, 1.0, 1, 0);
        double u = 1.0 + e * e;
        switch (k) {
        case 1: out[i] = M_PI * u; break;
        case 2: out[i] = 2.0 * M_PI * M_PI * e * u; break;
        case 3: out[i] = 2.0 * M_PI * M_PI * M_PI * u * (1.0 + 3.0 * e * e); break;
        case 4: out[i] = 8.0 * M_PI * M_PI * M_PI * M_PI * e * u * (2.0 + 3.0 * e * e);
            break;
        case 5: {
            double p2 = M_PI * M_PI, e2 = e * e;
            out[i] = 8.0 * p2 * p2 * M_PI * u * (2.0 + 15.0 * e2 + 15.0 * e2 * e2);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cauchit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::cauchit_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_rhobit_fwd_cpp(NumericVector theta, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double t = theta[i], u = 1.0 - t * t;
        switch (k) {
        case 1: out[i] = 1.0 / u; break;
        case 2: out[i] = 2.0 * t / (u * u); break;
        case 3: out[i] = (2.0 + 6.0 * t * t) / (u * u * u); break;
        case 4: out[i] = 24.0 * t * (1.0 + t * t) / (u * u * u * u); break;
        case 5: {
            double u2 = u * u, t2 = t * t;
            out[i] = (24.0 + 240.0 * t2 + 120.0 * t2 * t2) / (u2 * u2 * u);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_rhobit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::rhobit_inv(eta[i], k);
    return out;
}

// [[Rcpp::export]]
NumericVector lk_softplus_fwd_cpp(NumericVector theta, double a, int k) {
    R_xlen_t n = theta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double z = a * theta[i];
        double e = std::exp(-z), u = -std::expm1(-z);
        switch (k) {
        case 1: out[i] = 1.0 / u; break;
        case 2: out[i] = -a * e / (u * u); break;
        case 3: out[i] = a * a * e * (1.0 + e) / (u * u * u); break;
        case 4: out[i] = -a * a * a * e * (1.0 + e * (4.0 + e)) / (u * u * u * u);
            break;
        case 5: {
            double a2 = a * a, u2 = u * u;
            out[i] = a2 * a2 * e * (1.0 + e * (11.0 + e * (11.0 + e))) / (u2 * u2 * u);
            break;
        }
        default: out[i] = NA_REAL;
        }
    }
    return out;
}
