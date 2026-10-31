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
NumericVector lk_logit_inv_scaled_cpp(NumericVector eta, int k, double s) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = lf7::logit_inv_scaled(eta[i], k, s);
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
        case 6: {
            double t6 = t * t * t * t * t * t, q6 = q * q * q * q * q * q;
            out[i] = -120.0 / t6 + 120.0 / q6;
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
        case 6: {
            double e2 = e * e, p3 = phi * phi * phi;
            out[i] = e * (127.0 + e2 * (326.0 + 120.0 * e2)) / (p3 * p3);
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
        case 6: {
            double L2 = L * L, L3 = L2 * L, v2 = v * v;
            out[i] = -(120.0 + L * (360.0 + L * (510.0 + L * (450.0 + L * (274.0 + 120.0 * L))))) /
                (v2 * v2 * v2 * L3 * L3);
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
    // the denominator is (theta log theta)^k, formed from w = theta log theta
    // and not as theta^k times log(theta)^k: theta^k is subnormal for
    // theta near 1e-155 at k = 2, where the value is still a number
    for (R_xlen_t i = 0; i < n; ++i) {
        double t = theta[i], l = std::log(t), w = t * l;
        double w2 = w * w;
        switch (k) {
        case 1: out[i] = -1.0 / w; break;
        case 2: out[i] = (1.0 + l) / w2; break;
        case 3: out[i] = -(2.0 + 3.0 * l + 2.0 * l * l) / (w2 * w); break;
        case 4: {
            double l2 = l * l;
            out[i] = (6.0 + 12.0 * l + 11.0 * l2 + 6.0 * l * l2) / (w2 * w2);
            break;
        }
        case 5: {
            double l2 = l * l;
            out[i] = -(24.0 + 60.0 * l + 70.0 * l2 + 50.0 * l2 * l +
                       24.0 * l2 * l2) / (w2 * w2 * w);
            break;
        }
        case 6: {
            double w3 = w2 * w;
            out[i] = (120.0 + l * (360.0 + l * (510.0 + l * (450.0 + l * (274.0 + 120.0 * l))))) /
                (w3 * w3);
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
        case 6: {
            double p2 = M_PI * M_PI, e2 = e * e;
            out[i] = 16.0 * p2 * p2 * p2 * u * e * (17.0 + e2 * (60.0 + 45.0 * e2));
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
        // u = (1 - t)(1 + t) and not 1 - t * t: the product t * t is rounded to
        // a spacing of 1e-16 near 1, which is the size of u close to the ends
        double t = theta[i], u = (1.0 - t) * (1.0 + t);
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
        case 6: {
            double u2 = u * u, t2 = t * t;
            out[i] = 240.0 * t * (3.0 + t2 * (10.0 + 3.0 * t2)) / (u2 * u2 * u2);
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
    // order k >= 2 is (-1)^(k-1) exp(-z) N_k(e) q^(k-1) / u, with e = exp(-z),
    // z = a theta, u = 1 - e, q = a / u and N_k the Eulerian numerators
    // (1; 1, 1; 1, 4, 1; 1, 11, 11, 1; 1, 26, 66, 26, 1). The powers of u and
    // of a are never formed alone: a / u is about 1 / theta, so neither the
    // denominator underflows nor the numerator overflows before the value
    // does, and exp(-z) enters through exp_times() so that a value that is
    // still a number past z = 745 is not flushed by the first factor.
    for (R_xlen_t i = 0; i < n; ++i) {
        double z = a * theta[i];
        double e = std::exp(-z), u = -std::expm1(-z);
        double N;
        switch (k) {
        case 1: out[i] = 1.0 / u; continue;
        case 2: N = -1.0; break;
        case 3: N = 1.0 + e; break;
        case 4: N = -(1.0 + e * (4.0 + e)); break;
        case 5: N = 1.0 + e * (11.0 + e * (11.0 + e)); break;
        case 6: N = -(1.0 + e * (26.0 + e * (66.0 + e * (26.0 + e)))); break;
        default: out[i] = NA_REAL; continue;
        }
        double q = a / u, qp = q;
        for (int j = 2; j < k; ++j) qp *= q;
        out[i] = lf7::exp_times(-z, N * qp / u);
    }
    return out;
}
