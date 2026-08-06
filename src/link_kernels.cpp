#include <Rcpp.h>
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

// the logistic derivative polynomials in p: shared by the logit (as is),
// the doubly bounded link (scaled by the width) and the softplus (an
// antiderivative of the logistic, so shifted one place down)
static inline double logistic_poly(double p, int k) {
    double pq = p * (1.0 - p);
    switch (k) {
    case 1: return pq;
    case 2: return pq * (1.0 - 2.0 * p);
    case 3: return pq * (1.0 + p * (-6.0 + 6.0 * p));
    default: return pq * (1.0 + p * (-14.0 + p * (36.0 - 24.0 * p)));
    }
}

// [[Rcpp::export]]
NumericVector lk_logit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double p = 1.0 / (1.0 + std::exp(-eta[i]));
        out[i] = logistic_poly(p, k);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_logistic_poly_cpp(NumericVector p, int k) {
    R_xlen_t n = p.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) out[i] = logistic_poly(p[i], k);
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
        default: {
            double t4 = t * t * t * t, q4 = q * q * q * q;
            out[i] = -6.0 / t4 + 6.0 / q4;
        }
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
        default: out[i] = (7.0 * e + 6.0 * e * e * e) / (phi * phi * phi * phi);
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_probit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double e = eta[i], phi = R::dnorm4(e, 0.0, 1.0, 0);
        switch (k) {
        case 1: out[i] = phi; break;
        case 2: out[i] = -e * phi; break;
        case 3: out[i] = (e * e - 1.0) * phi; break;
        default: out[i] = (3.0 * e - e * e * e) * phi;
        }
    }
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
        default: {
            double L2 = L * L;
            out[i] = -(6.0 * L * L2 + 11.0 * L2 + 12.0 * L + 6.0) /
                (v * v * v * v * L2 * L2);
        }
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cloglog_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double w = std::exp(eta[i]);
        double E = std::exp(eta[i] - w);
        switch (k) {
        case 1: out[i] = E; break;
        case 2: out[i] = E * (1.0 - w); break;
        case 3: out[i] = E * (1.0 + w * (-3.0 + w)); break;
        default: out[i] = E * (1.0 + w * (-7.0 + w * (6.0 - w)));
        }
    }
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
        default: {
            double l2 = l * l;
            out[i] = (6.0 + 12.0 * l + 11.0 * l2 + 6.0 * l * l2) /
                (t * t * t * t * l2 * l2);
        }
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_loglog_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double z = std::exp(-eta[i]);
        double E = std::exp(-z);
        switch (k) {
        case 1: out[i] = E * z; break;
        case 2: out[i] = E * z * (z - 1.0); break;
        case 3: out[i] = E * z * (1.0 + z * (-3.0 + z)); break;
        default: out[i] = E * z * (-1.0 + z * (7.0 + z * (-6.0 + z)));
        }
    }
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
        default: out[i] = 8.0 * M_PI * M_PI * M_PI * M_PI * e * u * (2.0 + 3.0 * e * e);
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_cauchit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double e = eta[i], u = 1.0 + e * e;
        switch (k) {
        case 1: out[i] = 1.0 / (M_PI * u); break;
        case 2: out[i] = -2.0 * e / (M_PI * u * u); break;
        case 3: out[i] = 2.0 * (3.0 * e * e - 1.0) / (M_PI * u * u * u); break;
        default: out[i] = 24.0 * e * (1.0 - e * e) / (M_PI * u * u * u * u);
        }
    }
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
        default: out[i] = 24.0 * t * (1.0 + t * t) / (u * u * u * u);
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector lk_rhobit_inv_cpp(NumericVector eta, int k) {
    R_xlen_t n = eta.size();
    NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double t = std::tanh(eta[i]), t2 = t * t;
        switch (k) {
        case 1: out[i] = 1.0 - t2; break;
        case 2: out[i] = -2.0 * t * (1.0 - t2); break;
        case 3: out[i] = -2.0 + 8.0 * t2 - 6.0 * t2 * t2; break;
        default: out[i] = t * (16.0 + t2 * (-40.0 + 24.0 * t2));
        }
    }
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
        default: out[i] = -a * a * a * e * (1.0 + e * (4.0 + e)) / (u * u * u * u);
        }
    }
    return out;
}
