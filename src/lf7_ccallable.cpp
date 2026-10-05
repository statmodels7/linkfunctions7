#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <cfloat>
#include <cmath>
#include "link_points.h"

// The scalar C entry points of the fast route piano_parallel.txt section 2a
// describes: a package that consumes them (modelterms7's score-driven
// filter) resolves them once with R_GetCCallable and its loop then calls
// plain function pointers, touching no R API. An unknown link answers -1
// and the consumer keeps its R callbacks, so coverage is a speed property
// and never a correctness one.
//
// Two generations of entry. lf7_scalar_id()/lf7_inv12() take the link's
// NAME and cover the identity and the log; lf7_class_id()/lf7_inv12p()
// take the S7 CLASS name and the link's own parameters (link_scalar_route()
// returns both) and cover every link class the package ships.
//
// Every formula here MIRRORS the R method it stands for, expression by
// expression: the identity's inverse is eta itself, the log's is
// exp_floored() -- pmax(exp(eta), exp_floor) with exp_floor =
// (24 / double.xmax)^0.25, the derived guard of the exponential links --
// the transcendental links' derivatives are the point functions of
// link_points.h that the vector kernels loop over, the distribution
// functions are R's own (plogis, pnorm, pcauchy), and the clamp is
// link_bounds_clamp() read on one value.

extern "C" {

int lf7_scalar_id(const char* name) {
    if (std::strcmp(name, "identity") == 0) return 0;
    if (std::strcmp(name, "log") == 0) return 1;
    return -1;
}

// h = linkinv(eta), h1 = dlinkinv, h2 = d2linkinv, the three the chain rule
// of one component wants
void lf7_inv12(int id, double eta, double* h, double* h1, double* h2) {
    switch (id) {
    case 0:
        *h = eta; *h1 = 1.0; *h2 = 0.0;
        break;
    case 1: {
        // exp_floored(): the floor is the bound the fourth forward
        // derivative of the log link derives, written exactly as the R
        // constant computes it
        static const double floor_ = std::pow(24.0 / DBL_MAX, 0.25);
        double v = std::exp(eta);
        if (!ISNAN(v) && v < floor_) v = floor_;
        *h = v; *h1 = v; *h2 = v;
        break;
    }
    default:
        *h = R_NaN; *h1 = R_NaN; *h2 = R_NaN;
    }
}

// link_bounds_clamp() on one value: infinities to the largest finite double
// of their sign first, then EXACT equality with a finite bound moved to the
// nearest representable value strictly inside
double lf7_clamp(double th, double lwr, double upr) {
    if (th == R_PosInf) th = DBL_MAX;
    else if (th == R_NegInf) th = -DBL_MAX;
    if (R_FINITE(lwr) && !ISNAN(th) && th == lwr) {
        th = (lwr == 0.0) ? DBL_MIN : lwr + std::fabs(lwr) * DBL_EPSILON;
    }
    if (R_FINITE(upr) && !ISNAN(th) && th == upr) {
        th = (upr == 0.0) ? -DBL_MIN : upr - std::fabs(upr) * DBL_EPSILON;
    }
    return th;
}

// the id of a link class is its position here
static const char* const lf7_classes[] = {
    "IdentityLink",       // 0
    "LogLink",            // 1
    "LogitLink",          // 2
    "ProbitLink",         // 3
    "ClogLogLink",        // 4
    "LogLogLink",         // 5
    "CauchitLink",        // 6
    "RhobitLink",         // 7
    "SqrtLink",           // 8
    "InverseLink",        // 9
    "InverseSqLink",      // 10
    "PowerLink",          // 11, par = lambda
    "SoftplusLink",       // 12, par = a
    "DoublyBoundedLink",  // 13, par = (lwr, width)
    "LowerBoundedLink",   // 14, par = lwr
    "UpperBoundedLink"    // 15, par = upr
};

int lf7_class_id(const char* cls) {
    const int n = sizeof(lf7_classes) / sizeof(lf7_classes[0]);
    for (int i = 0; i < n; ++i) {
        if (std::strcmp(cls, lf7_classes[i]) == 0) return i;
    }
    return -1;
}

// h = linkinv(eta) before the generic's clamp, h1 = dlinkinv, h2 =
// d2linkinv, each as the class's R method computes it; par holds the
// link's own parameters in the order link_scalar_route() returns them
void lf7_inv12p(int id, const double* par, double eta, double* h,
                double* h1, double* h2) {
    // const_like() and na_from(): NA where eta is NA
    const bool na = ISNAN(eta);
    switch (id) {
    case 0:
        *h = eta; *h1 = na ? NA_REAL : 1.0; *h2 = na ? NA_REAL : 0.0;
        break;
    case 1:
        *h = lf7::exp_floored(eta); *h1 = *h; *h2 = *h;
        break;
    case 2:
        *h = R::plogis(eta, 0.0, 1.0, 1, 0);
        *h1 = lf7::logit_inv(eta, 1); *h2 = lf7::logit_inv(eta, 2);
        break;
    case 3:
        *h = R::pnorm(eta, 0.0, 1.0, 1, 0);
        *h1 = lf7::probit_inv(eta, 1); *h2 = lf7::probit_inv(eta, 2);
        break;
    case 4: {
        // pmax(pmin(-expm1(-exp(eta)), 1 - eps), exp_floor)
        static const double floor_ = std::pow(24.0 / DBL_MAX, 0.25);
        double v = -std::expm1(-std::exp(eta));
        if (!ISNAN(v)) {
            if (v > 1.0 - DBL_EPSILON) v = 1.0 - DBL_EPSILON;
            if (v < floor_) v = floor_;
        }
        *h = v;
        *h1 = lf7::cloglog_inv(eta, 1); *h2 = lf7::cloglog_inv(eta, 2);
        break;
    }
    case 5:
        *h = std::exp(-std::exp(-eta));
        *h1 = lf7::loglog_inv(eta, 1); *h2 = lf7::loglog_inv(eta, 2);
        break;
    case 6:
        *h = R::pcauchy(eta, 0.0, 1.0, 1, 0);
        *h1 = lf7::cauchit_inv(eta, 1); *h2 = lf7::cauchit_inv(eta, 2);
        break;
    case 7:
        *h = std::tanh(eta);
        *h1 = lf7::rhobit_inv(eta, 1); *h2 = lf7::rhobit_inv(eta, 2);
        break;
    case 8:
        *h = lf7::rpow(eta, 2.0); *h1 = 2.0 * eta;
        *h2 = na ? NA_REAL : 2.0;
        break;
    case 9:
        *h = 1.0 / eta; *h1 = -1.0 / lf7::rpow(eta, 2.0);
        *h2 = 2.0 / lf7::rpow(eta, 3.0);
        break;
    case 10:
        *h = 1.0 / std::sqrt(eta);
        *h1 = -1.0 / (2.0 * lf7::rpow(eta, 1.5));
        *h2 = 3.0 / (4.0 * lf7::rpow(eta, 2.5));
        break;
    case 11: {
        double k = 1.0 / par[0];
        *h = lf7::rpow(eta, 1.0 / par[0]);
        *h1 = na ? NA_REAL : k * lf7::rpow(eta, k - 1.0);
        *h2 = na ? NA_REAL : k * (k - 1.0) * lf7::rpow(eta, k - 2.0);
        break;
    }
    case 12: {
        double a = par[0];
        // pmax(0, eta) + log1p(exp(-abs(a * eta))) / a
        double m = na ? NA_REAL : (eta > 0.0 ? eta : 0.0);
        *h = m + std::log1p(std::exp(-std::fabs(a * eta))) / a;
        *h1 = R::plogis(a * eta, 0.0, 1.0, 1, 0);
        *h2 = a * lf7::logit_inv(a * eta, 1);
        break;
    }
    case 13:
        *h = par[0] + par[1] * R::plogis(eta, 0.0, 1.0, 1, 0);
        *h1 = par[1] * lf7::logit_inv(eta, 1);
        *h2 = par[1] * lf7::logit_inv(eta, 2);
        break;
    case 14: {
        double e = lf7::exp_floored(eta);
        *h = par[0] + e; *h1 = e; *h2 = e;
        break;
    }
    case 15: {
        double e = lf7::exp_floored(eta);
        *h = par[0] - e; *h1 = -e; *h2 = -e;
        break;
    }
    default:
        *h = R_NaN; *h1 = R_NaN; *h2 = R_NaN;
    }
}

} // extern "C"

// the class entries, for the same tests
// [[Rcpp::export]]
Rcpp::List lf7_class_probe(std::string cls, Rcpp::NumericVector par,
                           Rcpp::NumericVector eta) {
    int id = lf7_class_id(cls.c_str());
    int n = eta.size();
    Rcpp::NumericVector h(n), h1(n), h2(n);
    for (int i = 0; i < n; ++i) {
        double a, b, c;
        lf7_inv12p(id, par.begin(), eta[i], &a, &b, &c);
        h[i] = a; h1[i] = b; h2[i] = c;
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["h"] = h,
                              Rcpp::_["h1"] = h1, Rcpp::_["h2"] = h2);
}

// exposed to this package's own tests, so the twin comparison against the R
// methods lives where the formulas do
// [[Rcpp::export]]
Rcpp::List lf7_scalar_probe(std::string name, Rcpp::NumericVector eta,
                            Rcpp::NumericVector bounds) {
    int id = lf7_scalar_id(name.c_str());
    int n = eta.size();
    Rcpp::NumericVector h(n), h1(n), h2(n), hc(n);
    for (int i = 0; i < n; ++i) {
        double a, b, c;
        lf7_inv12(id, eta[i], &a, &b, &c);
        h[i] = a; h1[i] = b; h2[i] = c;
        hc[i] = lf7_clamp(a, bounds[0], bounds[1]);
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["h"] = h,
                              Rcpp::_["h1"] = h1, Rcpp::_["h2"] = h2,
                              Rcpp::_["clamped"] = hc);
}

// [[Rcpp::init]]
void lf7_register_ccallable(DllInfo* dll) {
    R_RegisterCCallable("linkfunctions7", "lf7_scalar_id",
                        (DL_FUNC) lf7_scalar_id);
    R_RegisterCCallable("linkfunctions7", "lf7_inv12", (DL_FUNC) lf7_inv12);
    R_RegisterCCallable("linkfunctions7", "lf7_clamp", (DL_FUNC) lf7_clamp);
    R_RegisterCCallable("linkfunctions7", "lf7_class_id",
                        (DL_FUNC) lf7_class_id);
    R_RegisterCCallable("linkfunctions7", "lf7_inv12p", (DL_FUNC) lf7_inv12p);
}
