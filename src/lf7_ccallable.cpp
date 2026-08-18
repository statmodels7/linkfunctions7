#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <cfloat>
#include <cmath>

// The scalar C entry points of the fast route piano_parallel.txt section 2a
// describes: a package that consumes them (modelterms7's score-driven
// filter) resolves them once with R_GetCCallable and its loop then calls
// plain function pointers, touching no R API. Only the links whose inverse
// is a short closed formula are covered; an unknown name answers -1 and the
// consumer keeps its R callbacks, so coverage is a speed property and never
// a correctness one.
//
// Every formula here MIRRORS the R method it stands for, expression by
// expression, and a twin test compares the two bit for bit: the identity's
// inverse is eta itself, the log's is exp_floored() -- pmax(exp(eta),
// exp_floor) with exp_floor = (24 / double.xmax)^0.25, the derived guard of
// the exponential links -- and the clamp is link_bounds_clamp() read on one
// value.

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

} // extern "C"

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
}
