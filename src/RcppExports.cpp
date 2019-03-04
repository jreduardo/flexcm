// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_logz
NumericVector compute_logz(NumericVector loglambda, double nu);
RcppExport SEXP _flexcm_compute_logz(SEXP loglambdaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type loglambda(loglambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_logz(loglambda, nu));
    return rcpp_result_gen;
END_RCPP
}
// compute_logk
NumericVector compute_logk(NumericVector mu, NumericVector lmu, double phi, double lphi);
RcppExport SEXP _flexcm_compute_logk(SEXP muSEXP, SEXP lmuSEXP, SEXP phiSEXP, SEXP lphiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lmu(lmuSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type lphi(lphiSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_logk(mu, lmu, phi, lphi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flexcm_compute_logz", (DL_FUNC) &_flexcm_compute_logz, 2},
    {"_flexcm_compute_logk", (DL_FUNC) &_flexcm_compute_logk, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_flexcm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
