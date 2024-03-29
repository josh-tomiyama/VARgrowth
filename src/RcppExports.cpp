// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4BoundedLinearTrendModel_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4DeterministicLinearTrendModel_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4LinearTrendModel_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4LinearTrendModel2_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4LinearTrendModel_noncentered_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4example_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4BoundedLinearTrendModel_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4BoundedLinearTrendModel_mod, 0},
    {"_rcpp_module_boot_stan_fit4DeterministicLinearTrendModel_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4DeterministicLinearTrendModel_mod, 0},
    {"_rcpp_module_boot_stan_fit4LinearTrendModel_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4LinearTrendModel_mod, 0},
    {"_rcpp_module_boot_stan_fit4LinearTrendModel2_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4LinearTrendModel2_mod, 0},
    {"_rcpp_module_boot_stan_fit4LinearTrendModel_noncentered_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4LinearTrendModel_noncentered_mod, 0},
    {"_rcpp_module_boot_stan_fit4example_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4example_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_VARgrowth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
