// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bounding_box
Rcpp::List bounding_box(Rcpp::IntegerVector runs, Rcpp::StringVector values, Rcpp::IntegerVector index1, Rcpp::StringVector ref_chr1, Rcpp::IntegerVector ref_start1, Rcpp::IntegerVector ref_end1, Rcpp::IntegerVector index2, Rcpp::StringVector ref_chr2, Rcpp::IntegerVector ref_start2, Rcpp::IntegerVector ref_end2, Rcpp::LogicalVector reflect);
RcppExport SEXP _GenomicInteractions_bounding_box(SEXP runsSEXP, SEXP valuesSEXP, SEXP index1SEXP, SEXP ref_chr1SEXP, SEXP ref_start1SEXP, SEXP ref_end1SEXP, SEXP index2SEXP, SEXP ref_chr2SEXP, SEXP ref_start2SEXP, SEXP ref_end2SEXP, SEXP reflectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type runs(runsSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type index1(index1SEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type ref_chr1(ref_chr1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_start1(ref_start1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_end1(ref_end1SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type index2(index2SEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type ref_chr2(ref_chr2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_start2(ref_start2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ref_end2(ref_end2SEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type reflect(reflectSEXP);
    rcpp_result_gen = Rcpp::wrap(bounding_box(runs, values, index1, ref_chr1, ref_start1, ref_end1, index2, ref_chr2, ref_start2, ref_end2, reflect));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GenomicInteractions_bounding_box", (DL_FUNC) &_GenomicInteractions_bounding_box, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_GenomicInteractions(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}