#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP DLMtool_bhnoneq_LL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP DLMtool_doprojPI_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP DLMtool_movfit_Rcpp(SEXP, SEXP, SEXP);
extern SEXP DLMtool_optQ_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP DLMtool_projOpt_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"DLMtool_bhnoneq_LL",   (DL_FUNC) &DLMtool_bhnoneq_LL,    8},
    {"DLMtool_doprojPI_cpp", (DL_FUNC) &DLMtool_doprojPI_cpp, 22},
    {"DLMtool_movfit_Rcpp",  (DL_FUNC) &DLMtool_movfit_Rcpp,   3},
    {"DLMtool_optQ_cpp",     (DL_FUNC) &DLMtool_optQ_cpp,     17},
    {"DLMtool_projOpt_cpp",  (DL_FUNC) &DLMtool_projOpt_cpp,  16},
    {NULL, NULL, 0}
};

void R_init_DLMtool(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
