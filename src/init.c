#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DLMtool_bhnoneq_LL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_doprojPI_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_genLenComp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_LSRA_MCMC_sim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_LSRA_opt_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_movfit_Rcpp(SEXP, SEXP, SEXP);
extern SEXP _DLMtool_optQ_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_popdynCPP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_popdynOneTScpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DLMtool_projOpt_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_DLMtool_bhnoneq_LL",     (DL_FUNC) &_DLMtool_bhnoneq_LL,      8},
  {"_DLMtool_doprojPI_cpp",   (DL_FUNC) &_DLMtool_doprojPI_cpp,   23},
  {"_DLMtool_genLenComp",     (DL_FUNC) &_DLMtool_genLenComp,      9},
  {"_DLMtool_LSRA_MCMC_sim",  (DL_FUNC) &_DLMtool_LSRA_MCMC_sim,  21},
  {"_DLMtool_LSRA_opt_cpp",   (DL_FUNC) &_DLMtool_LSRA_opt_cpp,   10},
  {"_DLMtool_movfit_Rcpp",    (DL_FUNC) &_DLMtool_movfit_Rcpp,     3},
  {"_DLMtool_optQ_cpp",       (DL_FUNC) &_DLMtool_optQ_cpp,       17},
  {"_DLMtool_popdynCPP",      (DL_FUNC) &_DLMtool_popdynCPP,      23},
  {"_DLMtool_popdynOneTScpp", (DL_FUNC) &_DLMtool_popdynOneTScpp, 13},
  {"_DLMtool_projOpt_cpp",    (DL_FUNC) &_DLMtool_projOpt_cpp,    17},
  {NULL, NULL, 0}
};

void R_init_DLMtool(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
