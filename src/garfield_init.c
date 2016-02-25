// This file was automatically generated by 'Kmisc::registerFunctions()'

#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

extern int garfield_prep(char **prune_dir, char **clump_dir, char **maf_dir, char **pv_dir, char **annot_dir, char **excl, char **prep_file);
extern int garfield_perm(char **input, char **link, char **out_file, char **p_thresh, char **pt_thresh, char **npermut_par, char **nannot_par, char **nqmaf_par, char **nqntag_par, char **nqtssd_par, char **optim_mode_par, char **minit_par, char **thresh);


//SEXP garfield_prep(SEXP bampathSEXP, SEXP grSEXP, SEXP mapqualSEXP, SEXP binsizeSEXP, SEXP shiftSEXP, SEXP ssSEXP, SEXP maskSEXP, SEXP pe_midSEXP, SEXP maxfraglengthSEXP, SEXP maxgapSEXP);
//SEXP garfield_perm(SEXP bampathSEXP, SEXP grSEXP, SEXP mapqualSEXP, SEXP maskSEXP, SEXP tspanSEXP, SEXP maxfraglengthSEXP, SEXP maxgapSEXP);
static const R_ExternalMethodDef externalMethods[]  = {
  {"garfield_prep", (DL_FUNC) &garfield_prep, 7},
  {"garfield_perm", (DL_FUNC) &garfield_perm, 13},
  {NULL, NULL, 0}
};
//cMethods
void R_init_garfield(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, NULL, externalMethods);
  R_useDynamicSymbols(info, TRUE);
}
