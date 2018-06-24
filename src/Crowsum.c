#include "lfe.h"
//#ifdef _OPENMP
//#include <omp.h>
//#endif
SEXP Crowsum(SEXP Rmat, SEXP Rfactor, SEXP Rmean) {
  if(!IS_NUMERIC(Rmat)) error("Only numeric matrices accepted");
  if(!isInteger(Rfactor) && !isFactor(Rfactor)) error("Only factor or integer vector accepted");
  mybigint_t len = length(Rmat);
  mybigint_t cols = 0, rows=0;
  int *f = INTEGER(Rfactor);
  double *mat;
  int nlev;
  SEXP res;
  double *mres;
  int mean = INTEGER(AS_LOGICAL(Rmean))[0];
  int *table = NULL;
  //  int nthr = INTEGER(threads)[0];
  mat = REAL(Rmat);
  if(isMatrix(Rmat)) {
    cols = ncols(Rmat);
    rows = nrows(Rmat);
  } else {
    cols = 1;
    rows = len;
  }
  if(length(Rfactor) != rows) error("matrix/vector must have same length as factor");

  nlev = nlevels(Rfactor);

  for(int i = 0; i < rows; i++) {
    if(f[i] < 1 || ISNA(f[i])) error("Missing levels not supported");
    if(f[i] > nlev) error("Level for %d is %d, too large %d",i,f[i],nlev);
  }

  if(mean) {
    table = (int*) R_alloc(nlev,sizeof(int));
    for(int i = 0; i < nlev; i++) table[i] = 0;
    for(int i = 0; i < rows; i++) table[f[i]-1]++;
  }

  // Allocate resultant matrix
  PROTECT(res = allocMatrix(REALSXP, nlev, cols));

  SEXP dn;
  SEXP rdn = GET_DIMNAMES(Rmat);
  PROTECT(dn = allocVector(VECSXP,2));
  SET_VECTOR_ELT(dn,0,GET_LEVELS(Rfactor));
  if(!isNull(rdn)) SET_VECTOR_ELT(dn,1,VECTOR_ELT(rdn,1));
  SET_DIMNAMES(res,dn);
  UNPROTECT(1);

  mres = REAL(res);
  // Now, run through column by column, summing the levels
  //#pragma omp parallel for num_threads(nthr)
  memset(mres,0,nlev*cols*sizeof(double));
  mres--; // factor levels are 1-based
  for(int k = 0; k < cols; k++,mres+=nlev) {
    for(int i = 0; i < rows; i++) {
      mres[f[i]] += *mat++;
    }
  }
  if(mean) {
    mres = REAL(res);
    for(int k = 0; k < cols; k++, mres+=nlev)
      for(int i = 0; i < nlev; i++) mres[i] /= table[i];
  }
  UNPROTECT(1);
  return(res);
}
