#include "lfe.h"

#define USE_FC_LEN_T
#ifndef FCONE
# define FCONE
#endif

SEXP MY_scalecols(SEXP mat, SEXP vec) {
  if(!isMatrix(mat)) error("first argument should be a matrix");
  mybigint_t col = ncols(mat), row = nrows(mat);
  if(isMatrix(vec)) {
    if(row != nrows(vec)) error("Rows of matrix should be the same as rows of vector");
    double *cmat = REAL(mat);
    double *cvec = REAL(AS_NUMERIC(vec));
    for(mybigint_t j = 0; j < col; j++) {
      double *cc = &cmat[j*row];
      for(mybigint_t i = 0; i < row; i++) {
	double tmp = 0.0;
	for(int k = 0; k < ncols(vec); k++)
	  tmp += cc[i]*cvec[i+k*row];
	cc[i] = tmp;
      }
    }
  } else {
    if(row != LENGTH(vec)) error("length of vector %d is different from number of rows %d",LENGTH(vec),row);
    double *cmat = REAL(mat);
    double *cvec = REAL(AS_NUMERIC(vec));
    for(mybigint_t j = 0; j < col; j++) {
      double *cc = &cmat[j*row];
      for(mybigint_t i = 0; i < row; i++)
	cc[i] *= cvec[i];
    }
  }
  return mat;
}

/*
  compute X + beta x Y
  where X and Y are matrices of column vectors, and beta a vector
  of length the number of columns of X and Y. Each column of Y should
  be scaled by the corresponding entry of beta
  The same operation can be done in R with
  X + t(beta * t(Y)), so this is a utility to save some costly transposes.
  used in cgsolve()
 */
SEXP MY_pdaxpy(SEXP inX, SEXP inY, SEXP inbeta) {
  mybigint_t col = ncols(inX), row=nrows(inX);
  if(col != ncols(inY) || row != nrows(inY))
    error("X and Y should have the same shape");
  if(LENGTH(inbeta) != col)
    error("beta should have the same length as the number of columns of Y");
  double *X = REAL(inX);
  double *Y = REAL(inY);
  double *beta = REAL(inbeta);
  SEXP res;
  PROTECT(res = allocMatrix(REALSXP,row,col));
  double *pres = REAL(res);
  for(mybigint_t j= 0; j < col; j++) {
    double b = beta[j];
    double *out = &pres[j*row];
    double *xin = &X[j*row];
    double *yin = &Y[j*row];
    for(mybigint_t i=0; i < row; i++) {
      out[i] = xin[i] + b*yin[i];
    }
  }
  UNPROTECT(1);
  return(res);
}

/*
  compute inner product pairwise of the columns of matrices X and Y
  This is diag(crossprod(X,Y)).
 */
SEXP MY_piproduct(SEXP inX, SEXP inY) {
  mybigint_t col = ncols(inX), row=nrows(inX);
  if(col != ncols(inY) || row != nrows(inY))
    error("X and Y should have the same shape");
  double *X = REAL(inX);
  double *Y = REAL(inY);
  SEXP res;
  PROTECT(res = allocVector(REALSXP, col));
  double *pres = REAL(res);
  for(mybigint_t j= 0; j < col; j++) {
    double *xin = &X[j*row];
    double *yin = &Y[j*row];
    pres[j] = 0.0;
    for(mybigint_t i= 0; i < row; i++) {
      pres[j] += xin[i]*yin[i];
    }
  }
  UNPROTECT(1);
  return(res);
}
// copy-free dimnames<-

SEXP MY_setdimnames(SEXP obj, SEXP nm) {
  if(!isNull(obj)) setAttrib(obj, R_DimNamesSymbol, nm);
  return(R_NilValue);
}

/* Compute and return alpha * bread %*% meat %*% bread */

SEXP MY_sandwich(SEXP inalpha, SEXP inbread, SEXP inmeat) {
  double alpha = REAL(AS_NUMERIC(inalpha))[0];
  if(!isMatrix(inbread)) error("bread must be a matrix");
  if(!isMatrix(inmeat)) error("bread must be a matrix");
  if(ncols(inbread) != nrows(inbread)) 
    error("bread must be square matrix");

  if(ncols(inmeat) != nrows(inmeat)) 
    error("meat must be square matrix");

  if(ncols(inmeat) != ncols(inbread))
    error("bread and meat must have the same size");

  int N = ncols(inmeat);
  double *meat = REAL(inmeat);
  double *bread = REAL(inbread);
  SEXP ret;
  double *tmp;
  tmp = (double*) R_alloc(N*N,sizeof(double));
  PROTECT(ret = allocMatrix(REALSXP, N, N));
  double *out = REAL(ret);
  double zero = 0.0;
  double one = 1.0;

  F77_CALL(dsymm)("R","U",&N,&N,&one,bread,&N,meat,&N,&zero,tmp,&N FCONE FCONE);
  F77_CALL(dsymm)("L","U",&N,&N,&alpha,bread,&N,tmp,&N,&zero,out,&N FCONE FCONE);

  UNPROTECT(1);
  return ret;
}


// copy-free dsyrk C = beta*C + alpha * A' A
SEXP MY_dsyrk(SEXP inbeta, SEXP inC, SEXP inalpha, SEXP inA) {
  double beta = REAL(AS_NUMERIC(inbeta))[0];
  double alpha = REAL(AS_NUMERIC(inalpha))[0];
  if(!isMatrix(inC)) error("C must be a matrix");
  if(!isMatrix(inA)) error("A must be a matrix");

  if(ncols(inC) != nrows(inC)) {
    error("C must be a square matrix, it is %d x %d",nrows(inC), ncols(inC));
  }
  int N = nrows(inC);
  double *C = REAL(inC);
  if(ncols(inA) != ncols(inC)) {
    error("A (%d x %d) must have the same number of columns as C (%d x %d)",nrows(inA),ncols(inA),nrows(inC),nrows(inC));
  }
  int K = nrows(inA);
  double *A = REAL(inA);
  F77_CALL(dsyrk)("U","T",&N, &K, &alpha, A, &K, &beta, C, &N FCONE FCONE);
  // fill in the lower triangular part
  for(mybigint_t row=0; row < N; row++) {
    for(mybigint_t col=0; col < row; col++) {
      C[col*N + row] = C[row*N+col];
    }
  }
  return R_NilValue;
}


// debugging memory copy
SEXP MY_address(SEXP x) {
  char chr[30];
  sprintf(chr, "adr=%p, named=%d", (void*)x, NAMED(x));
  return(mkString(chr));
}

SEXP MY_named(SEXP x, SEXP n) {
  if(isNull(n)) {
    // return NAMED status
    SEXP res = allocVector(INTSXP, 1);
    PROTECT(res);
    INTEGER(res)[0] = NAMED(x);
    setAttrib(res,install("x"),x);
    UNPROTECT(1);
    return(res);
  }
  // set named status. Use this at your own peril. Seriously. Things may ...
  SET_NAMED(x,INTEGER(n)[0]);
  return(x);
}

/* Trickery to check for interrupts when using threads */
static void chkIntFn(void *dummy) {
  if(dummy==NULL){}; //avoid pedantic warning about unused parameter
  R_CheckUserInterrupt();
}

/* this will call the above in a top-level context so it won't 
   longjmp-out of context */

int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == 0);
}


/* More trickery, we can't printf in threads since the R API is not
   thread-safe.  So set up a message stack.  This is pushed in the threads
   and popped in the main thread. */

#define MSGLIM 256
static char *msgstack[MSGLIM];
static int msgptr;

void initmsg() {
  msgptr = 0;
}
/* Craft our own strdup, it's not supported everywhere */
static char *mystrdup(char *s) {
  char *sc = (char*)malloc(strlen(s)+1);
  if(sc != NULL) strcpy(sc,s);
  return(sc);
}
void pushmsg(char *s, LOCK_T lock) {
#ifdef NOTHREADS  
  REprintf(s);
#else
  LOCK(lock);
  if(msgptr < MSGLIM) {
    msgstack[msgptr++] = mystrdup(s);
  }
  UNLOCK(lock);
#endif
}

void printmsg(LOCK_T lock) {
#ifdef NOTHREADS
  return;
#else
  char *s;
  int i;
  LOCK(lock);
  for(i = 0; i < msgptr; i++) {
    s = msgstack[i];
    if(s != NULL) {
      REprintf(s);
      free(s);
    }
  }
  msgptr = 0;
  UNLOCK(lock);
#endif
}
