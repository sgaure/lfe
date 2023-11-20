#include "config.h"

/* different syntax in pthread_setname_np between platforms, disable for now */
// #undef HAVE_THREADNAME
#if defined(_WIN32) || defined(_WIN64) || defined(WIN64)
#define WIN
#undef HAVE_THREADNAME
#endif

#ifdef WIN
#include <windows.h>
#else

#ifndef NOTHREADS
#ifdef HAVE_THREADNAME
#define _GNU_SOURCE /* to find pthread_setname_np */
#endif

#include <semaphore.h>
#include <pthread.h>
#endif

#endif

#ifdef HAVE_THREADNAME
#ifdef __linux__
// Linux
int pthread_setname_np(pthread_t thread, const char *name);
#define STNAME(s) pthread_setname_np(pthread_self(), s)
#elif __APPLE__
// Mac OS X: must be set from within the thread (can't specify thread ID)
int pthread_setname_np(const char *);
#define STNAME(s) pthread_setname_np(s)
#elif __FreeBSD__
// FreeBSD & OpenBSD: function name is slightly different, and has no return value
void pthread_set_name_np(pthread_t tid, const char *name);
#define STNAME(s) pthread_set_name_np(pthread_self(), s)
#else
// NetBSD: name + arg work like printf(name, arg)
// int pthread_setname_np(pthread_t thread, const char *name, void *arg);
// #define STNAME(s) pthread_setname_np(pthread_self(), s)
// no thread name
#define STNAME(s)
#endif
#else
// Don't have thread name
#define STNAME(s)
#endif

#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rversion.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <R_ext/BLAS.h>

#if defined(R_VERSION) && R_VERSION >= R_Version(3, 0, 0)
typedef R_xlen_t mybigint_t;
#else
typedef int mybigint_t;
#endif
/* If the the number of G() terms times the number of observations
   exceeds the 2^31 4-byte integer limit, define mysize_t as size_t
   This will increase the memory usage, so we wait until it's needed.
*/

#ifdef HUGE_INT
typedef R_xlen_t mysize_t;
#else
typedef int mysize_t;
#endif

/* Locking macros */
#ifdef NOTHREADS
#define LOCK_T int *
#define LOCK(l)
#define UNLOCK(l)
#else

#ifdef WIN
#define LOCK_T HANDLE
#define LOCK(l) WaitForSingleObject(l, INFINITE)
#define UNLOCK(l) ReleaseMutex(l)
#else
#define LOCK_T pthread_mutex_t *
#define LOCK(l) (void)pthread_mutex_lock(l)
#define UNLOCK(l) (void)pthread_mutex_unlock(l)
#endif
#endif

/* My internal definition of a factor */
typedef struct
{
  /* group[i] is the level of observation i */
  int *group;
  /* invgpsize[j] is the 1/(size of level j) */
  double *invgpsize;
  double *gpsize;
  int *gpl;  /* group list */
  int *ii;   /* indices into gpl */
  double *x; /* optional interaction covariate */
  int nlevels;
  int oneiter;
} FACTOR;

/* Routines used in more than one source file */
FACTOR **makefactors(SEXP flist, int allowmissing, double *weights);
extern int checkInterrupt(void);
extern void initmsg(void);
void pushmsg(char *s, LOCK_T lock);
void printmsg(LOCK_T lock);

/* R interface routines */
SEXP MY_kaczmarz(SEXP flist, SEXP vlist, SEXP Reps, SEXP initial, SEXP Rcores);
SEXP MY_wwcomp(SEXP flist);
SEXP MY_conncomp(SEXP flist);
SEXP MY_demeanlist(SEXP vlist, SEXP flist, SEXP Ricpt, SEXP Reps,
                   SEXP scores, SEXP quiet, SEXP gkacc, SEXP Rmeans,
                   SEXP weights, SEXP Rscale, SEXP attrs);
SEXP MY_scalecols(SEXP mat, SEXP vec);
SEXP MY_pdaxpy(SEXP inX, SEXP inY, SEXP inbeta);
SEXP MY_piproduct(SEXP inX, SEXP inY);
SEXP MY_setdimnames(SEXP obj, SEXP nm);
SEXP MY_dsyrk(SEXP inbeta, SEXP inC, SEXP inalpha, SEXP inA);
SEXP MY_sandwich(SEXP inalpha, SEXP inbread, SEXP inmeat);
SEXP MY_address(SEXP x);
// SEXP MY_named(SEXP x, SEXP n);
SEXP inplace(SEXP x);
SEXP crowsum(SEXP Rmat, SEXP Rfactor, SEXP Rmean);

// SEXP MY_ppf(SEXP flist, SEXP Rtype);
