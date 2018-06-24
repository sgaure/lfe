/*
  $Id: kaczmarz.c 1915 2016-04-05 12:56:27Z sgaure $
*/
#include "lfe.h"
/* Need sprintf */
#include <stdio.h>  

typedef struct {
  int nowdoing;
  
  double **source;
  double **target;
  FACTOR **factors;
  int e;
  int N;
  int threadnum;
  int stop;
  int numvec;
  double eps;
  double *tol;
  //  double **newR;
  //  int **indices;
  mysize_t **work;
  LOCK_T lock;
#ifndef NOTHREADS
#ifdef WIN
#else
  int running;
#ifdef HAVE_SEM
  sem_t finished;
#endif
#endif
#endif
} KARG;

static double kaczmarz(FACTOR *factors[],int e, mysize_t N, double *R, double *x,
		       double eps, mysize_t *work, int *stop, LOCK_T lock) {

  /* The factors define a matrix D, we will solve Dx = R
     There are N rows in D, each row a_i contains e non-zero entries, one
     for each factor, the level at that position.
     We iterate on k, start with x=0, an iteration consists
     of adding a multiple of row i (with i = k %% N), the multiplying coefficient
     is (R[i] - (a_i,x))/e 

     To get better memory locality, we create an (e X N)-matrix with the non-zero indices

     Now, as an afterthought we have made possible to interact a covariate
     with each factor, it's in factors[i]->x
     The theory is exactly the same, but the projection is different.
     I.e. a_i are not 0 or one, and we should not divide by e, but by
     ||a_i||^2.

  */

  double norm2;
  double prevdiff,neweps;
  double c,sd;
  int iter=0;

  int ie;
  int newN;

  /* We should remove duplicates, at least when they're consecutive.
     Non-consecutives duplicates are just alternating projections and may
     not be too inefficient.  Besides, removing them is more work...
     If any factor is interacted with a covariate, we don't remove
     anything
  */

  int hasinteract = 0;
  for(int i = 0; i < e; i++) if(NULL != factors[i]->x) hasinteract = 1;
  mysize_t workpos = 0;

  /* Do the doubles first to keep alignment */
  double *newR = (double *) work;
  mysize_t *indices = (mysize_t *) &work[workpos += N*sizeof(double)/sizeof(mysize_t)];
  mysize_t *perm = (mysize_t *) &work[workpos += e*N];
  mysize_t *prev = (mysize_t *) &work[workpos += N];
  mysize_t *this = (mysize_t *) &work[workpos += e];

  if(!hasinteract) for(int i = 0; i < e; i++)  prev[i] = 0;
  newN = 0;
  ie = 0;
  for(mysize_t i = 0; i < N; i++) {
    perm[i] = i;
    for(int j = 0; j < e; j++) {
      this[j] = factors[j]->group[i];
    }
    if(hasinteract || (memcmp(this,prev,e*sizeof(int)) != 0) ) {
      int nlev=0;
      /* not duplicate, store in indices */
      
      for(int j = 0; j < e; j++) {
	indices[ie+j] = this[j]-1 + nlev;
	nlev += factors[j]->nlevels;
      }
      newR[newN] = R[i];
      newN++;
      ie += e;
      if(!hasinteract) memcpy(prev,this,e*sizeof(int));
    }
  }

  /* 
     At this point we should perhaps randomly shuffle the equations.
     We don't know when this ought to be done, but we've only seen it when
     there are many factors.
     We want to use unif_rand to be able to get reproducible results, 
     at least for single-threaded things, and keeping the random number
     generator in the same state when we return to R.
     The unif_rand isn't concurrency-safe, so protect with mutex.
     The Get/PutRNGstate is done in the main thread.
     Knuth-Fisher-Yates shuffle.
     The randomization done here should not have consequences for the result
     of the routine, only the speed.
  */
  LOCK(lock);
  if(e > 1) {
    for(mysize_t i = newN-1; i > 0; i--) {
      mysize_t j;
      /* Pick j between 0 and i inclusive */
      j = (mysize_t) floor((i+1) * unif_rand());
      if(j == i) continue;
      /* exchange newR[i] and newR[j]
	 as well as indices[i*e:i*e+e-1] and indices[j*e:j*e+e-1]
      */
      double dtmp = newR[j];
      newR[j] = newR[i];
      newR[i] = dtmp;
      mysize_t itmp = perm[j];
      perm[j] = perm[i];
      perm[i] = itmp;
      for(mysize_t k = 0; k < e; k++) {
	mysize_t itmp;
	itmp = indices[j*e+k];
	indices[j*e+k] = indices[i*e+k];
	indices[i*e+k] = itmp;
      }
    }
  }
  UNLOCK(lock);

  /* Then, do the Kaczmarz iterations */
  norm2 =0.0;
  for(mysize_t i = 0; i < newN; i++) norm2 += newR[i]*newR[i];
  norm2 = sqrt(norm2);
  prevdiff = 2*norm2;
  neweps = eps*norm2;


  do {
    mysize_t ie = 0; /* equals i*e; integer multiplication is slow, keep track instead */
    double diff = 0.0;
    for(mysize_t i = 0; i < newN; i++,ie+=e) {
      const mysize_t ip = perm[i];
      double upd = 0.0;
      double ai2 = 0.0;
      upd = newR[i];
      
      /* Subtract inner product */
      for(int j = 0; j < e; j++) {
	const mysize_t idx = indices[ie + j];
	const double *fx = factors[j]->x;
	const double w = (fx == NULL) ? 1.0 : fx[ip];
	upd -= x[idx]*w;
	ai2 += w*w;
      }
      /* Update */
      upd /= ai2;
      for(int j = 0; j < e; j++) {
	const mysize_t idx = indices[ie + j];
	const double *fx = factors[j]->x;
	const double w = (fx == NULL) ? upd : upd*fx[ip];
	x[idx] += w;
	diff += w*w;
      }
    }
    iter++;

    sd = sqrt(diff);
    c = sd/prevdiff;
    prevdiff = sd;
    if(c >= 1.0 && iter > 20) {
      char buf[256];
      sprintf(buf,"Kaczmarz failed in iter %d*%d, c=1-%.0e, delta=%.1e, eps=%.1e\n",iter,newN,1.0-c,sd,neweps);
      pushmsg(buf,lock);
      break;
    }
#ifdef NOTHREADS
    R_CheckUserInterrupt();
#else
    if(*stop != 0) return(0);
#endif
  } while(sd >= neweps*(1.0-c) && neweps > 1e-15);

  return(sd);
}

#ifdef WIN
DWORD WINAPI kaczmarz_thr(LPVOID varg) {
#else
static void *kaczmarz_thr(void *varg) {
#endif
  KARG *arg = (KARG*) varg;
  int myid;
  int vecnum;
  /* Get a thread id */
  LOCK(arg->lock);
  myid = arg->threadnum++;
  UNLOCK(arg->lock);
  while(1) {
    LOCK(arg->lock);
    vecnum = arg->nowdoing++;
    UNLOCK(arg->lock);

    if(vecnum >= arg->numvec) break;
#ifdef HAVE_THREADNAME
    char thrname[16];
    snprintf(thrname, 16, "Kz %d/%d",vecnum+1, arg->numvec);
    STNAME(thrname);
#endif

    (void) kaczmarz(arg->factors,arg->e,arg->N,
		    arg->source[vecnum],arg->target[vecnum],
		    arg->eps, arg->work[myid], &arg->stop, arg->lock);
  }
#ifndef NOTHREADS
#ifndef WIN
  LOCK(arg->lock);
  arg->running--;
#ifdef HAVE_SEM
  if(arg->running == 0) sem_post(&arg->finished);
#endif
  UNLOCK(arg->lock);
#endif
#endif
  return 0;
}

SEXP MY_kaczmarz(SEXP flist, SEXP vlist, SEXP Reps, SEXP initial, SEXP Rcores) {

  double eps = REAL(Reps)[0];
  /*  double *R = REAL(RR);*/
  FACTOR **factors;
  double *init = 0;
  int cores = INTEGER(Rcores)[0];
  int numfac;
  mysize_t N = 0;
  //  int i;
  int sumlev = 0;
  int listlen;
  int numvec;
  SEXP reslist;
  KARG arg;
  double **vectors, **target;
  int cnt;
  int numthr = 1;
  int protectcount=0;
#ifndef NOTHREADS
  int thr;
#ifdef WIN
  HANDLE *threads;
  DWORD *threadids;
  HANDLE lock;
#else
  pthread_t *threads;
  pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif
#endif

#ifndef NOTHREADS
  /* Set up the message stack */
  initmsg();
#endif
  if(!isNull(initial)) {
    init = REAL(initial);
  }
  PROTECT(flist = AS_LIST(flist));protectcount++;
  //  numfac = LENGTH(flist);

  factors = makefactors(flist, 0, NULL);
  numfac = 0;
  for(FACTOR **f = factors; *f != NULL; f++) numfac++;
  N = LENGTH(VECTOR_ELT(flist,0));
  for(int i = 0; i < numfac; i++)
    sumlev += factors[i]->nlevels;


  if(!isNull(initial) && LENGTH(initial) != sumlev)
    error("Initial vector must have length %d, but is %d\n",sumlev, LENGTH(initial));

  /* Then the vectors */

  PROTECT(vlist = AS_LIST(vlist)); protectcount++;
  listlen = LENGTH(vlist);
  PROTECT(reslist = NEW_LIST(listlen)); protectcount++;

  /* First, count the number of vectors in total */
  numvec = 0;
  for(int i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    /* Each entry in the list is either a vector or a matrix */
    if(!isMatrix(elt)) {
      if(LENGTH(elt) != N) 
	error("Vector length (%d) must be equal to factor length (%d)",LENGTH(elt),N);
      numvec++;
    } else {
      if(nrows(elt) != N)
	error("Vector length must be equal to factor length %d %d ",nrows(elt),N);
      numvec += ncols(elt);
    }
  }

  /* Allocate pointers to source vectors */
  vectors = (double **)R_alloc(numvec,sizeof(double*));
  /* Allocate pointers to result vectors */
  target = (double**) R_alloc(numvec,sizeof(double*));
  /* Loop through list again to set up result structure */
  cnt = 0;
  for(int i = 0; i < listlen; i++) {
    SEXP elt = VECTOR_ELT(vlist,i);
    if(!isReal(elt)) {
      elt = PROTECT(coerceVector(elt, REALSXP)); protectcount++;
    }
    if(!isMatrix(elt)) {
      /* It's a vector */
      SEXP resvec;
      vectors[cnt] = REAL(elt);
      PROTECT(resvec = allocVector(REALSXP,sumlev));
      target[cnt] = REAL(resvec);
      SET_VECTOR_ELT(reslist,i,resvec);
      UNPROTECT(1);
      cnt++;
    } else {
      /* It's a matrix */
      int cols = ncols(elt);
      SEXP mtx;
      /* Allocate a matrix */
      PROTECT(mtx = allocMatrix(REALSXP,sumlev,cols));
      SET_VECTOR_ELT(reslist,i,mtx);
      UNPROTECT(1);
      /* Set up pointers */
      for(int j = 0; j < cols; j++) {
	vectors[cnt] = REAL(elt) + j*N;
	target[cnt] = REAL(mtx) + j*sumlev;
	cnt++;
      }
    }
  }


  for(int cnt = 0; cnt < numvec; cnt++) {
    if(init != 0)
      for(int i = 0; i < sumlev; i++) target[cnt][i] = init[i];
    else
      for(int i = 0; i < sumlev; i++) target[cnt][i] = 0.0;
  }

  
  /* set up for threading */

  numthr = cores;
  if(numthr > numvec) numthr = numvec;
  if(numthr < 1) numthr = 1;
  GetRNGstate();
#ifndef NOTHREADS
#ifdef WIN
  lock = CreateMutex(NULL,FALSE,NULL);
  arg.lock = lock;
  threads = (HANDLE*) R_alloc(numthr,sizeof(HANDLE));
  threadids = (DWORD*) R_alloc(numthr,sizeof(DWORD));
#else
  threads = (pthread_t*) R_alloc(numthr,sizeof(pthread_t));
  arg.running = numthr;
#ifdef HAVE_SEM
  if(sem_init(&arg.finished,0,0) != 0) error("sem_init failed, errno=%d",errno);
#endif
  arg.lock = &lock;
#endif
#endif

  arg.factors = factors;
  arg.source = vectors;
  arg.target = target;
  arg.nowdoing = 0;
  arg.threadnum = 0;
  arg.e = numfac;
  arg.eps = eps;
  arg.numvec = numvec;
  arg.N = N;
  arg.stop = 0;
  arg.work = (mysize_t**) R_alloc(numthr, sizeof(mysize_t*));

  // when allocating the work, we use mysize_t, but parts of it is accessed as double which may be
  // larger. So we allocate some more (8*sizeof(mysize_t) more), and adjust the address so it's aligned on a double
  // When using it, make sure we do all the doubles first.

#ifdef NOTHREADS
  arg.work[0] = (mysize_t*) R_alloc(numfac*N + N*sizeof(double)/sizeof(mysize_t) + N + 2*numfac+8, sizeof(mysize_t));
  uintptr_t amiss = (uintptr_t) arg.work[0] % sizeof(double);
  if(amiss != 0) arg.work[0] = (mysize_t*) ((uintptr_t)arg.work[0] + sizeof(double)-amiss);
  kaczmarz_thr((void*)&arg);
#else
  /* Do it in separate threads */
  for(thr = 0; thr < numthr; thr++) {
    // allocate some thread-specific storage, we can't use R_alloc in a thread
    arg.work[thr] = (mysize_t*) R_alloc(numfac*N + N*sizeof(double)/sizeof(mysize_t) + N + 2*numfac+8, sizeof(mysize_t));
    uintptr_t amiss = (uintptr_t) arg.work[thr] % sizeof(double);
    if(amiss != 0) arg.work[thr] = (mysize_t*) ((uintptr_t)arg.work[thr] + sizeof(double)-amiss);

#ifdef WIN
    threads[thr] = CreateThread(NULL,0,kaczmarz_thr,&arg,0,&threadids[thr]);
    if(0 == threads[thr]) error("Failed to create kaczmarz thread");
#else
    int stat = pthread_create(&threads[thr],NULL,kaczmarz_thr,&arg);
    if(0 != stat) error("Failed to create kaczmarz thread, stat=%d",stat);
#endif
  }
  /* wait for completion */
  /* We want to check for interrupts regularly, and
     set a stop flag */
  while(1) {
    printmsg(arg.lock);
    if(arg.stop == 0 && checkInterrupt()) {
      REprintf("...stopping Kaczmarz threads...\n");
      arg.stop=1;
    }

#ifdef WIN
    if(WaitForMultipleObjects(numthr,threads,TRUE,3000) != WAIT_TIMEOUT) {
      for(thr = 0; thr < numthr; thr++) {
	CloseHandle(threads[thr]);
      }
      /* Print any remaining messages */
      printmsg(arg.lock);
      CloseHandle(lock);
      break;
    }
#else
    {
#ifndef HAVE_SEM
      struct timespec atmo = {0,50000000};
      /* Kludge in MacOSX because no timedwait */
      
      if(arg.stop == 0) nanosleep(&atmo,NULL);
      if(arg.stop == 1 || arg.running == 0) {
#else
      struct timespec tmo = {time(NULL)+3,0};
      if(arg.stop == 1 || sem_timedwait(&arg.finished,&tmo) == 0) {
#endif
	for(thr = 0; thr < numthr; thr++) {
	  (void)pthread_join(threads[thr], NULL);
	}
#ifdef HAVE_SEM
	sem_destroy(&arg.finished);
#endif
	/* Print any remaining messages */
	printmsg(arg.lock);
	pthread_mutex_destroy(arg.lock);
	break;
      }
    }
#endif
  }
#endif
  PutRNGstate();
  if(arg.stop == 1) error("Kaczmarz interrupted");
  UNPROTECT(protectcount);
  return(reslist);
}
