/*
 $Id: factor.c 1662 2015-03-20 15:04:22Z sgaure $
*/

#include "lfe.h"
static void invertfactor(FACTOR *f,int N) {
  int nlev = f->nlevels;
  int *curoff;
  int i;
  f->ii = (int*) R_alloc(nlev+1,sizeof(int));
  f->gpl = (int*) R_alloc(N,sizeof(int));

  memset(f->ii,0,sizeof(int)*(nlev+1));

  /* find sizes of groups */
  for(i = 0; i < N; i++) {
    int gp = f->group[i];
    if(gp < 1) error("Factors can not have missing levels");
    f->ii[gp]++;
  }

  /* cumulative */
  for(i = 1; i < nlev+1; i++) {
    f->ii[i] += f->ii[i-1];
  }

  curoff = Calloc(nlev+1,int);
  for(i = 0; i < N; i++) {
    int gp = f->group[i]-1;
    f->gpl[f->ii[gp]+curoff[gp]] = i;
    curoff[gp]++;
  }
  Free(curoff);
}

FACTOR** makefactors(SEXP flist, int allowmissing, double *weights) {
  FACTOR **factors;
  int numfac = LENGTH(flist);
  int N=0;
  int oneiter = 0;
  numfac = 0;
  for(int i = 0; i < LENGTH(flist); i++) {
    SEXP sf = VECTOR_ELT(flist,i);
    SEXP xattr = getAttrib(sf,install("x"));    
    if(isNull(xattr)) {
      numfac++;
      continue;
    } else if(LENGTH(flist) == 1) {
      SEXP ortho = getAttrib(xattr, install("ortho"));
      if(isLogical(ortho)) oneiter = LOGICAL(ortho)[0];
    }
    if(!isMatrix(xattr)) {
      numfac++;
      continue;
    }
    numfac += ncols(xattr);
  }

  if(!oneiter) {
    SEXP Roneiter = getAttrib(flist, install("oneiter"));
    if(isLogical(Roneiter)) oneiter = LOGICAL(Roneiter)[0];
  }

  factors = (FACTOR**) R_alloc(numfac+1,sizeof(FACTOR*));
  factors[numfac] = NULL;
  int truefac = 0;
  for(int i = 0; i < LENGTH(flist); i++) {
    int len;
    FACTOR *f;
    len = LENGTH(VECTOR_ELT(flist,i));
    if(i == 0) {
      N = len;
    } else if(len != N) {
      error("All factors must have the same length %d %d",len,N);
    }
    
    f = factors[truefac++] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    if(f->nlevels <= 0) error("factor %d in list has no levels\n",i+1);
    f->oneiter = oneiter;
    SEXP xattr = getAttrib(VECTOR_ELT(flist,i),install("x"));
    if(isNull(xattr)) {
      f->x = NULL;
    } else {
      if(isMatrix(xattr)) {
	if(nrows(xattr) != len) {
	  error("Factor interaction terms (%d) must have the same length (%d) as the factor",
		LENGTH(xattr),len);
	}
	truefac--;
	for(int j = 0; j < ncols(xattr); j++) {
	  FACTOR *g = factors[truefac++] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
	  g->group = f->group;
	  g->nlevels = f->nlevels;
	  g->oneiter = f->oneiter;
	  g->x = &REAL(xattr)[j*(mybigint_t)nrows(xattr)];
	}
      } else {
	if(LENGTH(xattr) != len) {
	  error("Factor interaction terms (%d) must have the same length (%d) as the factor",
		LENGTH(xattr),len);
	}
	f->x = REAL(xattr);
      }
    }
  }

  /* Make array for holding precomputed group levels 
     Now, what about entries which don't belong to a group
     I.e. NA entries, how is that handled, I wonder.
     seems to be negative.  Anyway, we fail on them. No, we don't */
  
  for(int i = 0; i < truefac; i++) {
    FACTOR *f = factors[i];
    f->gpsize = (double *)R_alloc(f->nlevels,sizeof(double));
    f->invgpsize = (double *)R_alloc(f->nlevels,sizeof(double));
    memset(f->gpsize,0,f->nlevels*sizeof(double));
    /* first count it */
    for(int j = 0; j < N; j++) {
	/* skip entries without a group, do we need that? */
	/* if(f->group[j] < 1) error("Factors can't have missing levels"); */
      if(f->group[j] > 0) {
	double w = (f->x == NULL) ? (weights==NULL ? 1.0 : weights[j]) :
	  (weights==NULL ? f->x[j] : f->x[j]*weights[j]);
	  f->gpsize[f->group[j]-1] += w*w;
      } else {
	if(!allowmissing) error("Factors can't have missing levels");
      }
    }
    /* then invert it, it's much faster to multiply than to divide */
    /* in the iterations */
    for(int j = 0; j < f->nlevels; j++) {
      f->invgpsize[j] = 1.0/f->gpsize[j];
    }
  }
  return(factors);
}


/*
  This one is a bit tricky.  We could do it quite elegantly recursively,
  but we would suffer from a deep call-stack.  Hence, we make our own
  stack and push/pop in a loop
 */
static int Components(int **vertices, FACTOR **factors, int K) {
  int stacklevel = 0;
  int *stack;
  int curfac, curvert, curcomp,candvert;
  int startfac=0,startvert=0;
  int ii=0,i;
  int numvert=0;

  /* How big a stack ? */
  /* The number of vertices */

  for(i = 0; i < K; i++) numvert += factors[i]->nlevels;

  /* Never used in threads, so use R's Calloc */
  stack = Calloc(numvert*4,int);

#define PUSH(x) stack[stacklevel++] = x
#define POP(x) x = stack[--stacklevel]
#define PUSHALL {PUSH(startvert); PUSH(startfac); PUSH(curfac); PUSH(ii);}
#define POPALL {POP(ii); POP(curfac); POP(startfac); POP(startvert);}

  curcomp = 1;
  candvert = 0;

  do {
    curvert = candvert;
    curfac = 0;
    /* Find the entire component */

    while(1) {
      /* At the top here, curfac,curvert is a candidate for marking off as 
	 a vertex in curcomp
	 For each iteration:  

	 If it's not marked, mark it, push on stack, go to first datapoint
	 If it's already marked, go to next datapoint for the vertex (incidence matrix)
 	 If data-points are exhausted, go to next factor, start over with datapoints.
	 If factors are exhausted, pop the stack
	 If final stack-frame, we're done with component.
      */
      
      if(vertices[curfac][curvert] == 0) {
	/* Mark new vertex, find incidence list */
	vertices[curfac][curvert] = curcomp;
	PUSHALL;
	startvert = curvert;
	startfac = curfac;
	curfac = (startfac+1)%K;
	ii = factors[startfac]->ii[startvert];
      } else {
	/* Been there, try next in group */
	ii++;
      }
      if(ii >= factors[startfac]->ii[startvert+1]) {
	/* No more, move to next factor */
	curfac = (curfac + 1) % K;
	if(curfac == startfac) {
	  /* This is where we began, pop */
	  /* No more neighbours, go back to previous */
	  POPALL;
	  /* Get outta here */
	  if(0 == stacklevel) break;
	} else {
	  /* start over with data-points */
	  ii = factors[startfac]->ii[startvert];	
	}
      }
      curvert = factors[curfac]->group[factors[startfac]->gpl[ii]]-1;
    }

    /* Find next component */
    while(candvert < factors[0]->nlevels && vertices[0][candvert] != 0) candvert++;
    curcomp++;
  } while(candvert < factors[0]->nlevels);

  Free(stack);
  return(curcomp-1);
#undef PUSH
#undef POP
#undef PUSHALL
#undef POPALL
}
 
// Algorithm from:
// "A Note on the Determination of Connectedness in an N-Way Cross Classification"
// D.L. Weeks and D.R. Williams, Technometrics, vol 6 no 3, August 1964
// There probably exists faster algorithms, this one is quite slow
static void wwcomp(FACTOR *factors[], int numfac, int N, int *newlevels) {
   int level = 0;
   int chead = 0;
   int *newlist = Calloc(N,int);
   int *oldlist = Calloc(N,int);
   int oldstart = 0;
   // For cache-efficiency, make a numfac x N matrix of the factors
   int *facmat = Calloc(numfac*N, int);

   for(mysize_t i = 0; i < N; i++) {
     int *obsp = &facmat[i*numfac];
     newlevels[i] = 0;
     oldlist[i] = i;
     for(int j = 0; j < numfac; j++) {
       obsp[j] = factors[j]->group[i];
     }
   }
   while(oldstart < N) {
     int newidx;
     // set component number for first node in component
     // find the next we haven't checked, it's oldlist[oldstart]
     // increase oldstart by one
     level++;
     chead = oldlist[oldstart++];
     // put it as the first element in our newlist
     newlist[0] = chead;
     newlevels[chead] = level;
     newidx = 1;
     // loop over the list of newly added nodes, including the head
     // Note that we may increase newidx during the loop
     for(int i = 0; i < newidx; i++) {
       mysize_t newnode = newlist[i];
       int *newp = &facmat[newnode*numfac];
       // search for observations with distance 1 from newnode, mark them with level
       for(int jidx = oldstart; jidx < N; jidx++) { 
	 mysize_t trynode = oldlist[jidx];
	 int *tryp = &facmat[trynode*numfac];
	 int dist = 0;
	 // compute distance
	 for(int fi = 0; fi < numfac && dist < 2; fi++) 
	   dist += (newp[fi] != tryp[fi]);
	 //dist += (factors[fi]->group[newnode] != factors[fi]->group[trynode]);
	 // if close, set its level, add it to the list, move the start node
	 // to the empty place in oldlist. 
	 if(dist < 2) {
	   newlevels[trynode] = level;
	   newlist[newidx++] = trynode;
	   oldlist[jidx] = oldlist[oldstart++];
	 }
       }
     }
   }
   Free(facmat);
   Free(newlist);
   Free(oldlist);
 }

/*
R entry-point for conncomp.  Takes a list of factors as input.
 */

SEXP MY_wwcomp(SEXP flist) {
   int numfac, N;
   FACTOR **factors;
   SEXP result;

  
  numfac = LENGTH(flist);
  if(numfac < 2) error("At least two factors must be specified");

  N = LENGTH(VECTOR_ELT(flist,0));
  for(int i = 0; i < numfac; i++) {
    if(N != LENGTH(VECTOR_ELT(flist,i))) 
      error("Factors must have the same length");
  }

  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  for(int i = 0; i < numfac; i++) {
    FACTOR *f;

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
  }
  PROTECT(result = allocVector(INTSXP,N));
  int *fac = INTEGER(result);
  wwcomp(factors, numfac, N, fac);
  // Now it's time to order the levels by decreasing size, so let's compute the sizes
  int levels = 0;
  for(int i = 0; i < N; i++) if(fac[i] > levels) levels = fac[i];
  double *levsize = (double*) R_alloc(levels, sizeof(double));
  int *index = (int*) R_alloc(levels, sizeof(int));
  for(int i = 0; i < levels; i++) {
    levsize[i] = 0.0;
    index[i] = i;
  }
  for(int i = 0; i < N; i++) levsize[fac[i]-1] = levsize[fac[i]-1]+1;
  revsort(levsize,index,levels);
  int *rindex = (int*) R_alloc(levels, sizeof(int));
  for(int i = 0; i < levels; i++) rindex[index[i]] = i;
  for(int i = 0; i < N; i++) {
    fac[i] = rindex[fac[i]-1]+1;
  }

  
  UNPROTECT(1);
  return result;
}
/* 
Then for finding connection components 
From R we take a list of factors, we return
a factor of the same length with the connection
components
*/

SEXP MY_conncomp(SEXP flist) {
  int numfac;
  int i;
  
  int N;
  FACTOR **factors;
  int *group;
  int **vertices;
  SEXP result;
  int *resgroup;
  int comps;
  int *levtrl;
  double *gpsiz;
  int *idx;

  numfac = LENGTH(flist);
  if(numfac < 2) error("At least two factors must be specified");
  N = LENGTH(VECTOR_ELT(flist,0));
  for(i = 0; i < numfac; i++) {
    if(N != LENGTH(VECTOR_ELT(flist,i))) 
      error("Factors must have the same length");
  }
  factors = (FACTOR**) R_alloc(numfac,sizeof(FACTOR*));
  PROTECT(flist = AS_LIST(flist));
  for(i = 0; i < numfac; i++) {
    FACTOR *f;

    factors[i] = (FACTOR*) R_alloc(1,sizeof(FACTOR));
    f = factors[i];
    f->group = INTEGER(VECTOR_ELT(flist,i));
    f->nlevels = LENGTH(getAttrib(VECTOR_ELT(flist,i),R_LevelsSymbol));
    if(f->nlevels == 0) error("Factor %s has zero levels", CHAR(STRING_ELT(GET_NAMES(flist),i)));
    invertfactor(f,N);
  }

  /* Create vertices */
  vertices = (int**) R_alloc(numfac,sizeof(int*));
  /* Create arrays for them */
  for(i = 0; i < numfac; i++) {
    vertices[i] = (int*) R_alloc(factors[i]->nlevels,sizeof(int));
    /* Assign no component to them*/
    memset(vertices[i],0,sizeof(int)*factors[i]->nlevels);
  }

  /* Do the stuff */
  comps = Components(vertices,factors,numfac);

  /* allocate result structure */
  PROTECT(result = allocVector(INTSXP,N));
  resgroup = INTEGER(result);
  group = factors[0]->group;
  for(i = 0; i < N; i++) {
    resgroup[i] = vertices[0][group[i]-1];
  }

  /* the levels should be ordered by decreasing size. How do we do this? 
     Hmm, we should have a look at revsort supplied by R.
     There must be an easier way, I'm clumsy today.   
  */

  gpsiz = Calloc(comps,double);
  idx = Calloc(comps,int);
  
  for(i = 0; i < comps; i++) idx[i] = i;
  for(i = 0; i < N; i++) {
    gpsiz[resgroup[i]-1]++;
  }

  revsort(gpsiz,idx,comps);
  Free(gpsiz);
  levtrl = Calloc(comps,int);
  for(i = 0; i < comps; i++) levtrl[idx[i]] = i+1;
  Free(idx);

  for(i = 0; i < N; i++) {
    resgroup[i] = levtrl[resgroup[i]-1];
  }
  Free(levtrl);

  UNPROTECT(2);
  return(result);
}
