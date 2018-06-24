#include "lfe.h"
/*
 $Id$

 Parallel Pothen-Fan algorithm
 as described in 
 "On parallel push–relabel based algorithms for bipartite maximum matching",
 Langguth et al, Parallel Computing 40(7):289-308, 2014
 http://dx.doi.org/10.1016/j.parco.2014.03.004
*/

typedef struct {
  int degree;
  int tdeg;
  int visited;
  int isR;
  int lookahead;
  int *nlist;  // pointer to start of neighbor list
} VERTEX;

typedef struct {
  int u,v;  // indices into a VERTEX array, u is row, v is columns
} PAIR;


// A matching.  I.e. a set of paths. 
typedef struct {
  // Indexed by row, value is column
  int *pairs;
  // Indexed by column, value is column
  int *ipairs;
  VERTEX *vlist;
  LOCK_T lock;
  int Nr,Nc;
} MATCHING;

typedef struct {
  int *stack;
  int sp;
} STACK;
#define PUSH(S,i) {S.stack[S.sp++] = i;}
#define POP(S) S.stack[--S.sp]
#define TOP(S) (S.sp == 0)


R_INLINE static void addmate(MATCHING M, int u, int v) {
  if(M.vlist[u].isR) {
    M.pairs[u] = v;
    M.ipairs[v-M.Nr] = u;
  } else {
    M.pairs[v] = u;
    M.ipairs[u-M.Nr] = v;
  }
}

R_INLINE static void rmmate(MATCHING M, int u) {
  if(M.vlist[u].isR) {M.ipairs[M.pairs[u]-M.Nr] = -1; M.pairs[u] = -1; return;}
  M.pairs[M.ipairs[u-M.Nr]] = -1; M.ipairs[u-M.Nr] = -1;
}
R_INLINE static int mate(MATCHING M, int u) {
  if(M.vlist[u].isR) return(M.pairs[u]);
  return(M.ipairs[u-M.Nr]);
}


R_INLINE static int isMatched(MATCHING M, int u) {
  if(M.vlist[u].isR) return(M.pairs[u] != -1);
  return(M.ipairs[u-M.Nr] != -1);
}

R_INLINE static int matched(MATCHING M, int u, int v) {
  if(M.vlist[u].isR) return(M.pairs[u] == v);
  return M.pairs[v] == u;
}

R_INLINE static int atomicadd(int *u, int v) {
  int ret;
  //#pragma omp critical
  ret = *u;
  (*u) += v;
  //#pragma end critical
  return ret;
}
R_INLINE static int testandset(int *u) {
  int v;
  //#pragma omp critical
  v = *u;
  *u = 1;
  //#pragma omp end critical
  return(v);
}

static void MatchAndUpdate(VERTEX *vlist, MATCHING M, int v) {
  if(atomicadd(&vlist[v].visited,1) != 0) return;
  for(int i = 0; i < vlist[v].degree; i++) {
    int u = vlist[v].nlist[i];
    if(atomicadd(&vlist[u].visited,1) != 0) continue;
    addmate(M,v,u);
    for(int j = 0; j < vlist[u].degree; j++) {
      int w = vlist[u].nlist[j];
      if(atomicadd(&vlist[w].tdeg, -1) == 2) {
	MatchAndUpdate(vlist, M, w);
      }
    }
    break;
  }
}


static void InitMatch(VERTEX *vlist, MATCHING M, int N) {
  int *Q1, *Qo;
  int Q1len = 0, Qolen = 0;
  for(int i = 0; i < N; i++) {
    vlist[i].visited = 0;
    vlist[i].tdeg = vlist[i].degree;
    if(vlist[i].isR) continue;
    if(vlist[i].degree == 1)
      Q1len++;
    else
      Qolen++;
  }
  Q1 = (int *) R_alloc(Q1len+1, sizeof(int));
  Qo = (int *) R_alloc(Qolen+1, sizeof(int));
  int q1 = 0, qo = 0;
  for(int i = 0; i < N; i++) {
    if(vlist[i].isR) continue;
    if(vlist[i].degree == 1) Q1[q1++] = i; else Qo[qo++] = i;
  }
  Q1[Q1len] = -1;
  Qo[Qolen] = -1;
  for(int *u = Q1; *u != -1; u++) {
    MatchAndUpdate(vlist, M, *u);
  }
  for(int *v = Qo; *v != -1; v++) {
    MatchAndUpdate(vlist, M, *v);
  }
}


static int DFS(VERTEX *vlist, MATCHING M, int v, STACK stack) {
  int idx;
  while(1) {
    //  top:
    for(int i = vlist[v].lookahead; i < vlist[v].degree; i++) {
      int u = vlist[v].nlist[i];
      vlist[v].lookahead++;
      if(!isMatched(M, u)) {
	if(testandset(&vlist[u].visited) == 0) {
	  if(TOP(stack)) {PUSH(stack,v); PUSH(stack,0);}
	  PUSH(stack,u); PUSH(stack,0);
	  // Here the path is on the stack, return success
	  return 1;
	}
      }
    }
    idx = 0;
    while(1) {
      int u = vlist[v].nlist[idx];
      if(testandset(&vlist[u].visited) == 0) {
	PUSH(stack,v); PUSH(stack,idx);
	v = mate(M,u);
	break;
	//	goto top;
      }
      // next neighbour
      while(++idx == vlist[v].degree) {
	// failed, pop the stack until non-fail
	if(TOP(stack)) return 0;
	idx = POP(stack); v = POP(stack);
      }
    }
  }
}

// ppf takes an array of edges as input
// vertices are numbered from 0 to N-1, where N is the number of vertices
MATCHING ppf(PAIR *edges, int m) {
  int nthreads = 1; 
  int Nr = 0, Nc = 0;
  // Find the number of rows and columns
  for(int i = 0; i < m; i++) {
    if(edges[i].u > Nr) Nr = edges[i].u;
    if(edges[i].v > Nc) Nc = edges[i].v;
  }
  Nr++; Nc++;
  int N = Nr+Nc;
  // Allocate the vertex array
  VERTEX *vlist = (VERTEX*) R_alloc(N, sizeof(VERTEX));

  // Fill in the vertex list
  for(int i = 0; i < N; i++) {
    vlist[i].degree = 0;
  }
  // fill in the degrees
  for(int i = 0; i < m; i++) {
    vlist[edges[i].u].degree++;
    vlist[Nr+edges[i].v].degree++;
  }
  // find total number of neighbours
  int totalnb = 0;
  for(int i = 0; i < N; i++) {
    totalnb += vlist[i].degree;
  }
  // allocate neighbour list
  int *nlist = (int *) R_alloc(totalnb, sizeof(int));
  // fill in a pointer in each vertex
  int *curnlist = nlist;
  for(int i = 0; i < N; i++) {
    vlist[i].nlist = curnlist;
    vlist[i].tdeg = 0;
    curnlist += vlist[i].degree;
  }
  // fill in neighbour list for each vertex
  for(int i = 0; i < m; i++) {
    int R = edges[i].u;
    int C = Nr+edges[i].v;
    vlist[R].nlist[vlist[R].tdeg++] = C;
    vlist[C].nlist[vlist[C].tdeg++] = R;
    vlist[R].isR = 1;
    vlist[C].isR = 0;
  }
  /*
  for(int i = 0; i < Nr; i++) {
    Rprintf("vert %d/%d: R=%d, deg=%d\n",i,Nr-1,vlist[i].isR,vlist[i].degree);
    Rprintf("neighbours: ");
    for(int j = 0; j < vlist[i].degree; j++) Rprintf("%d ",vlist[i].nlist[j]);
    Rprintf("\n");
  }
  */

  // There, we have filled in the VERTEX list.
  // We need to allocate the matching structure.
  // It can't be larger than min(Nr,Nc)
  MATCHING match;

  /*
#ifdef WIN
  match.lock = CreateMutex(NULL,FALSE,NULL);
  if(match.lock == NULL) {
    error("Couldn't create mutex (error=%d)", GetLastError());
  }
  match.lock = lock;
#else
  pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
  match.lock = &lock;
#endif
  */
  match.pairs = (int*) R_alloc(Nr, sizeof(int));
  match.Nr = Nr;
  match.Nc = Nc;
  for(int i = 0; i < Nr; i++) match.pairs[i] = -1;
  match.ipairs = (int*) R_alloc(Nc, sizeof(int));
  for(int i = 0; i < Nc; i++) match.ipairs[i] = -1;
  match.vlist = vlist;
  InitMatch(vlist, match, N);
  for(int i = Nr; i < N; i++) {
    vlist[i].lookahead = 0;
  }

  /*
  for(int i=0; i < Nr; i++) {
    Rprintf("match %d %d\n", i, match.pairs[i]);
  }
  for(int i=0; i < Nc; i++) {
    Rprintf("imatch %d %d\n", i+Nr, match.ipairs[i]);
  }
  */
  int pathfound;
  STACK *SL;
  int ss;

  SL = (STACK*) R_alloc(nthreads,sizeof(STACK));
  ss = (Nr < Nc)?( (m < Nr)?m:Nr ): Nc;
  for(int i = 0; i < nthreads; i++)  SL[i].stack = (int *) R_alloc(2*ss,sizeof(int));
  do {
    pathfound = 0;
    for(int u = 0; u < Nr; u++) {
      vlist[u].visited = 0;
    }
    // This could be done in parallel

    //    #pragma omp parallel for num_threads(nthreads) reduction(|:pathfound)

    for(int v = Nr; v < N; v++) { 
      //      int thread = omp_get_thread_num();
      //      STACK ST = SL[thread];
      STACK ST = SL[0];
      if(!isMatched(match, v)) continue;
      ST.sp = 0;
      int ok = DFS(vlist, match, v, ST);
      if(ok && ST.sp > 0) {
	pathfound = 1;
	/* augment match with the path in the stack 
	   Edges in ST which occurs in M should be removed from M
	   Other edges in ST should be added to M
	*/
	//	LOCK(match.lock);  // hmm, we don't need the lock, I think
	int u,v;
	(void)POP(ST); u = POP(ST);
	while(!TOP(ST)) {
	  (void)POP(ST); v = POP(ST);
	  if(matched(match,u,v))
	    rmmate(match,u);
	  else
	    addmate(match,u,v);
	  u = v;
	};
	//	UNLOCK(match.lock);
      }
    }
  } while(pathfound != 0);
  return match;
}


static void markresidual(MATCHING match) {
  VERTEX *vlist = match.vlist;
  STACK stack;
  int N = match.Nr + match.Nc;
  stack.sp = 0;
  stack.stack = (int *) R_alloc(2*N,sizeof(int));
  for(int i = 0; i < N; i++) vlist[i].visited = 0;
  for(int i = 0; i < match.Nr; i++) {
    // These are edges from the virtual source to the row partition.
    // Use edge only if there's capacity. I.e. if it's matched it must have degree > 1
    // If it's unmatched it must have degree > 0
    if(vlist[i].degree == 0 || (vlist[i].degree == 1 && isMatched(match,i))) continue;
    if(testandset(&vlist[i].visited)) continue;
    Rprintf("mr i=%d\n", i);
    // loop through the neigbours 
    int u = i;
    int idx = 0;
    while(1) {
      int v = vlist[u].nlist[idx];
      // Only consider it if it's not matched

      if(matched(match,u,v) || testandset(&vlist[v].visited)) {
	while(++idx >= vlist[u].degree) {
	  if(TOP(stack)) goto done;
	  idx = POP(stack);  u = POP(stack);
	  Rprintf("pop u=%d idx=%d\n",u,idx);
	}
      } else {
	Rprintf("mr v=%d (from %d)\n", v, u);
	//	Rprintf("push u=%d idx=%d\n",u,idx);
	PUSH(stack, u); PUSH(stack,idx);
	u = v; idx = 0;
      }
    }
  done:;
  }
}

SEXP MY_ppf(SEXP flist, SEXP Rtype) {
  if(LENGTH(Rtype) != 1) error("Rtype should have length 1");

  int type = INTEGER(AS_INTEGER(Rtype))[0];
  if(type < 1 || type > 2) error("Invalid type");
  if(LENGTH(flist) != 2) error("list of factors must have length 2");
  int N = LENGTH(VECTOR_ELT(flist,0));
  if(N != LENGTH(VECTOR_ELT(flist,1)))
    error("Factors must have the same length");
  int *f1 = INTEGER(VECTOR_ELT(flist,0));
  int *f2 = INTEGER(VECTOR_ELT(flist,1));
  // We should remove duplicates. Hmm, no
  PAIR *edges = (PAIR*) R_alloc(N,sizeof(PAIR));
  for(int i = 0; i < N; i++) {
    if(f1[i] <= 0 || f2[i] <= 0) error("Factors should not have missing levels");
    edges[i].u = f1[i]-1;
    edges[i].v = f2[i]-1;
  }
  MATCHING M = ppf(edges,N);

  // We return a 2 x K matrix, where K is the number of edges
  // Each column is an edge
  // Find size of matching K
  SEXP ret;
  if(type == 2) {
    // return matrix of minimum edge cut
    
    // Find the residual graph. Edges between a residual vertex and a matched vertex
    // is the cut set
    Rprintf("mark resid\n");
    markresidual(M);
    Rprintf("done resid\n");
    MATCHING mincut;
    mincut.pairs = (int*) R_alloc(M.Nr, sizeof(int));
    mincut.Nr = M.Nr;
    for(int i = 0; i < M.Nr; i++) mincut.pairs[i] = -1;
    mincut.ipairs = (int*) R_alloc(M.Nc, sizeof(int));
    for(int i = 0; i < M.Nc; i++) mincut.ipairs[i] = -1;
    mincut.vlist = M.vlist;

    // Loop through the edges. Those with one end in the residual graph and one outside
    // are the ones to use
    for(int i = 0; i < N; i++) {
      if(M.vlist[edges[i].u].visited ^ M.vlist[edges[i].v+M.Nr].visited) {
	addmate(mincut,edges[i].u,edges[i].v+M.Nr);
      }
    }
    M = mincut;
    type = 1;
  }

  int K = 0;
  for(int i = 0; i < M.Nr; i++) {
    K += (M.pairs[i] >= 0);
  }
  if(type == 1) {
    // return matrix of pairs
    PROTECT(ret = allocMatrix(INTSXP, 2, K));
    int *mat = INTEGER(ret);
    K = 0;
    for(int u = 0; u < M.Nr; u++) {
      int v = M.pairs[u];
      v -= M.Nr;
      if(v < 0) continue;
      mat[2*K] = u+1;
      mat[2*K + 1] = v+1;
      K++;
    }
    UNPROTECT(1);
  } 
  return ret;
}
