/* -------------------------------------------------------------------
 * CHOLUPEXCH
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * Update of Cholesky factor after special row/column permutation
 * of system matrix.
 *   A = R' R, A_ = E' A E, E spec. permut. matrix
 * Here, R is upper triangular. The factor must be passed in R and is
 * overwritten by R_. The method works by computing an orthonormal
 * U s.t. U R E = R'.
 * E is determined by K, L, JOB (integers). Here, 1 <= K < L <= n
 * (size of A). E reorders rows/columns as follows:
 * If JOB==1, the new ordering is ...,K-1,L,K,...,L-1,L+1,... In this
 * case, U = U(L-K)*...*U(1), where U(i) is a Givens rotation in the
 * plane (L-i,L-i+1). The Givens form is [c(i) s(i); -s(i) c(i)].
 * If JOB==2, the new ordering is ...,K-1,K+1,...,L,K,L+1,... In this
 * case, U = U(L-K)*...*U(1), where U(i) is a Givens rotation in the
 * plane (K+i-1,K+i).
 *
 * Dragging along: If X is passed, it is replaced by X_ = U X. Note
 * that if R' X = B, then (R_)' X_ = E' B
 *
 * The method is just a wrapper around the LINPACK routine DCHEX. For
 * more information, see
 * @book{Dongarra:79,
 *   author      = {Dongarra, J. and Moler, C. and Bunch, J. and Stewart, G.},
 *   title       = {{LINPACK} User's Guide},
 *   publisher   = {SIAM Society for Industrial and Applied Mathematics},
 *   year        = {1979}
 * }
 * NOTE: There is a bug in DCHEX, causing elements of DIAG(R_) to be
 * negative. We include a workaround here.
 *
 * Input:
 * - R:     Factor R, overwritten by R_
 * - K:     Describes E (s.a.)
 * - L:     Describes E (s.a.)
 * - JOB:   Describes E (s.a.)
 * - X:     Drag-along matrix, replaced by X_ [def: []]
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include <math.h>
#include "mex.h"
#include "mex_helper.h" /* Helper functions */
#include "blas_headers.h"

char errMsg[200];

/* LINPACK DCHEX declaration */
extern void BLASFUNC(dchex) (double* r,int* ldr,int* p,int* k,int* l,
			     double* z,int* ldz,int* nz,double* c,double* s,
			     int* job);

/* Main function CHOLUPEXCH */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int n,k,l,job,nz=0,farg1,j;
  double* rfact,*xmat=0,*cvec,*svec;
  double temp;

  /* Read arguments */
  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  if (!mxIsDouble(prhs[0]) || (n=mxGetM(prhs[0]))!=mxGetN(prhs[0]))
    mexErrMsgTxt("Wrong argument R");
  rfact=mxGetPr(prhs[0]);
  if ((k=getScalInt(prhs[1],"K"))<1)
    mexErrMsgTxt("Wrong argument K");
  l=getScalInt(prhs[2],"L");
  if (l<=k || l>n)
    mexErrMsgTxt("Wrong argument L");
  job=getScalInt(prhs[3],"JOB");
  if (job<1 || job>2)
    mexErrMsgTxt("Wrong argument JOB");
  if (nrhs>4) {
    if (!mxIsEmpty(prhs[4])) {
      if (!mxIsDouble(prhs[4]) || mxGetM(prhs[4])!=n)
	    mexErrMsgTxt("Wrong argument X");
      nz=mxGetN(prhs[4]); xmat=mxGetPr(prhs[4]);
    }
  }
  cvec=(double*) mxMalloc(n*sizeof(double));
  svec=(double*) mxMalloc(n*sizeof(double));

  /* Call DCHEX */
  BLASFUNC(dchex) (rfact,&n,&n,&k,&l,xmat,&n,&nz,cvec,svec,&job);
  /* There is a strange bug in DCHEX. In some cases,
     R(j,j) < 0 for some k<=j<=l, the whole corr. row has to be multiplied
     by -1 to get the correct Cholesky factor. Here is a workaround. */
  for (j=k; j<=l; j++)
    if (rfact[(j-1)*(n+1)]<0) {
      farg1=n-j+1; temp=-1.0;
      BLASFUNC(dscal) (&farg1,&temp,rfact+((j-1)*(n+1)),&n);
      if (nz!=0)
	BLASFUNC(dscal) (&nz,&temp,xmat+(j-1),&n);
    }

  /* Deallocate */
  mxFree((void*) cvec); mxFree((void*) svec);
}
