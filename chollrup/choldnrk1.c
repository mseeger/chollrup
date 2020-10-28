/* -------------------------------------------------------------------
 * CHOLDNRK1
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * If
 *   A = L*L', A_ = A - v*v' = L_ L_', A, L n-by-n,
 * where L is lower triangular, the method computes L_ from L. This
 * is called Cholesky rank one downdate. L or L' (upper triangular) can
 * be passed, only the relevant triagle is accessed. L (or L') is
 * passed in L, v in VEC.
 * We require p = L\v. If ISP==true, VEC contains p rather than v.
 * Otherwise, p is computed locally, stored in WORKV.
 * NOTE: For the present implementation, the method is more efficient
 * when a lower triangular matrix is used.
 *
 * Dragging along:
 * If Z (r-by-n) is given, so must be the r-vector y. In this case,
 * we overwrite Z by Z_, where
 *   Z_ L_' = Z L' - y v'.
 *
 * Working array:
 * Requires a working vector of size >= max(n,r), passed in WORKV.
 * NOTE: The same vector can be passed for VEC and WORKV, in which case
 * VEC is overwritten in an undefined way. If r > n, VEC can be of size
 * r, containing v (or p) in the first n components.
 * The method uses n Givens rotations, param. by angles c_k, s_k. The
 * change L -> L' is specified by these 2*n numbers. They are written
 * into CVEC, SVEC.
 * NOTE: The original routine LINPACK dchdd can in theory produce negative
 * values in DIAG(L), something we want to avoid here. If this happens,
 * the corr. column of L_ is flipped. In the present implementation, this
 * is not reported back, so the change L -> L' cannot always be
 * reconstructed from CVEC, SVEC alone.
 *
 * Input:
 * - L:     Factor L (or L'), overwritten by L_ (or L_'). Must be
 *          lower (upper) triangular, str. code UPLO
 * - VEC:   Vector v. Can have size >n, only first n elem. are used
 * - CVEC:  Vector [n]. c_k ret. here
 * - SVEC:  Vector [n]. s_k ret. here
 * - WORKV: Working vector of size max(n,r). Can be same as VEC
 * - ISP:   S.a. Def.: false
 * - Z:     Dragging along matrix [r-by-n]. Optional
 * - Y:     Dragging along vector [r]. Iff Z is given
 *
 * Return:
 * - STAT:  0 (OK), 1 (Numerical error)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#include <math.h>
#include "mex.h"
#include "mex_helper.h"
#include "blas_headers.h"

char errMsg[200];

/*
 * The method is adapted from LINPACK dchdd. We did the following
 * modifications:
 * - Using BLAS drot in order to avoid any explicit O(n^2) loops
 * - Keeping diag(L_) positive, by flipping columns of L_ whenever a
 *   negative element pops up
 * See the TR
 *   M. Seeger
 *   Low Rank Updates for the Cholesky Decomposition
 *   Available at: www.kyb.tuebingen.mpg.de/bs/people/seeger/papers/
 */

/* Main function CHOLDNRK1 */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,n,r=0,stp,sz,ione=1,retcode=0,nxi,npos;
  double qs,cval,sval,c1,c2;
  fst_matrix lmat,zmat;
  const double* vvec;
  double* cvec,*svec,*wkvec,*yvec;
  bool islower,isp=false;
  double* tbuff,*zcol;
  const char* diag="N";
  char trans[2];
  int* flind=0;

  /* Read arguments */
  if (nrhs<5)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs>1)
    mexErrMsgTxt("Too many return arguments");
  parseBLASMatrix(prhs[0],"L",&lmat,-1,-1);
  if ((n=lmat.n)!=lmat.m ||
      (!(islower=(UPLO(lmat.strcode)=='L')) && UPLO(lmat.strcode)!='U'))
    mexErrMsgTxt("L must be lower/upper triangular (use UPLO str. code!)");
  if (getVecLen(prhs[1],"VEC")<n) mexErrMsgTxt("VEC too short");
  vvec=mxGetPr(prhs[1]);
  if (getVecLen(prhs[2],"CVEC")!=n || getVecLen(prhs[3],"SVEC")!=n)
    mexErrMsgTxt("CVEC, SVEC have wrong size");
  cvec=mxGetPr(prhs[2]); svec=mxGetPr(prhs[3]);
  if (nrhs>5) {
    isp=(getScalInt(prhs[5],"ISP")!=0);
    if (nrhs>6) {
      if (nrhs<8) mexErrMsgTxt("Need both Z, Y");
      parseBLASMatrix(prhs[6],"Z",&zmat,-1,n);
      r=zmat.m;
      if (getVecLen(prhs[7],"Y")!=r) mexErrMsgTxt("Y has wrong size");
      yvec=mxGetPr(prhs[7]);
    }
  }
  i=getVecLen(prhs[4],"WORKV");
  if (i<n || i<r) mexErrMsgTxt("WORKV too short");
  wkvec=mxGetPr(prhs[4]);

  /* Compute p (if not given) */
  BLASFUNC(dcopy) (&n,vvec,&ione,wkvec,&ione);
  if (!isp) {
    trans[1]=0;
    trans[0]=islower?'N':'T';
    BLASFUNC(dtrsv) (&UPLO(lmat.strcode),trans,diag,&n,lmat.buff,&lmat.stride,
		     wkvec,&ione);
  }
  /* Generate Givens rotations */
  qs=1.0-BLASFUNC(ddot) (&n,wkvec,&ione,wkvec,&ione);
  if (qs<=0.0)
    retcode=1;
  else {
    qs=sqrt(qs);
    for (i=n-1; i>=0; i--) {
      BLASFUNC(drotg) (&qs,wkvec+i,cvec+i,svec+i);
      /* 'qs' must remain positive */
      if (qs<0.0) {
	qs=-qs; cvec[i]=-cvec[i]; svec[i]=-svec[i];
      }
    }
    /* NOTE: 'qs' should be 1 now */

    /* Update L. If there are any flips of L_ cols, we alloc. 'flind' are
       store their pos. there */
    fillVec(wkvec,n,0.0);
    stp=islower?1:lmat.stride;
    for (i=n-1,sz=0,tbuff=lmat.buff+((n-1)*(lmat.stride+1)); i>=0; i--) {
      /* BAD: Slower for upper triangular! */
      sz++;
      if (*tbuff<=0.0) {
	retcode=1; break;
      }
      BLASFUNC(drot) (&sz,wkvec+i,&ione,tbuff,&stp,cvec+i,svec+i);
      /* Do not want negative elements on diagonal */
      if (*tbuff<0.0) {
	if (flind==0) {
	  /* Does this ever happen?
	     Allocate 'flind'. Size n, to make sure. Grows from the right */
	  flind=(int*) mxMalloc(n*sizeof(int));
	  npos=n;
	}
	flind[--npos]=i;
	qs=-1.0;
	BLASFUNC(dscal) (&sz,&qs,tbuff,&stp);
      } else if (*tbuff==0.0) {
	retcode=1; break;
      }
      tbuff-=(lmat.stride+1);
    }
    /* NOTE: Should have v in 'wkvec' now */

    /* Dragging along */
    if (r>0 && retcode==0) {
      BLASFUNC(dcopy) (&r,yvec,&ione,wkvec,&ione);
      nxi=(flind!=0)?flind[npos]:-1;
      for (i=0; i<n; i++) {
	zcol=zmat.buff+(i*zmat.stride);
	cval=cvec[i]; sval=svec[i];
	qs=-sval;
	BLASFUNC(daxpy) (&r,&qs,wkvec,&ione,zcol,&ione);
	if (nxi==i) {
	  if (++npos<n) nxi=flind[npos];
	  c1=-1.0/cval; c2=sval;
	} else {
	  c1=1.0/cval; c2=-sval;
	}
	BLASFUNC(dscal) (&r,&c1,zcol,&ione);
	if (i<n-1) {
	  BLASFUNC(dscal) (&r,&cval,wkvec,&ione);
	  BLASFUNC(daxpy) (&r,&c2,zcol,&ione,wkvec,&ione);
	}
      }
    }
  }

  if (flind!=0) mxFree((void*) flind);
  if (nlhs==1) {
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[0]))=(double) retcode;
  }
}
