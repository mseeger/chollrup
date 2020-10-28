/* -------------------------------------------------------------------
 * CHOLUPRK1
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * If
 *   A = L*L', A_ = A + v*v' = L_ L_', A, L n-by-n,
 * where L is lower triangular, the method computes L_ from L. This
 * is called Cholesky rank one update. L or L' (upper triangular) can
 * be passed, only the relevant triagle is accessed. L (or L') is
 * passed in L, v in VEC.
 * NOTE: For the present implementation, the method is more efficient
 * when a lower triangular matrix is used.
 *
 * Dragging along:
 * If Z (r-by-n) is given, so must be the r-vector y. In this case,
 * we overwrite Z by Z_, where
 *   Z_ L_' = Z L' + y v'.
 *
 * Working array:
 * The method uses n Givens rotations, param. by angles c_k, s_k. The
 * change L -> L' is specified by these 2*n numbers. They are written
 * into CVEC, SVEC.
 * Requires a working vector of size >= max(n,r), passed in WORKV.
 * NOTE: The same vector can be passed for VEC and WORKV, in which case
 * VEC is overwritten in an undefined way. If r > n, VEC can be of size
 * r, containing v in the first n components.
 *
 * Input:
 * - L:     Factor L (or L'), overwritten by L_ (or L_'). Must be
 *          lower (upper) triangular, str. code UPLO
 * - VEC:   Vector v. Can have size >n, only first n elem. are used
 * - CVEC:  Vector [n]. c_k ret. here
 * - SVEC:  Vector [n]. s_k ret. here
 * - WORKV: Working vector of size max(n,r). Can be same as VEC
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
 * The method is adapted from LINPACK dchud. We did the following
 * modifications:
 * - Using BLAS drot in order to avoid any explicit O(n^2) loops
 * - Keeping diag(L_) positive, by flipping angles c_k, s_k whenever
 *   a negative element pops up
 * See the TR
 *   M. Seeger
 *   Low Rank Updates for the Cholesky Decomposition
 *   Available at: www.kyb.tuebingen.mpg.de/bs/people/seeger/papers/
 */

/* Main function CHOLUPRK1 */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,n,r=0,stp,sz,ione=1,retcode=0;
  double temp;
  fst_matrix lmat,zmat;
  const double* vvec;
  double* cvec,*svec,*wkvec,*yvec;
  bool islower;
  double* tbuff;

  /* Read arguments */
  if (nrhs<5)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs>1)
    mexErrMsgTxt("Too many return arguments");
  parseBLASMatrix(prhs[0],"L",&lmat,-1,-1);
  /*
  sprintf(errMsg,"%d, %d, '%s', %c\n",lmat.m,lmat.n,lmat.strcode,
	  UPLO(lmat.strcode));
  mexPrintf(errMsg);
  */
  if ((n=lmat.n)!=lmat.m ||
      (!(islower=(UPLO(lmat.strcode)=='L')) && UPLO(lmat.strcode)!='U'))
    mexErrMsgTxt("L must be lower/upper triangular (use UPLO str. code!)");
  if (getVecLen(prhs[1],"VEC")<n) mexErrMsgTxt("VEC too short");
  vvec=mxGetPr(prhs[1]);
  if (getVecLen(prhs[2],"CVEC")!=n || getVecLen(prhs[3],"SVEC")!=n)
    mexErrMsgTxt("CVEC, SVEC have wrong size");
  cvec=mxGetPr(prhs[2]); svec=mxGetPr(prhs[3]);
  if (nrhs>5) {
    if (nrhs<7) mexErrMsgTxt("Need both Z, Y");
    parseBLASMatrix(prhs[5],"Z",&zmat,-1,n);
    r=zmat.m;
    if (getVecLen(prhs[6],"Y")!=r) mexErrMsgTxt("Y has wrong size");
    yvec=mxGetPr(prhs[6]);
  }
  i=getVecLen(prhs[4],"WORKV");
  if (i<n || i<r) mexErrMsgTxt("WORKV too short");
  wkvec=mxGetPr(prhs[4]);

  /* Generate Givens rotations, update L */
  BLASFUNC(dcopy) (&n,vvec,&ione,wkvec,&ione);
  stp=islower?1:lmat.stride;
  for (i=0,sz=n,tbuff=lmat.buff; i<n-1; i++) {
    /* drotg(a,b,c,s): J = [c s; -s c], s.t. J [a; b] = [r; 0]
       a overwritten by r, b by some other information (NOT 0!) */
    if (*tbuff==0.0 && wkvec[i]==0.0) {
      retcode=1; break;
    }
    BLASFUNC(drotg) (tbuff,wkvec+i,cvec+i,svec+i);
    /* Do not want negative elements on factor diagonal */
    if ((temp=*tbuff)<0.0) {
      *tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
    } else if (temp==0.0) {
      retcode=1; break;
    }
    /* drot(x,y,c,s): J = [c s; -s c]. [x_i; y_i] overwritten by
       J [x_i; y_i], for all i
       BAD: Slower for upper triangular! */
    sz--;
    BLASFUNC(drot) (&sz,tbuff+stp,&stp,wkvec+(i+1),&ione,cvec+i,svec+i);
    tbuff+=(lmat.stride+1);
  }
  if (retcode==0 && (*tbuff!=0.0 || wkvec[n-1]!=0.0)) {
    BLASFUNC(drotg) (tbuff,wkvec+(n-1),cvec+i,svec+i);
    if ((temp=*tbuff)<0.0) {
      *tbuff=-temp; cvec[i]=-cvec[i]; svec[i]=-svec[i];
    } else if (temp==0.0) retcode=1;
  } else retcode=1;

  /* Dragging along */
  if (r>0 && retcode==0) {
    BLASFUNC(dcopy) (&r,yvec,&ione,wkvec,&ione);
    for (i=0; i<n; i++)
      BLASFUNC(drot) (&r,zmat.buff+(i*zmat.stride),&ione,wkvec,&ione,cvec+i,
		      svec+i);
  }

  if (nlhs==1) {
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[0]))=(double) retcode;
  }
}
