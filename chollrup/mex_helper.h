/* -------------------------------------------------------------------
 * Helper functions for MEX functions
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef MEX_HELPER_H
#define MEX_HELPER_H

#include <math.h>
#include "mex.h"

char errMsg[];

/*
 * The BLAS/LAPACK function xxx is called as xxx_ in Linux, but as
 * xxx in Windows. Uncomment the corresponding definition here.
 */
/* Version for Linux */
#define BLASFUNC(NAME) NAME ## _
/* Version for Windows */
/* #define BLASFUNC(NAME) NAME */

/*
 * Types for matrix arguments to support BLAS convention.
 * See 'parseBLASMatrix'.
 */
typedef struct {
  double* buff;
  int m,n;
  int stride;
  char strcode[4];
} fst_matrix;

/*
 * Macros picking out specific positions from structure code array
 */
#define UPLO(arr) (arr)[0]

#define DIAG(arr) (arr)[2]

/*
 * Exported functions
 */

double getScalar(const mxArray* arg,const char* name)
{
  if (!mxIsDouble(arg) || mxGetM(arg)!=1 || mxGetN(arg)!=1) {
    sprintf(errMsg,"Expect double scalar for %s",name);
    mexErrMsgTxt(errMsg);
  }
  return *mxGetPr(arg);
}

int getScalInt(const mxArray* arg,const char* name)
{
  double val,temp;

  if (!mxIsDouble(arg) || mxGetM(arg)!=1 || mxGetN(arg)!=1) {
    sprintf(errMsg,"Expect scalar for %s",name);
    mexErrMsgTxt(errMsg);
  }
  val=*mxGetPr(arg);
  if ((temp=floor(val))!=val) {
    sprintf(errMsg,"Expect integer for %s",name);
    mexErrMsgTxt(errMsg);
  }
  return (int) temp;
}

int getVecLen(const mxArray* arg,const char* name)
{
  int n;

  if (!mxIsDouble(arg)) {
    sprintf(errMsg,"Expect real vector for %s",name);
    mexErrMsgTxt(errMsg);
  }
  if (mxGetM(arg)==0 || mxGetN(arg)==0)
    return 0;
  if ((n=mxGetM(arg))==1)
    n=mxGetN(arg);
  else if (mxGetN(arg)!=1) {
    sprintf(errMsg,"Expect real vector for %s",name);
    mexErrMsgTxt(errMsg);
  }
  return n;
}

/*
 * NOTE: The string returned is allocated here using 'mxMalloc',
 * it has to be dealloc. by the user using 'mxFree'.
 */
const char* getString(const mxArray* arg,const char* name)
{
  int len;
  char* buff;

  if (!mxIsChar(arg) || mxGetM(arg)!=1) {
    sprintf(errMsg,"Expect char row vector for %s",name);
    mexErrMsgTxt(errMsg);
  }
  len=mxGetN(arg)+1;
  buff=(char*) mxMalloc(len*sizeof(char));
  mxGetString(arg,buff,len);

  return buff;
}

void checkMatrix(const mxArray* arg,const char* name,int m,int n)
{
  if (!mxIsDouble(arg)) {
    sprintf(errMsg,"Expect real matrix for %s",name);
    mexErrMsgTxt(errMsg);
  }
  if (m!=-1 && mxGetM(arg)!=m) {
    sprintf(errMsg,"Expect %d rows for %s",m,name);
    mexErrMsgTxt(errMsg);
  }
  if (n!=-1 && mxGetN(arg)!=n) {
    sprintf(errMsg,"Expect %d columns for %s",n,name);
    mexErrMsgTxt(errMsg);
  }
}

/*
 * Matrix argument can come with additional BLAS attributes. To pass these,
 * 'arg' must be a cell vector { BUFF, [YS XS M N], {SCODE} }. Here, BUFF
 * is a normal (buffer) matrix, YS, XS, M, N are integers, SCODE is a
 * string (optional). The matrix is BUFF(YS:(YS+M-1),XS:(XS+N-1)).
 *
 * Structure codes:
 * If SCODE is given, it must be a string of length 2. Pos.:
 * - 0: UPLO field, values 'U' (upper), 'L' (lower)
 * - 1: DIAG field, values 'U' (unit tri.), 'N' (non unit tri.)
 *      If UPLO spec., the def. value for DIAG is 'N'
 * A field value ' ' means: not specified.
 * These are passed to BLAS routines if required. The 'strcode' field
 * contains the codes separ. by 0, i.e. the C string for the codes
 * attached to each other.
 */
void parseBLASMatrix(const mxArray* arg,const char* name,
  fst_matrix* mat,int m,int n)
{
  int bm,bn,ys,xs,am,an,csz;
  const mxArray* bmat,*szvec,*scdvec;
  const double* iP;
  char sbuff[3];

  mat->strcode[0]=mat->strcode[2]=' ';
  mat->strcode[1]=mat->strcode[3]=0;
  if (!mxIsCell(arg)) {
    /* No cell array: Must be normal matrix */
    checkMatrix(arg,name,m,n);
    mat->buff=mxGetPr(arg);
    mat->stride=mat->m=mxGetM(arg); mat->n=mxGetN(arg);
  } else {
    if ((csz=mxGetM(arg)*mxGetN(arg))<2) {
      sprintf(errMsg,"Array %s has wrong size",name);
      mexErrMsgTxt(errMsg);
    }
    bmat=mxGetCell(arg,0);
    checkMatrix(bmat,name,-1,-1);
    bm=mxGetM(bmat); bn=mxGetN(bmat);
    if (getVecLen(mxGetCell(arg,1),name)!=4) {
      sprintf(errMsg,"Index vector in %s has wrong size",name);
      mexErrMsgTxt(errMsg);
    }
    iP=mxGetPr(mxGetCell(arg,1));
    ys=((int) iP[0])-1; xs=((int) iP[1])-1;
    am=(int) iP[2]; an=(int) iP[3];
    if (ys<0 || xs<0 || am<0 || an<0 || ys+am>bm || xs+an>bn) {
      sprintf(errMsg,"Index vector in %s wrong",name);
      mexErrMsgTxt(errMsg);
    }
    if ((m!=-1 && am!=m) || (n!=-1 && an!=n)) {
      sprintf(errMsg,"Matrix %s has wrong size",name);
      mexErrMsgTxt(errMsg);
    }
    mat->buff=mxGetPr(bmat)+(xs*bm+ys);
    mat->m=am; mat->n=an; mat->stride=bm;
    if (csz>2) {
      /* Structure codes */
      scdvec=mxGetCell(arg,2);
      if (!mxIsChar(scdvec) || mxGetM(scdvec)!=1 ||
	  mxGetN(scdvec)!=2) {
	sprintf(errMsg,"Structure code string in %s wrong",name);
	mexErrMsgTxt(errMsg);
      }
      mxGetString(scdvec,sbuff,3);
      if ((sbuff[0]!='U' && sbuff[0]!='L' && sbuff[0]!=' ') ||
	  (sbuff[1]!='U' && sbuff[1]!='N' && sbuff[1]!=' ')) {
	sprintf(errMsg,"Structure code string in %s wrong",name);
	mexErrMsgTxt(errMsg);
      }
      if (sbuff[0]!=' ' && sbuff[1]==' ')
	sbuff[1]='N'; /* def. value */
      if ((mat->m!=mat->n && sbuff[0]!=' ') ||
	  (sbuff[1]!=' ' && sbuff[0]==' ')) {
	sprintf(errMsg,"Structure code string in %s inconsistent",name);
	mexErrMsgTxt(errMsg);
      }
      mat->strcode[0]=sbuff[0]; mat->strcode[2]=sbuff[1];
    }
  }
}

void fillVec(double* vec,int n,double val)
{
  int i;
  for (i=0; i<n; i++) vec[i]=val;
}

bool isUndef(double val)
{
  return isnan(val) || (isinf(val)!=0);
}

#endif
