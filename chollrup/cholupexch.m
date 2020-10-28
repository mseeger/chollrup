%CHOLUPEXCH Update Cholesky factor for special permutation
%  CHOLUPEXCH(R,K,L,JOB,{X=[]})
%
%  ATTENTION: We use the undocumented fact that the content of
%  matrices passed as arguments to a MEX function can be overwritten
%  like in a proper call-by-reference. This is not officially
%  supported and may not work in future Matlab versions!
%
%  Update of Cholesky factor after special row/column permutation
%  of system matrix.
%    A = R' R, A_ = E' A E, E spec. permut. matrix
%  Here, R is upper triangular. The factor must be passed in R and is
%  overwritten by R_. The method works by computing an orthonormal
%  U s.t. U R E = R'.
%  E is determined by K, L, JOB (integers). Here, 1 <= K < L <= n
%  (size of A). E reorders rows/columns as follows:
%  If JOB==1, the new ordering is ...,K-1,L,K,...,L-1,L+1,... In this
%  case, U = U(L-K)*...*U(1), where U(i) is a Givens rotation in the
%  plane (L-i,L-i+1). The Givens form is [c(i) s(i); -s(i) c(i)].
%  If JOB==2, the new ordering is ...,K-1,K+1,...,L,K,L+1,... In this
%  case, U = U(L-K)*...*U(1), where U(i) is a Givens rotation in the
%  plane (K+i-1,K+i).
%
%  Dragging along: If X is passed, it is replaced by X_ = U X. Note
%  that if R' X = B, then (R_)' X_ = E' B
%
%  The method is just a wrapper around the LINPACK routine DCHEX
%  NOTE: There is a bug in DCHEX, causing elements of DIAG(R_) to be
%  negative. We include a workaround here.
