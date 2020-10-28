%CHOLUPRK1 Update Cholesky factor for rank-one modification
%  STAT=CHOLUPRK1(L,VEC,CVEC,SVEC,WORKV,{Z,Y})
%
%  ATTENTION: We use the undocumented fact that the content of
%  matrices passed as arguments to a MEX function can be overwritten
%  like in a proper call-by-reference. This is not officially
%  supported and may not work in future Matlab versions!
%
%  If
%    A = L*L', A_ = A + v*v' = L_ L_', A, L n-by-n,
%  where L is lower triangular, the method computes L_ from L. This
%  is called Cholesky rank one update. L or L' (upper triangular) can
%  be passed, only the relevant triagle is accessed. L (or L') is
%  passed in L, v in VEC.
%  NOTE: For the present implementation, the method is more efficient
%  when a lower triangular matrix is used.
%
%  Dragging along:
%  If Z (r-by-n) is given, so must be the r-vector y. In this case,
%  we overwrite Z by Z_, where
%    Z_ L_' = Z L' + y v'.
%
%  Working array:
%  The method uses n Givens rotations, param. by angles c_k, s_k. The
%  change L -> L' is specified by these 2*n numbers. They are written
%  into CVEC, SVEC.
%  Requires a working vector of size >= max(n,r), passed in WORKV.
%  NOTE: The same vector can be passed for VEC and WORKV, in which case
%  VEC is overwritten in an undefined way. If r > n, VEC can be of size
%  r, containing v in the first n components.
