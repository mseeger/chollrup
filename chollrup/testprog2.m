n=100;

% Create matrix A with controlled spectrum
maxlam=2; minlam=0.1;
[q,r]=qr(randn(n,n));
a=muldiag(q,rand(n,1)*(maxlam-minlam)+minlam)*q';
lfact=chol(a)';

% Test rank-2 updates
for i=1:1
  uvec=randn(n,1); vvec=randn(n,1);
  if i==1
    uvec=zeros(n,1); uvec(20)=1; uvec(100)=-1; vvec=-uvec;
  end
  fprintf(1,'vtu=%f\n',uvec'*vvec);
  b=randn(3*n,n);
  x1=b/(lfact'); x2=b*lfact;
  l=lfact; l(1,1)=l(1,1)+1; l(1,1)=l(1,1)-1;
  if choluprk2(l,uvec,vvec,0,0,x1,x2)~=0
    error('Numerical error in CHOLUPRK2!');
  end
  tvec=a*uvec; tscal=uvec'*tvec;
  l_2=chol(a+tvec*vvec'+vvec*tvec'+tscal*vvec*vvec')';
  fprintf(1,'Max. dist. L: %f\n',max(max(abs(l-l_2))));
  x1_2=b/(l_2'); x2_2=b*l_2;
  fprintf(1,'Max. dist. X1: %f\nMax. dist. X2: %f\n', ...
	  max(max(abs(x1-x1_2))),max(max(abs(x2-x2_2))));
end
