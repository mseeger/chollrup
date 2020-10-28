n=1000;
maxlam=2; minlam=0.1;
cvec=zeros(n,1); svec=zeros(n,1);
wkvec=zeros(3*n,1);

% Test positive updates
fprintf(1,'Testing pos. updates\n');
for i=1:10
  % Create matrix A with controlled spectrum
  [q,r]=qr(randn(n,n));
  a=muldiag(q,rand(n,1)*(maxlam-minlam)+minlam)*q';
  lfact=chol(a)';
  % Update
  b=randn(3*n,n); y=randn(3*n,1);
  z=b/lfact';
  vec=randn(n,1);
  if choluprk1({lfact,[1 1 n n],'L '},vec,cvec,svec,wkvec,z,y)~=0
    error('Numerical error in CHOLUPRK1!');
  end
  l_2=chol(a+vec*vec')';
  fprintf(1,'Max. dist. L: %e\n',max(max(abs(lfact-l_2))));
  z_2=(b+y*vec')/l_2';
  fprintf(1,'Max. dist. Z: %e\n',max(max(abs(z-z_2))));
end

% Test negative updates
fprintf(1,'Testing neg. updates\n');
for i=1:10
  % Create matrix A with controlled spectrum
  [q,r]=qr(randn(n,n));
  a=muldiag(q,rand(n,1)*(maxlam-minlam)+minlam)*q';
  vec=randn(n,1);
  lfact=chol(a+vec*vec')';
  % Downdate
  b=randn(3*n,n); y=randn(3*n,1);
  z=(b+y*vec')/lfact';
  isp=double(rand>=0.5);
  if isp
    vec=lfact\vec;
  end
  if choldnrk1({lfact,[1 1 n n],'L '},vec,cvec,svec,wkvec, ...
	       isp,z,y)~=0
    error('Numerical error in CHOLDNRK1!');
  end
  l_2=chol(a)';
  fprintf(1,'Max. dist. L: %e\n',max(max(abs(lfact-l_2))));
  z_2=b/l_2';
  fprintf(1,'Max. dist. Z: %e\n',max(max(abs(z-z_2))));
end
