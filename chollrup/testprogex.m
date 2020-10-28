n=100;

% Create matrix A with controlled spectrum
maxlam=2; minlam=0.1;
[q,r]=qr(randn(n,n));
a=muldiag(q,rand(n,1)*(maxlam-minlam)+minlam)*q';
rfact=chol(a);

% Test exchange updates
for i=1:20
  k=floor(rand*(n-1))+1;
  l=k+1+floor(rand*(n-k));
  b=randn(n,3*n);
  x=(rfact')\b;
  r=rfact; r(1,1)=r(1,1)+1; r(1,1)=r(1,1)-1;
  cholupexch(r,k,l,2,x);
  ind=[1:(k-1) (k+1):l k (l+1):n];
  apr=a(ind,ind);
  r_2=chol(apr);
  x_2=(r_2')\b(ind,:);
  fprintf(1,'k=%d, l=%d\n',k,l);
  fprintf(1,'Max. dist. R: %f\n',max(max(abs(r-r_2))));
  fprintf(1,'Max. dist. X: %f\n',max(max(abs(x-x_2))));
  if max(max(abs(r-r_2)))>(1e-6)
    imagesc(abs(r-r_2));
    ind=find(abs(r-r_2)>(1e-6))
    r(ind)
    r_2(ind)
    diag(r)
    pause;
  end
end
