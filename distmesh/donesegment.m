function ds=donesegment(p,pv)
% 求点p=[x,y]到给定线段pv=[x1,y1;x2,y2]的距离,p可以是向量。
e=ones(size(p,1),1);

v=diff(pv,1);
w=p-e*pv(1,:);

c1=sum(w.*v(e,:),2);
c2=sum(v(e,:).^2,2);

ds=0*e;

ix=c1<=0;
ds(ix)=sqrt(sum((p(ix,:)-pv(1*ones(sum(ix),1),:)).^2,2));

ix=c1>=c2;
ds(ix)=sqrt(sum((p(ix,:)-pv(2*ones(sum(ix),1),:)).^2,2));

ix=c1>0 & c2>c1;
nix=sum(ix);
if nix>0
  Pb=ones(nix,1)*pv(1,:)+c1(ix)./c2(ix)*v;
  ds(ix)=sqrt(sum((p(ix,:)-Pb).^2,2));
end
