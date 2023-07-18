function d = dexpr(p,fin,nit,alpha)

% Copyright (C) 2004-2012 Per-Olof Persson. 
% @ Modified by Terence YUE Yu

if nargin<3, nit = 20; end
if nargin<4, alpha = 0.1; end

syms x y;

fin = eval(fin);

fx = matlabFunction(diff(fin,x),'Vars',{x,y});    
fy = matlabFunction(diff(fin,y),'Vars',{x,y});   
fxx = matlabFunction(diff(fin,x,2),'Vars',{x,y}); 
fyy = matlabFunction(diff(fin,y,2),'Vars',{x,y});  
fxy = matlabFunction(diff(diff(fin,x),y),'Vars',{x,y});  
f = matlabFunction(fin,'Vars',{x,y});   

x0=p(:,1);
y0=p(:,2);
x=x0;
y=y0;
for it=1:nit
  cf=f(x,y);
  cfx=fx(x,y);
  cfy=fy(x,y);
  cfxx=fxx(x,y);
  cfxy=fxy(x,y);
  cfyy=fyy(x,y);

  F1=cf;
  F2=(x-x0).*cfy-(y-y0).*cfx;
  J11=cfx;
  J12=cfy;
  J21=cfy+(x-x0).*cfxy-(y-y0).*cfxx;
  J22=-cfx-(y-y0).*cfxy+(x-x0).*cfyy;
  
  detJ=J11.*J22-J12.*J21;
  detJ(detJ==0)=inf;
  
  x=x-alpha*(J22.*F1-J21.*F2)./detJ;
  y=y-alpha*(-J12.*F1+J11.*F2)./detJ;
end

d=sqrt((x-x0).^2+(y-y0).^2).*sign(feval(f,x0,y0));
