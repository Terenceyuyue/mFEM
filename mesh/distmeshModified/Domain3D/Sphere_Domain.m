function [fd,fh,BdBox,pfix] = Sphere_Domain

xc = 0; yc = 0; zc = 0; r = 1;

fd = @(p) dsphere(p,xc,yc,zc,r);
fh = @huniform;
BdBox = [-1 1 -1 1 -1 1];
pfix = [];
