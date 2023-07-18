function [fd,fh,BdBox,pfix] = Cylinder_Domain

xc = 0; yc = 0; r = 1;
z1 = 0; z2 = 4;

fd = @(p) dcylinder(p,xc,yc,r,z1,z2);
fh = @huniform;
BdBox = [xc-r xc+r yc-r yc+r z1 z2];
pfix = [];
