function [fd,fh,BdBox,pfix] = Circle_Domain

xc = 0; yc = 0; r = 1;

fd = @(p) dcircle(p,xc,yc,r);
fh = @huniform; 
BdBox = [-1 1 -1 1];
pfix = [];