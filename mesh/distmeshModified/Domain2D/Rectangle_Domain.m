function [fd,fh,BdBox,pfix] = Rectangle_Domain

a = -1; b = 1; c = -1; d = 1;

fd = @(p) drectangle(p,a,b,c,d);
fh = @huniform;
BdBox = [a b c d];
pfix = [a c; b c; b d; a d];
