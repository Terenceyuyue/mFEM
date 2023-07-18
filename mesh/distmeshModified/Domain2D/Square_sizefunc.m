function [fd,fh,BdBox,pfix] = Square_sizefunc

fd = @(p) drectangle(p,0,1,0,1);
fh = @(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
    0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
BdBox = [0 1 0 1];
pfix = [0,0;1,0;0,1;1,1];
