function [fd,fh,BdBox,pfix] = Polygon_Domain

fd = @dpoly;
fh = @huniform;
BdBox = [-1 2 -1 1];
pfix = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
            1.6  0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
        
