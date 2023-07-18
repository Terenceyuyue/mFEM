function [fd,fh,BdBox,pfix] = Cube_Domain

x1 = 0; x2 = 1;
y1 = 0; y2 = 1;
z1 = 0; z2 = 1;

fd = @(p) dcube(p,x1,x2,y1,y2,z1,z2);
fh = @huniform;
BdBox = [x1 x2 y1 y2 z1 z2];
pfix = [
    x1, y1, z1;
    x1, y1, z2;
    x1, y2, z1;
    x1, y2, z2;
    x2, y1, z1;
    x2, y1, z2;
    x2, y2, z1;
    x2, y2, z2 ];
