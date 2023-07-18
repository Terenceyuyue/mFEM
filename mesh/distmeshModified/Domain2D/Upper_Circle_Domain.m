
% Ref: 上半圆区域
% dline: 左右点连线 (线的左侧区域是区域内部，左侧按逆时针理解)

function [fd,fh,BdBox,pfix] = Upper_Circle_Domain
BdBox = [-1 1 0 1];
fd = @DistFnc;
fh = @huniform;
pfix = [-1 0; 1 0];


    function Dist = DistFnc(p)
        d1 = dcircle(p,0,0,1);
        d2 = dline(p,0,0,1,0);
        Dist = dintersect(d2,d1);
    end

end