% Ref: 上半圆环区域
% dline: 左右点连线 (线的左侧区域是区域内部，左侧按逆时针理解)

function [fd,fh,BdBox,pfix] = Upper_Circle_Circle_Domain
BdBox = [-1 1 0 1];
fd = @DistFnc;
fh = @huniform;
pfix = [-1 0; -0.5 0; 0.5 0; 1 0];

    function Dist = DistFnc(p)
        d1 = dcircle(p,0,0,1);
        d2 = dcircle(p,0,0,0.5);
        d3 = dline(p,0,0,1,0);
        Dist = dintersect(d3,ddiff(d1,d2));
    end

end