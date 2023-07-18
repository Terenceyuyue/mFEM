function [fd,fh,BdBox,pfix,h0] = NACA0012_airfoil

% 机翼的边界函数
% y= +- 0.6*[0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4]
% 机翼是上下对称的，+ 是上部分，- 是下部分

hlead = 0.01; htrail = 0.04; hmax = 2; 
xc = 2; yc = 0;
r = 4;
a = .12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1015];

fd = @DistFnc;
fh = @(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);

fixx = 1-htrail*cumsum(1.3.^(0:4)');
fixy = a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
pfix = [[xc+[-1,1,0,0]*r; 0,0,r*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
BdBox = [xc-r xc+r -r r];
h0 = min([hlead,htrail,hmax]);


    function Dist = DistFnc(p)
        d1 = dcircle(p,xc,yc,r);
        d2 = (abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1);
        Dist = ddiff(d1,d2);
    end


end

