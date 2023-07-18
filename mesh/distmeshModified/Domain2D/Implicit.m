function [fd,fh,BdBox,pfix] = Implicit

BdBox = [-5*pi/2 5*pi/2 -5 2];
fd = @DistFunc;
fh = @huniform;
pfix = 5*pi/2*[-1,0; 1,0];

    function Dist = DistFunc(p)
        d1=dexpr(p,'y-cos(x)');
        d2=dexpr(p,'-(y-(-5+5/(5/4*2*pi)^4*x^4))');
        Dist=dintersect(d1,d2);
    end

end