function [fd,fh,BdBox,pfix] = Superellipse

BdBox = [-1.05 1 -1.05 1];
fd = @DistFunc;
fh = @hFunc;
pfix = [];

    function Dist = DistFunc(p)
        d1 = dexpr(p,'(x^4+y^4)^(1/4)-1');
        d2 = dexpr(p,'(x^4+y^4)^(1/4)-0.5');
        Dist = ddiff(d1,d2);
    end

    function h = hFunc(p)
        h = dexpr(p,'(x^4+y^4)^(1/4)');
    end

end