function [fd,fh,BdBox,pfix] = Circle_Circle_Domain

fd = @fd1;
fh = @(p) min(4*sqrt(sum(p.^2,2))-1,2);
BdBox = [-1 1 -1 1];
pfix = [];

    function d = fd1(p)
        d1 = dcircle(p,0,0,1);
        d2 = dcircle(p,0,0,0.4);
        d = ddiff(d1,d2);
    end
        

end