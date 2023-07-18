
function [fd,fh,BdBox,pfix] = Suspension_Domain
BdBox = [-2 24 -2 24];
fd = @fd1;
fh = @huniform;
pfix = [ 2 2; 2 16; 20 2.5];

    function Dist = fd1(p)
        d1  = drectangle(p,0,18.885,0,14.56);
        d2  = dline(p,18.885,1.3030,4,0);
        d3  = dline(p,3.92,14.56,6.1699,6.88);
        d4  = dline(p,9.8651,4.0023,18.885,3.70);
        d5  = dline(p,4,0,0,4);
        d13 = dline(p,0,14,3.92,14.56);
        d14 = dcircle(p,10,8,4);
        d15 = dline(p,9.8651,4.0023,6.1699,6.88);
        d   = ddiff(ddiff(ddiff(ddiff(d1,d2),d5),d13),...
            dunion(ddiff(dintersect(d3,d4),d15),d14));
        d6  = dcircle(p,2,2,2);
        d7  = dcircle(p,4,2,2);
        d8  = dcircle(p,2,4,2);
        d   = dunion(d,dunion(d6,dunion(d7,d8)));
        d9  = dcircle(p,2,14,2);
        d10 = dcircle(p,2,16,2);
        d11 = dcircle(p,18.885,2.5,1.2);
        d12 = dcircle(p,20,2.5,1.2);
        Dist = dunion(d,dunion(d9,dunion(d10,dunion(d11,d12))));
    end



end