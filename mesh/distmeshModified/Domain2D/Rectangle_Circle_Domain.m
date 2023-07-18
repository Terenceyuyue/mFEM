function [fd,fh,BdBox,pfix] = Rectangle_Circle_Domain
% Square with hole

fd = @fd1;
fh = @(p) 0.05+0.3*dcircle(p,0,0,0.5); 
BdBox = [-1 1 -1 1];
pfix = [-1,-1;-1,1;1,-1;1,1];


    function d = fd1(p)
        d1 = drectangle(p,-1,1,-1,1);
        d2 = dcircle(p,0,0,0.5);
        d = ddiff(d1,d2);
    end

end