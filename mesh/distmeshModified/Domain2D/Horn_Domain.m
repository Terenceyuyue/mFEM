
function [fd,fh,BdBox,pfix] = Horn_Domain

fd = @fd1;
fh = @fh1;
BdBox = [-1 1 0 1];
pfix = [-1,0;-.95,0;.1,0;1,0];

    function d = fd1(p)        
        d1 = dcircle(p,0,0,1);
        d2 = dcircle(p,-0.4,0,0.55);
        d3 = dline(p,0,0,1,0);
        d = dintersect(d3,ddiff(d1,d2));
    end

    function d = fh1(p)
        d1 = dcircle(p,0,0,1);  d2 =  dcircle(p,-.4,0,.55);        
        h1 = (0.15-0.2*d1); h2 = (0.06+0.2*d2); h3 = (d2-d1)/3;        
        d = min(min(h1,h2),h3);
    end

      
        
end