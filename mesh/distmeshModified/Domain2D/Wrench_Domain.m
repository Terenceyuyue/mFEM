
function [fd,fh,BdBox,pfix] = Wrench_Domain
BdBox = [-0.3 2.5 -0.5 0.5];
fd = @DistFnc;
fh = @huniform;
pfix = [];

    function Dist = DistFnc(p)
        d1 = dline(p,0,0.3,0,-0.3);
        d2 = dline(p,0,-0.3,2,-0.5);
        d3 = dline(p,2,-0.5,2,0.5);
        d4 = dline(p,2,0.5,0,0.3);
        d5 = dcircle(p,0,0,0.3);
        d6 = dcircle(p,2,0,0.5);
        douter = dunion(d6,dunion(d5,...
            dintersect(d4,dintersect(d3,dintersect(d2,d1)))));
        d7 = dcircle(p,0,0,0.175);
        d8 = dcircle(p,2,0,0.3);
        din = dunion(d8,d7);
        Dist = ddiff(douter,din);
    end

end