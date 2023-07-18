function [fd,fh,BdBox,pfix] = Pie_hole

c=(sqrt(119)-9)/20;
pfix=[.9+c,c;cos(pi/12),sin(pi/12)];
pfix=[pfix;pfix];
pfix(3:4,2)=-pfix(3:4,2);
pfix=[pfix;0,0;0.9,0];
BdBox = [0 1 -1 1];

fd = @DistFunc;
fh = @hFunc;


    function d=DistFunc(p)
        d1=dcircle(p,0,0,1);
        d2=drectangle(protate(p,pi/12),-1,1,0,1);
        d3=drectangle(protate(p,-pi/12),-1,1,-1,0);
        d4=drectangle(protate(pshift(p,-.9,0),-pi/4),0,.2,0,.2);
        d5=dcircle(p,.6,0,.1);
        d=ddiff(ddiff(ddiff(ddiff(d1,d2),d3),d4),d5);
    end

    function h=hFunc(p)        
        h1=0.005+0.2*sqrt(sum(p.^2,2));
        h2=0.02+0.2*(sqrt((p(:,1)-.6).^2+p(:,2).^2)-.1);
        h3=0.005+0.2*sqrt((p(:,1)-.9).^2+p(:,2).^2);
        h=min(min(min(h1,h2),h3),0.03);
    end

end