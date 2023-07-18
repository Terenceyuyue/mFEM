function [fd,fh,BdBox,pfix] = Cylinder_Hole_Domain


fd = @DistFunc;
fh = @hFunc;
BdBox = [-1 1 -1 1 -1 1];
pfix = [];


    function d = DistFunc(p)
        r = sqrt ( p(:,1).^2 + p(:,2).^2 );
        z = p(:,3);
        
        d1 = r - 1.0;
        d2 = z - 1.0;
        d3 = - z - 1.0;
        d4 = sqrt(d1.^2 + d2.^2);
        d5 = sqrt(d1.^2 + d3.^2);
        
        d = dintersect( dintersect ( d1, d2 ), d3 );
        ix = ( 0.0 < d1 ) & ( 0.0 < d2 );
        d(ix) = d4(ix);
        ix = ( 0.0 < d1 ) & ( 0.0 < d3 );
        d(ix) = d5(ix);
        
        d = ddiff( d, dsphere ( p, 0.0, 0.0, 0.0, 0.5 ) );
    end

    function h = hFunc(p)
        h1 = 4.0*sqrt( sum( p.^2, 2 ) ) - 1.0;        
        h = min ( h1, 2.0 );
    end

end