% Distmesh2d examples.
% Copyright (C) Terence YUE Yu.

clc;clear;close all

Domain = {
    % Domain                        % h
    @Lshape_Domain,                 0.20; 
    @Circle_Domain,                 0.20;
    @Rectangle_Domain,              0.20;       
    @Ellipse_Domain,                0.20;
    @Polygon_Domain,                0.10;
    @Rectangle_Circle_Domain,       0.05;
    @Circle_Circle_Domain,          0.05;
    @Horn_Domain,                   0.01;
    @Suspension_Domain,             1.00;   
    @Upper_Circle_Domain,           0.20;
    @Upper_Circle_Circle_Domain,    0.10;
    @Wrench_Domain,                 0.05;
    @Square_sizefunc,               0.01;
    @NACA0012_airfoil,              0.01;
    @Pie_hole,                      0.005;
    @Superellipse,                  0.08;
    @Implicit,                      0.75;
    };
% very slow: Superellipse, Implicit

for id = 1  % very slow: Superellipse, Implicit (id = 16,17)
    
    func = Domain{id,1};
    str = func2str(func);
    [fd,fh,BdBox,pfix] = feval(func);
    
    fprintf('\n ------ %s ------ \n', str);
    
    h0 = Domain{id,2};
    
    if mycontains(str, 'Polygon')
        [node,elem] = distmesh2d(fd,fh,h0,BdBox,pfix,pfix);
    else
        [node,elem] = distmesh2d(fd,fh,h0,BdBox,pfix);
    end
    
    clf, showmesh(node,elem)
    str1 = strrep(str,'_','\_');
    title(sprintf('%s',str1))
    pause(1)
    
    
end


