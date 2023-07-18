
clc;clear;close all

Domain = {@Sphere_Domain, @Cube_Domain, @Cylinder_Domain, @Cylinder_Hole_Domain};
h = [0.2 0.2 0.2 0.1];

for id = 1:4
    
    func = Domain{id};
    str = func2str(func);
    [fd,fh,BdBox,pfix] = feval(func);
    
    fprintf('\n ------ %s ------ \n', str);
    
    h0 = h(id);
    
    [node,elem] = distmesh3d(fd,fh,h0,BdBox,pfix); 
    
    clf
    options.FaceAlpha = 1; 
    if mycontains(str, 'Hole')
        options.expr = 'y>0';
    end
    showmesh(node,elem,options);
    str1 = strrep(str,'_','\_');
    title(sprintf('%s',str1))
    pause(0.1)
    
    
end

