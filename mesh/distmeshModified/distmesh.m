function [p,t] = distmesh(fd,h0,BdBox,pfix)

% A simplified version of distmesh2d.m

% Copyright (C) 2004-2012 Per-Olof Persson.
% @ Modified by Terence YUE Yu

if nargin==3, pfix = []; end

geps = 0.001*h0; Fscale = 1.2; dt = 0.2;
deps = sqrt(eps)*h0;  dptol = 1e-3;

% 1. Create initial distribution in bounding box (equilateral triangles)
x = BdBox(1):h0:BdBox(2); y = BdBox(3):h0*sqrt(3)/2:BdBox(4);
[x,y] = meshgrid(x,y);
x(2:2:end,:) = x(2:2:end,:) + h0/2;                  % Shift even rows
p = [x(:),y(:)];                                     % List of node coordinates

% 2. Remove points outside the region
p = p(fd(p)<geps,:);                 % Keep only d<0 points
if ~isempty(pfix), p = setdiff(p,pfix,'rows'); end  % Remove duplicated nodes
p = [pfix; p];
N = size(p,1);

iter = 0;  MaxIter = 1e3;
while iter < MaxIter
    iter = iter+1;
    % 3. Get the truss and bars by the Delaunay algorithm
    t = delaunay(p);                                  
    pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;  % Compute centroids
    t = t(fd(pc)<-geps,:);                         % Keep interior triangles
    totalEdge = sort([t(:,[2 3]); t(:,[3 1]); t(:,[1 2])],2);
    bars = unique(totalEdge,'rows');               % Bars as node pairs
        
    % 4. Update the points
    barvec = p(bars(:,2),:)-p(bars(:,1),:);          % List of bar vectors
    L = sqrt(sum(barvec.^2,2));                      % Bar lengths
    L0 = Fscale*sqrt(sum(L.^2)/length(L));           % Desired lengths    
    f = max(L0-L,0);                                 % Bar forces (scalars)
    F = f./L*[1,1].*barvec;                          % Bar forces (x,y components)    
    % Force resultant
    Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(L))*[1,2,1,2],[-F,F],N,2)); 
    Ftot(1:size(pfix,1),:) = 0;                      % Force = 0 at fixed points
    p = p + dt*Ftot;
    
    % 5. Bring outside points back to the boundary
    d = fd(p); ix = d>0;                 % Find points outside (d>0)
    n1 = (fd([p(ix,1)+deps,p(ix,2)])-d(ix))/deps; 
    n2 = (fd([p(ix,1),p(ix,2)+deps])-d(ix))/deps; 
    nn = n1.^2 + n2.^2;
    n1 = n1./nn; n2 = n2./nn;
    p(ix,:) = p(ix,:)-[d(ix).*n1,d(ix).*n2]; % Project back to the boundary
    
    % 6. Termination criterion: All interior nodes move less than dptol (scaled)
    if max(sqrt(sum(dt*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
end

% final mesh
t = delaunay(p);  
pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;      
t = t(fd(pc)<-geps,:);   