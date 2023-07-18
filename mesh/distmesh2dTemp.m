function [p,t] = distmesh2dTemp(pseeds,fd,fh,h0,BdBox,pfix,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%
%      p:         Node positions (N x 2)
%      t:         Triangle indices (NT x 3)
%      fd:        Distance function d(x,y)
%      fh:        Scaled mesh size function h(x,y)
%      h0:        Initial edge length
%      BdBox:     Bounding box [xmin xmax ymin ymax]
%      pfix:      Fixed node positions (nfix x 2)
%      parameters (varargin):   Additional parameters passed to fd and fh
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%      fd = @(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
%      fh = @(p) 0.05+0.3*dcircle(p,0,0,0.5);
%      h0 = 0.1;
%      BdBox = [-1 1 -1 1];
%      pfix = [-1,-1; -1,1; 1,-1; 1,1];
%      [p,t] = distmesh2d(fd,fh,h0,BdBox,pfix);

%   Copyright (C) 2004-2012 Per-Olof Persson, modified by Terence Yu.

if nargin==4, pfix = []; end

geps = 0.001*h0; Fscale = 1.2; dt = 0.2;
deps = sqrt(eps)*h0;  dptol = 1e-3;

p = pseeds;
if ~isempty(pfix), p = setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
pfix = unique(pfix,'rows');
p = [pfix; p];
N = size(p,1);

iter = 0;  MaxIter = 1e3;
pold = inf; ttol = 0.1;
tic;
while iter < MaxIter
    iter = iter+1;
    % 3. Get the truss and bars by the Delaunay algorithm (only for large
    % movement)
    if max(sqrt(sum((p-pold).^2,2))/h0)>ttol 
        pold = p;
        t = delaunay(p);
        pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;      % Compute centroids
        t = t(feval(fd,pc,varargin{:})<-geps,:);           % Keep interior triangles
        totalEdge = sort([t(:,[2,3]); t(:,[3,1]); t(:,[1,2])],2);
        bars = unique(totalEdge,'rows');                   % Bars as node pairs
    end
    
    % 4. Update the points
    barvec = p(bars(:,2),:)-p(bars(:,1),:);              % List of bar vectors
    L = sqrt(sum(barvec.^2,2));                          % Bar lengths
    pm = (p(bars(:,1),:)+p(bars(:,2),:))/2;
    hbars = feval(fh,pm,varargin{:});
    L0 = Fscale*hbars*sqrt(sum(L.^2)/sum(hbars.^2));   % Desired lengths    
    f = max(L0-L,0);                                   % Bar forces (scalars)
    F = barvec.*repmat(f./L,1,2);                      % Bar forces (x,y components)
    % Force resultant
    Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(L))*[1,2,1,2],[-F,F],N,2)); 
    Ftot(1:size(pfix,1),:) = 0;                      % Force = 0 at fixed points
    p = p + dt*Ftot;
    
    % 5. Bring outside points back to the boundary
    d = feval(fd,p,varargin{:}); ix=d>0;             % Find points outside (d>0)
    n1 = (feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; 
    n2 = (feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; 
    nn = (n1.^2 + n2.^2); % better performance without sqrt
    nn(nn<deps) = 1;
    n1 = n1./nn; n2 = n2./nn; 
    p(ix,:) = p(ix,:)-[d(ix).*n1,d(ix).*n2];    % Project back to the boundary
    
    % 6. Termination criterion: All interior nodes move less than dptol (scaled)
    if max(sqrt(sum(dt*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
    if toc>60*2, break; end  
end

% final mesh
t = delaunay(p);  
pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;      
t = t(feval(fd,pc,varargin{:})<-geps,:);  
