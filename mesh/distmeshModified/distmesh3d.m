function [p,t] = distmesh3d(fd, fh, h0, BdBox, pfix, varargin)
% DISTMESH3D 3-D Mesh Generator using Distance Functions.
%
%  Copyright (C) 2004 Per-Olof Persson, modified by Terence Yu.


dim = 3;
geps = 0.1*h0; dptol = 0.001;  dt = 0.1;
deps = sqrt (eps)*h0;
Fscale = 1 + 0.4 / 2^(dim-1);

% 1. Create the initial distribution in the bounding box.
cbox = cell(1,dim);
for i = 1 : dim
    cbox{i} = BdBox(2*i-1) : h0 : BdBox(2*i);
end
pp = cell(1,dim);
[pp{:}] = ndgrid(cbox{:});
p = zeros(numel(pp{1}),dim);
for i = 1 : dim
    p(:,i) = pp{i}(:);
end

% 2. Remove points outside the region, apply the rejection method.
p = p (feval(fd,p,varargin{:}) < geps, :);          % Keep only d<0 points
r0 = 1./feval(fh,p,varargin{:});                    % Probability to keep point
p = p(rand(size(p,1),1)<r0.^dim/max(r0)^dim,:);     % Rejection method
if ~isempty(pfix), p = setdiff(p,pfix,'rows'); end
pfix = unique(pfix,'rows');
p = [pfix; p];

N = size (p,1);

iter = 0; MaxIter = 200;
pold = inf; ttol = 0.1;
while iter < MaxIter
    iter = iter + 1;
    % 3. Get the truss and bars by the Delaunay algorithm (only for large
    % movement)
    if max(sqrt(sum((p-pold).^2,2))/h0)>ttol
        pold = p;
        t = delaunayn(p);
        pc = zeros(size(t,1),dim);  % centroid of each simplex.
        for j = 1 : dim+1
            pc = pc + p(t(:,j),:)/(dim+1);
        end
        t = t(feval(fd, pc, varargin{:})<-geps, :); % Keep interior triangles
        %  4. Describe each edge by a unique pair of nodes.
        pairs = nchoosek(1:dim+1,2);  %  [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
        % nchoosek(10,3) 为组合数 C_10^3, 而 nchoosek(1:4,2) 表示从向量
        % 1:4 中选两个数（忽略顺序）所得矩阵
        npairs= size(pairs,1); NT = size(t,1);
        totalEdge = zeros(npairs*NT,2);
        for i = 1:npairs
            totalEdge((i-1)*NT+1:i*NT,:) = t(:,pairs(i,:));
        end
        totalEdge = sort(totalEdge,2);
        bars = unique(totalEdge,'rows');         % Bars as node pairs
    end
    
    % 4. Update the points
    barvec = p(bars(:,2),:)-p(bars(:,1),:);            % List of bar vectors
    L = sqrt(sum(barvec.^2,2));                        % Bar lengths
    pm = (p(bars(:,1),:)+p(bars(:,2),:))/2;
    hbars = feval(fh,pm);
    L0 = Fscale*hbars*(sum(L.^dim)/sum(hbars.^dim))^(1/dim); % Desired lengths
    f = max(L0-L,0);                                         % Bar forces (scalars)
    F = barvec.*repmat(f./L,1,dim);                          % Bar forces (x,y,z components)
    % Force resultant
    Ftot = full(sparse(bars(:,[ones(1,dim),2*ones(1,dim)]), ...
        ones(size(L))*[1:dim,1:dim], [-F,F],N,dim));
    Ftot(1:size(pfix,1),:) = 0;
    p = p + dt*Ftot;
    
    % 5. Bring outside points back to the boundary
    %  simply a gradient step (twice).
        for steps = 1:2
            d = feval(fd,p,varargin{:}); ix=d>0;       % Find points outside (d>0)
            nvec = zeros (sum(ix),dim);      %n1,n2,...
            for j = 1:dim  % nj
                pdx = p(ix,:); pdx(:,j) = pdx(:,j)+deps;
                nvec(:,j) = (feval(fd,pdx,varargin{:})-d(ix))/deps;
            end
            nn = (sum(nvec.^2,2)); % better performance without sqrt
            nn(nn<deps) = 1;
            nvec = nvec./repmat(nn,1,dim);     
            p(ix,:) = p(ix,:) - d(ix) * ones(1,dim).*nvec;
        end
        
    % 6. Termination criterion: All interior nodes move less than dptol (scaled)
    if max(sqrt(sum(dt*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end    
end

% final mesh
t = delaunayn(p);
pc = zeros(size(t,1),dim);  
for j = 1 : dim+1
    pc = pc + p(t(:,j),:)/(dim+1);
end
t = t(feval(fd, pc, varargin{:})<-geps, :); 