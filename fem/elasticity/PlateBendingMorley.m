function [w,he] = PlateBendingMorley(node,elem,pde,bdStruct)
%PlateBendingMorley solves plate bending problem using Morley element,
% a complete quadratic triangular element, where,
% the global sign is resolved by using signed stiffness matrix and load
% vector.
%
%       -D_{ij} M_{ij}(w) + cw = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%
%   References
%   K. Feng and Z.C. Shi. "Mathematical Theory of Elastic Structures", Springer-Verlag, 1996.
%
% Copyright (C)  Terence Yu. 

para = pde.para; f = pde.f; D = para.D; nu = para.nu;

%% elem2dof
% numbers
N = size(node,1); NT = size(elem,1); 
% edge and elem2edge
allEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdge = sort(allEdge,2);
[edge,~,totalJ] = unique(totalEdge, 'rows');
elem2edge  = reshape(totalJ, NT, 3);
% elem2dof
elem2dof = [elem, elem2edge+N];
% dof numbers
NE = size(edge,1);
Ndof = 6;  NNdof = N + NE;

%% Preparation for the computation
% area
z1 = node(elem(:,1),:); 
z2 = node(elem(:,2),:); 
z3 = node(elem(:,3),:);
e1 = z2-z3; e2 = z3-z1; e3 = z1-z2; % ei = [xi, etai]
area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
% xi, eta
xi = [e1(:,1), e2(:,1), e3(:,1)]; 
eta = [e1(:,2), e2(:,2), e3(:,2)]; 
% sij
sb = zeros(NT,3,3);
for i = 1:3
    for j = 1:3
        sb(:,i,j) = -xi(:,i).*xi(:,j) - eta(:,i).*eta(:,j);
    end
end
% elementwise edge length
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z2-z1).^2,2)); L = he(elem2edge);
% elementwise sign of basis functions
sgnedge = sign([elem(:,3)-elem(:,2), elem(:,1)-elem(:,3), elem(:,2)-elem(:,1)]);
bdEdgeIdx = bdStruct.bdEdgeIdx;
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnedge(sgnbd) = 1;
sgnbase = ones(NT,Ndof); sgnbase(:,4:6) = sgnedge;
% coefficients in the basis functions
ind = [1 2 3; 2 3 1; 3 1 2];  % rotation index
it = ind(:,1); jt = ind(:,2); kt = ind(:,3);
c0 = zeros(NT,3); c1 = c0; c2 = c0;
for i = 1:3
    j = ind(i,2); k = ind(i,3);
    c0(:,i) = 1./(2*area.^2);
    c1(:,i) = sb(:,i,j)./L(:,j).^2;
    c2(:,i) = sb(:,k,i)./L(:,k).^2;
end
c3 = c1+c2;
% Gauss quadrature rule
[lambda,weight] = quadpts(3);
nQuad = length(weight);

%% Second derivatives of basis functions
b11 = zeros(NT,6); b22 = b11; b12 = b11; % [phi_i, i=1:6]
% i = 1,2,3
b11(:,it) = c0.* (eta(:,it).^2 - c1.*eta(:,jt).^2 - c2.*eta(:,kt).^2);
b22(:,it) = c0.* (xi(:,it).^2 - c1.*xi(:,jt).^2 - c2.*xi(:,kt).^2);
b12(:,it) = -c0.* (xi(:,it).*eta(:,it) - c1.*xi(:,jt).*eta(:,jt) - c2.*xi(:,kt).*eta(:,kt));
% i = 4,5,6
ci = 1./(repmat(area,1,3).*L(:,it));
b11(:,3+it) = ci.*eta(:,it).^2;
b22(:,3+it) = ci.*xi(:,it).^2;
b12(:,3+it) = -ci.*xi(:,it).*eta(:,it);

%% First stiffness matrix and sign matrix
K = zeros(NT,Ndof^2);  sgnK = zeros(NT,Ndof^2);
s = 1;
for i = 1:Ndof
    for j = 1:Ndof
        K(:,s) = b11(:,i).*b11(:,j) + b22(:,i).*b22(:,j) ...
            + nu*(b11(:,i).*b22(:,j) + b22(:,i).*b11(:,j)) ...
            + 2*(1-nu)*b12(:,i).*b12(:,j);
        sgnK(:,s) = sgnbase(:,i).*sgnbase(:,j);
        s = s+1;
    end
end
K = D*repmat(area,1,Ndof^2).*K;

%% Second stiffness matrix and load vector
G = zeros(NT,Ndof^2); F = zeros(NT,Ndof);
if isnumeric(para.c), cf = @(xy) para.c+0*xy(:,1); end
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof);
    % i = 1,2,3    
    base(:,it) = repmat(lambda(p,it).^2,NT,1) ...
               + c3.*repmat(lambda(p,jt).*lambda(p,kt),NT,1) ...
               + c2.*repmat(lambda(p,kt).*lambda(p,it),NT,1) ...
               + c1.*repmat(lambda(p,it).*lambda(p,jt),NT,1);
    % i = 4,5,6
    ci = 2*repmat(area,1,3)./L(:,it);
    base(:,3+it) = ci.*repmat(lambda(p,it).*(lambda(p,it)-1),NT,1);  
    % Second stiffness matrix
    s = 1;
    for i = 1:Ndof
        for j = 1:Ndof
            gs = cf(pxy).*base(:,i).*base(:,j);
            G(:,s) = G(:,s) + weight(p)*gs;
            s = s+1;
        end
    end
    % load vector
    F = F + weight(p)*repmat(f(pxy),1,Ndof).*base;
end
G = repmat(area,1,Ndof^2).*G; 
F = repmat(area,1,Ndof).*F;
sgnF = ones(NT,Ndof); sgnF(:,4:6) = sgnedge; % sign vector

%% Assemble stiffness matrix and load vector
ii = reshape(repmat(elem2dof, Ndof,1), [], 1);
jj = repmat(elem2dof(:), Ndof, 1);
kk = sparse(ii,jj,(K(:)+G(:)).*sgnK(:),NNdof,NNdof);
ff = accumarray(elem2dof(:), F(:).*sgnF(:), [NNdof 1]);

%% Dirichlet boundary conditions
bdNodeIdxD = bdStruct.bdNodeIdxD;
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Dw;

fixedNode = [bdNodeIdxD; bdEdgeIdx+N];
isBdNode = false(NNdof,1); isBdNode(fixedNode) = true;
freeDof = ~isBdNode;

pD = node(bdNodeIdxD,:); wD = g_D(pD);
z1 = node(bdEdgeD(:,1),:); 
z2 = node(bdEdgeD(:,2),:); 
zc = (z1+z2)./2;
e = z1-z2;  % e = z2-z1
ne = [-e(:,2),e(:,1)];  % scaled
ne = ne./repmat(he(bdEdgeIdx),1,2);
wnD = sum(Dw(zc).*ne,2);
w = zeros(NNdof,1); w(fixedNode) = [wD;wnD];
ff = ff - kk*w;

%% Solver
w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);