function [u,info] = PoissonCR(node,elem,pde,bdStruct)
% PoissonCR solves Poisson equation with Crouzeix-Raviart element (2D).
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.
%
% Copyright (C) Terence Yu.

%% elem2dof
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(allEdge,2);
[edge, ~, totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = reshape(totalJ,NT,3); 
elem2dof = elem2edge;

Ndof = 3; NNdof = size(edge,1); % NE
f = pde.f;

%% Assemble stiffness matrix
[Dphi,area] = gradbasis(node,elem);
K = zeros(NT,Ndof^2); % straighten
s = 1;
for i = 1:Ndof
    for j = 1:Ndof
        K(:,s) = 4*sum(Dphi(:,:,i).*Dphi(:,:,j),2).*area;
        s = s+1;
    end
end
ii = reshape(repmat(elem2dof, Ndof,1), [], 1);
jj = repmat(elem2dof(:), Ndof, 1);
kk = sparse(ii,jj,K(:),NNdof,NNdof);

%% Assemble load vector
% Gauss quadrature rule
[lambda,weight] = quadpts(3);
F = zeros(NT,3); % straighten
phi = 1-2*lambda; % basis functions for CR element
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fxy = f(pxy);  fv = fxy*phi(p,:); % [f*phi1, f*phi2, f*phi3] at (xp,yp)
    F = F + weight(p)*fv;
end
F = repmat(area,1,3).*F;  % F = area.*F;
ff = accumarray(elem2dof(:), F(:),[NNdof 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    % sparse assembling index
    elem2dof1 = bdStruct.bdElem2edgeN; 
    % Gauss quadrature rule
    [lambdaN,weightN] = quadpts1(3); ng = length(weightN);
    % basis function
    phi1(:,3) = 1-2*lambdaN(:,2)';
    phi1(:,1) = ones(ng,1);
    phi1(:,2) = 1-2*lambdaN(:,1)'; 
    % nvec
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); 
    e = z1-z2;  
    nvec = [-e(:,2), e(:,1)];   % scaled
    % assemble
    FN = zeros(size(e,1),size(phi1,2));
    for p = 1:ng
        pz = lambdaN(p,1)*z1 + lambdaN(p,2)*z2;
        Dnu = sum(pde.Du(pz).*nvec,2);
        FN = FN + weightN(p)*Dnu*phi1(p,:);
    end
    ff = ff + accumarray(elem2dof1(:), FN(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
isBdDof = false(NNdof,1);
fixedNode = bdStruct.bdEdgeIdxD;
bdEdgeD = bdStruct.bdEdgeD;
isBdDof(fixedNode) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
zm = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;
bdval = pde.g_D(zm);    
u = zeros(NNdof,1); u(bdDof) = bdval;
ff = ff - kk*u;

%% Set up solver type
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

info.elem2edge = elem2edge;