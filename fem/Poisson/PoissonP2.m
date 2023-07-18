function u = PoissonP2(node,elem,pde,bdStruct,option)
% Poisson solves Poisson equation with P2 Lagrange element (2D).
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.
%
% Copyright (C) Long Chen, modified by Terence Yu.

%% Input check
if ~exist('option','var'), option = []; end

%% Get elem2dof
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;
elem2edge = auxT.elem2edge;
% numbers
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N + NE; Ndof = 6; %global and local d.o.f. numbers
% elem2dof
elem2dof = [elem, elem2edge + N];

%% Assemble stiffness matrix
% lambda and Dlambda
quadOrder = 4;
[lambda, weight] = quadpts(quadOrder); nG = length(weight);
[Dlambda,area] = gradbasis(node,elem);
% stiffness matrix
K = zeros(NT,Ndof^2); % straighten
for p = 1:nG
    % Dphi at quadrature points
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);
    Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
    s = 1;
    for i = 1:Ndof
        for j = 1:Ndof
            K(:,s) = K(:,s) + weight(p)*sum(Dphip(:,:,i).*Dphip(:,:,j),2).*area;
            s = s+1;
        end
    end
end
ii = reshape(repmat(elem2dof, Ndof,1), [], 1);
jj = repmat(elem2dof(:), Ndof, 1);
kk = sparse(ii,jj,K(:),NNdof,NNdof);

%% Assemble load vector
% basis functions
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
% load vector
F = zeros(NT,Ndof); % straighten
for p = 1:nG
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    F = F + weight(p)*pde.f(pxy)*phi(p,:);
end
F = repmat(area,1,Ndof).*F;  % F = area.*F;
ff = accumarray(elem2dof(:), F(:),[NNdof 1]);

%% Assemble Neumann boundary conditions
bdEdgeIdxN = bdStruct.bdEdgeIdxN; bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    % Sparse assembling index
    elem1 = [bdEdgeN, bdEdgeIdxN + N]; ndof = 3;
    % Gauss quadrature rule
    [lambda,weight] = quadpts1(quadOrder); ng = length(weight);
    % basis function
    phi1(:,3) = 4*lambda(:,1).*lambda(:,2);
    phi1(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi1(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    % nvec
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); nel = size(bdEdgeN,1);
    e = z1-z2;  he = sqrt(sum(e.^2,2));
    nvec = [-e(:,2)./he, e(:,1)./he];    
    % assemble
    FN = zeros(nel,ndof);
    for p = 1:ng
        pz = lambda(p,1)*z1 + lambda(p,2)*z2;
        Dnu = sum(pde.Du(pz).*nvec,2);
        FN = FN + weight(p)*Dnu*phi1(p,:);
    end
    FN = repmat(he,1,ndof).*FN;
    ff = ff + accumarray(elem1(:), FN(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
bdNodeIdxD = bdStruct.bdNodeIdxD; 
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
id = [bdNodeIdxD; bdEdgeIdxD+N];
g_D = pde.g_D; 
bdEdgeD = bdStruct.bdEdgeD;
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:);   
zc = (z1+z2)/2;  
nodeD = node(bdNodeIdxD,:);  
wD = g_D(nodeD); wc = g_D(zc); 
u = zeros(NNdof,1); u(bdDof) = [wD; wc];
ff = ff - kk*u;

%% Set up solver type
if isempty(option) || ~isfield(option,'solver')  % no option.solver
    if NNdof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % mg-Vcycle solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
switch solver
    case 'direct'
        u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'mg'
        disp('Multigrid V-cycle Preconditioner with Gauss-Seidel Method');
        fprintf('\n');
        if ~isfield(option,'J'), option.J = 2; end  % at least two levels
        J = option.J; 
        [Pro,Res] = uniformtransferoperatorP2(elem,J);
        A = speye(NNdof); A(freeDof,freeDof) = kk(freeDof,freeDof);
        b = u; b(freeDof) = ff(freeDof);
        u = mgVcycle(A,b,Pro,Res); % multigrid Vcycle
end
