function [u,info] = Poisson(node,elem,pde,bdStruct,option)
% Poisson solves Poisson equation with P1 linear element (2D).
%
%    - Delta u = f   in \Omega, with
%    Dirichlet boundary conditions u = g_D on \Gamma_D,
%    Neumann boundary conditions   grad(u)*n = g_N on \Gamma_N.
%
% Copyright (C) Long Chen, modified by Terence Yu.

%% Input check
if ~exist('option','var'), option = []; end

N = size(node,1); NT = size(elem,1); Ndof = 3;
f = pde.f;

%% Assemble stiffness matrix
[Dphi,area] = gradbasis(node,elem);
K = zeros(NT,Ndof^2); % straighten
s = 1;
for i = 1:Ndof
    for j = 1:Ndof
        K(:,s) = sum(Dphi(:,:,i).*Dphi(:,:,j),2).*area;
        s = s+1;
    end
end
ii = reshape(repmat(elem, Ndof,1), [], 1);
jj = repmat(elem(:), Ndof, 1);
kk = sparse(ii,jj,K(:),N,N);

%% Assemble load vector
% Gauss quadrature rule
[lambda,weight] = quadpts(2);
F = zeros(NT,3); % straighten
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fxy = f(pxy);  fv = fxy*lambda(p,:); % [f*phi1, f*phi2, f*phi3] at (xp,yp)
    F = F + weight(p)*fv;
end
F = repmat(area,1,3).*F;  % F = area.*F;
ff = accumarray(elem(:), F(:),[N 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2; % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    Du = pde.Du;
    gradu1 = Du(z1); gradu2 = Du(z2);
    F1 = sum(ne.*gradu1,2)./2; F2 = sum(ne.*gradu2,2)./2;
    FN = [F1,F2];
    ff = ff + accumarray(bdEdgeN(:), FN(:),[N 1]);
end

%% Apply Dirichlet boundary conditions
bdNodeIdxD = bdStruct.bdNodeIdxD; g_D = pde.g_D;
isBdNode = false(N,1); isBdNode(bdNodeIdxD) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdDof,:);
u = zeros(N,1); u(bdDof) = g_D(nodeD);
ff = ff - kk*u;

%% Set up solver type
if isempty(option) || ~isfield(option,'solver')  % no option.solver
    if N <= 2e3  % Direct solver for small size systems
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
        [Pro,Res] = uniformtransferoperator(elem,J);
        A = speye(N); A(freeDof,freeDof) = kk(freeDof,freeDof);
        b = u; b(freeDof) = ff(freeDof);
        u = mgVcycle(A,b,Pro,Res); % multigrid Vcycle
end

%% Store information
info.kk = kk; 
