function u = elasticityNavierCR(node,elem,pde,bdStruct)
%elasticityNavierCR solves linear elasticity equation of Navier form using Crouzeix-Raviart element 
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C) Terence Yu.

%% elem2dof of ui
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(allEdge,2);
[edge, ~, totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = reshape(totalJ,NT,3); 
elem2dof = elem2edge;

NE = size(edge,1); 
Ndof = 3;  % local dof number of ui

mu = pde.mu; lambda = pde.lambda; f = pde.f;

%% Compute (Dibase,Djbase)
[Dphi,area] = gradbasis(node,elem);
Dbase = cell(2,2);
for i = 1:2
    for j = 1:2
        k11 = 4*Dphi(:,i,1).*Dphi(:,j,1).*area;
        k12 = 4*Dphi(:,i,1).*Dphi(:,j,2).*area;
        k13 = 4*Dphi(:,i,1).*Dphi(:,j,3).*area;
        k21 = 4*Dphi(:,i,2).*Dphi(:,j,1).*area;
        k22 = 4*Dphi(:,i,2).*Dphi(:,j,2).*area;
        k23 = 4*Dphi(:,i,2).*Dphi(:,j,3).*area;
        k31 = 4*Dphi(:,i,3).*Dphi(:,j,1).*area;
        k32 = 4*Dphi(:,i,3).*Dphi(:,j,2).*area;
        k33 = 4*Dphi(:,i,3).*Dphi(:,j,3).*area;
        K = [k11,k12,k13,k21,k22,k23,k31,k32,k33]; % stored in rows
        Dbase{i,j} = K(:); % straighten
    end
end

%% Get sparse assembling index
ii = reshape(repmat(elem2dof, Ndof,1), [], 1);
jj = repmat(elem2dof(:), Ndof, 1);
ii11 = ii;     jj11 = jj;  ii12 = ii;    jj12 = jj+NE;
ii21 = ii+NE;  jj21 = jj;  ii22 = ii+NE; jj22 = jj+NE;

%% Assemble stiffness matrix
% (grad u,grad v)
ss11 = Dbase{1,1}+Dbase{2,2};  ss22 = ss11;
ii = [ii11; ii22]; jj = [jj11; jj22]; ss = [ss11; ss22];
A = sparse(ii,jj,ss,2*NE,2*NE);
A = mu*A;

% (div u,div v)
ss11 = Dbase{1,1};            ss12 = Dbase{1,2};
ss21 = Dbase{2,1};            ss22 = Dbase{2,2};
ii = [ii11; ii12; ii21; ii22];
jj = [jj11; jj12; jj21; jj22];
ss = [ss11; ss12; ss21; ss22];
B = sparse(ii,jj,ss,2*NE,2*NE);
B = (lambda+mu)*B;

% stiffness matrix
kk = A + B;

%% Assemble load vector
% Gauss quadrature rule
[lambda,weight] = quadpts(3);
F1 = zeros(NT,3); F2 = zeros(NT,3);
phi = 1-2*lambda; % basis functions for CR element
for p = 1:length(weight)
    pxy = lambda(p,1)*node(elem(:,1),:) ...
       + lambda(p,2)*node(elem(:,2),:) ...
       + lambda(p,3)*node(elem(:,3),:);
    fxy = f(pxy); % fxy = [f1xy,f2xy]
    F1 = F1 + weight(p)*fxy(:,1)*phi(p,:);
    F2 = F2 + weight(p)*fxy(:,2)*phi(p,:);
end
F1 = repmat(area,1,3).*F1;  % F = area.*F;
F2 = repmat(area,1,3).*F2;
ff = accumarray([elem2dof(:);elem2dof(:)+NE], [F1(:);F2(:)], [2*NE 1]);

%% Apply Dirichlet boundary conditions
isBdDof = false(2*NE,1);
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
fixedNode = [bdEdgeIdxD; bdEdgeIdxD+NE];
isBdDof(fixedNode) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
bdEdgeD = bdStruct.bdEdgeD;
zm = (node(bdEdgeD(:,1),:)+node(bdEdgeD(:,2),:))/2;   
u = zeros(2*NE,1); u(bdDof) = pde.g_D(zm);
ff = ff - kk*u;

%% Set up solver type
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);