function u = elasticityVector(node,elem,pde,bdStruct)
%elasticityVector  Conforming P1 FEM of linear elasticity equation
% Programming in the vectorized finite element space
%
%       u = [u1, u2]
%       -div (sigma) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D
%       Neumann boundary condition   \sigma*n = g  on \Gamma_N
%       \sigma = (sigma_{ij}): stress tensor, 1<=i,j<=2
%
% Copyright (C) Terence Yu.

N = size(node,1); NT = size(elem,1); Ndof = 3; Ndof2 = 2*Ndof;
mu = pde.mu; lambda = pde.lambda; f = pde.f;

%% Compute (Dibase,Djbase)
[Dphi,area] = gradbasis(node,elem);
Dbase = cell(3,3);
for i = 1:2
    for j = 1:2
        Dbase{1,1}{i,j} = Dphi(:,i,1).*Dphi(:,j,1).*area; % b_st ---> (phi_s, phi_t)
        Dbase{1,2}{i,j} = Dphi(:,i,1).*Dphi(:,j,2).*area;
        Dbase{1,3}{i,j} = Dphi(:,i,1).*Dphi(:,j,3).*area;
        Dbase{2,1}{i,j} = Dphi(:,i,2).*Dphi(:,j,1).*area;
        Dbase{2,2}{i,j} = Dphi(:,i,2).*Dphi(:,j,2).*area;
        Dbase{2,3}{i,j} = Dphi(:,i,2).*Dphi(:,j,3).*area;
        Dbase{3,1}{i,j} = Dphi(:,i,3).*Dphi(:,j,1).*area;
        Dbase{3,2}{i,j} = Dphi(:,i,3).*Dphi(:,j,2).*area;
        Dbase{3,3}{i,j} = Dphi(:,i,3).*Dphi(:,j,3).*area;
    end
end

%% Get sparse assembling index
elem2dof = zeros(NT,Ndof2);
elem2dof(:,[1 3 5]) = 2*elem-1; elem2dof(:,[2 4 6]) = 2*elem;
ii = reshape(repmat(elem2dof, Ndof2,1), [], 1);
jj = repmat(elem2dof(:), Ndof2, 1);

%% Assemble stiffness matrix
K = zeros(NT,Ndof2^2);
s1 = 1; s2 = s1+6; 
for s = 1:3
    for t = 1:3
        k11 = (lambda+2*mu)*Dbase{s,t}{1,1} + mu*Dbase{s,t}{2,2};
        k12 = lambda*Dbase{s,t}{1,2} + mu*Dbase{s,t}{2,1};
        k21 = lambda*Dbase{s,t}{2,1} + mu*Dbase{s,t}{1,2};
        k22 = (lambda+2*mu)*Dbase{s,t}{2,2} + mu*Dbase{s,t}{1,1};    
        K(:,s1) = k11; K(:,s1+1) = k12; s1 = s1+2;
        K(:,s2) = k21; K(:,s2+1) = k22; s2 = s2+2;
    end
    s1 = s1+Ndof2; s2 = s2+Ndof2;
end
kk = sparse(ii,jj,K(:),2*N,2*N);

%% Assemble load vector
% Gauss quadrature rule
[lambda,weight] = quadpts(2);
F = zeros(NT,Ndof2); 
for iel = 1:NT
    vK = node(elem(iel,:),:); % vertices of K
    xy = lambda*vK;  fxy = f(xy); % fxy = [f1xy,f2xy]
    fv1 = fxy.*[lambda(:,1),lambda(:,1)]; % (f,phi1)
    fv2 = fxy.*[lambda(:,2),lambda(:,2)]; % (f,phi2)
    fv3 = fxy.*[lambda(:,3),lambda(:,3)]; % (f,phi3)
    
    F(iel,1:2) = area(iel)*weight*fv1;
    F(iel,3:4) = area(iel)*weight*fv2;
    F(iel,5:6) = area(iel)*weight*fv3;
end
ff = accumarray(elem2dof(:),F(:),[2*N,1]);

%% Assembl Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(ne.*Sig1(:,[1,3]),2)./2; F12 = sum(ne.*Sig2(:,[1,3]),2)./2;
    F21 = sum(ne.*Sig1(:,[3,2]),2)./2; F22 = sum(ne.*Sig2(:,[3,2]),2)./2;
    FN = [F11,F12,F21,F22];
    elem2dofN = [2*bdEdgeN(:)-1; 2*bdEdgeN(:)];
    ff = ff + accumarray(elem2dofN, FN(:),[2*N 1]);
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(N,1); isBdNode(bdNodeIdx) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
nodeD = node(bdNode,:);
bdDof = [2*bdNode-1; 2*bdNode]; freeDof = [2*freeNode-1;2*freeNode];
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);