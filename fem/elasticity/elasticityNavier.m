function u = elasticityNavier(node,elem,pde,bdStruct,option)
%elasticityNavier solves linear elasticity equation of Navier form using P1 element 
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C) Terence Yu.

%% Input check
if ~exist('option','var'), option = []; end

N = size(node,1); NT = size(elem,1); Ndof = 3;
mu = pde.mu; lambda = pde.lambda; f = pde.f;

%% Compute (Dibase,Djbase)
[Dphi,area] = gradbasis(node,elem);
Dbase = cell(2,2);
for i = 1:2
    for j = 1:2
        k11 = Dphi(:,i,1).*Dphi(:,j,1).*area;
        k12 = Dphi(:,i,1).*Dphi(:,j,2).*area;
        k13 = Dphi(:,i,1).*Dphi(:,j,3).*area;
        k21 = Dphi(:,i,2).*Dphi(:,j,1).*area;
        k22 = Dphi(:,i,2).*Dphi(:,j,2).*area;
        k23 = Dphi(:,i,2).*Dphi(:,j,3).*area;
        k31 = Dphi(:,i,3).*Dphi(:,j,1).*area;
        k32 = Dphi(:,i,3).*Dphi(:,j,2).*area;
        k33 = Dphi(:,i,3).*Dphi(:,j,3).*area;
        K = [k11,k12,k13,k21,k22,k23,k31,k32,k33]; % stored in rows
        Dbase{i,j} = K(:); % straighten
    end
end

%% Get sparse assembling index
ii = reshape(repmat(elem, Ndof,1), [], 1);
jj = repmat(elem(:), Ndof, 1);
ii11 = ii;   jj11 = jj;  ii12 = ii;   jj12 = jj+N;
ii21 = ii+N; jj21 = jj;  ii22 = ii+N; jj22 = jj+N;

%% Assemble stiffness matrix
% (grad u,grad v)
ss11 = Dbase{1,1}+Dbase{2,2};  ss22 = ss11;
ii = [ii11; ii22]; jj = [jj11; jj22]; ss = [ss11; ss22];
A = sparse(ii,jj,ss,2*N,2*N);
A = mu*A;

% (div u,div v)
ss11 = Dbase{1,1};            ss12 = Dbase{1,2};
ss21 = Dbase{2,1};            ss22 = Dbase{2,2};
ii = [ii11; ii12; ii21; ii22];
jj = [jj11; jj12; jj21; jj22];
ss = [ss11; ss12; ss21; ss22];
B = sparse(ii,jj,ss,2*N,2*N);
B = (lambda+mu)*B;

% stiffness matrix
kk = A + B;

%% Assemble load vector
% Gauss quadrature rule
[lambda,weight] = quadpts(3);
F1 = zeros(NT,3); F2 = zeros(NT,3);
for p = 1:length(weight)
    pxy = lambda(p,1)*node(elem(:,1),:) ...
       + lambda(p,2)*node(elem(:,2),:) ...
       + lambda(p,3)*node(elem(:,3),:);
    fxy = f(pxy); % fxy = [f1xy,f2xy]
    F1 = F1 + weight(p)*fxy(:,1)*lambda(p,:);
    F2 = F2 + weight(p)*fxy(:,2)*lambda(p,:);
end
F1 = repmat(area,1,3).*F1;  % F = area.*F;
F2 = repmat(area,1,3).*F2;
ff = accumarray([elem(:);elem(:)+N], [F1(:);F2(:)], [2*N 1]);

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  bdNodeIdxD = bdStruct.bdNodeIdxD;
id = [bdNodeIdxD; bdNodeIdxD+N]; 
isBdNode = false(2*N,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdxD,:);
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set up solver type
if isempty(option) || ~isfield(option,'solver')  % no option.solver
    if 2*N <= 2e3  % Direct solver for small size systems
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
        if ~isfield(option,'J'), option.J = 2; end % at least two levels
        J = option.J; 
        [Pro,Res] = uniformtransferoperator(elem,J);  
        fun = @(mat) blkdiag(mat,mat);
        Pro = cellfun(fun, Pro, 'UniformOutput', false);
        Res = cellfun(fun, Res, 'UniformOutput', false);
        A = speye(2*N); A(freeDof,freeDof) = kk(freeDof,freeDof);
        b = u; b(freeDof) = ff(freeDof);
        u = mgVcycle(A,b,Pro,Res); % multigrid Vcycle
end