function u = FEM1D(node,elem,pde,bdStruct,option)
%FEM1D solves the two-point boundary value problem 
%  
%     -au'' + bu' + cu = f  in \Omega = (a,b), with
%     Dirichlet boundary condition u = g_D  on \Gamma_D = {a}, {b} or {a,b},
%     Neumann boundary condition   u = u' on \Gamma_N = {a,b} - \Gamma_D
%
% Copyright (C) Terence Yu. 

%% Input check
if ~exist('option','var'), option = []; end

N = size(node,1);  Ndof = 2;

%% Assemble stiffness matrix
% All element matrices 
para = pde.para;
a = para.a; b = para.b; c = para.c;
x1 = node(elem(:,1));  x2 = node(elem(:,2));
h = x2-x1;
k11 = a./h+b/2*(-1)+c*h./6*2;      
k12 = a./h*(-1)+b/2+c*h./6;
k21 = a./h*(-1)+b/2*(-1)+c*h./6;   
k22 = a./h+b/2+c*h./6*2;
K = [k11,k12,k21,k22]; % stored in rows
% stiffness matrix
ii = reshape(repmat(elem, Ndof,1), [], 1);
jj = repmat(elem(:), Ndof, 1);
kk = sparse(ii,jj,K(:),N,N);

%% Assemble load vector
xc = (x1+x2)./2; 
F1 = pde.f(xc).*h./2; F2 = F1; F = [F1,F2];
ff = accumarray(elem(:), F(:),[N 1]);

%% Assemble Neumann boundary conditions
Neumann = bdStruct.Neumann; 
if ~isempty(Neumann)
    nvec = 1;
    if find(elem(:,1)==Neumann), nvec = -1; end
    Dnu = pde.Du(node(Neumann,:))*nvec;
    ff(Neumann) = ff(Neumann) + a*Dnu;
end

%% Apply Dirichlet boundary conditions
Dirichlet = bdStruct.Dirichlet; g_D = pde.g_D;
isBdNode = false(N,1); isBdNode(Dirichlet) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
u = zeros(N,1); u(bdDof) = g_D(node(Dirichlet));
ff = ff - kk*u;

%% Set up solver type
if isempty(option) || ~isfield(option,'solver')  % no option.solver
    if N <= 1e3  % Direct solver for small size systems
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
        [Pro,Res] = uniformtransferoperator1(elem,J);
        A = speye(N); A(freeDof,freeDof) = kk(freeDof,freeDof);
        b = u; b(freeDof) = ff(freeDof);
        u = mgVcycle(A,b,Pro,Res); % multigrid Vcycle
end
