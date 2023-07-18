function u = FEM1DVcycle(node,elem,pde,bdStruct,Pro,Res)
%FEM1DVcycle solves the 1-D partial diffential equation by Multigrid
% V-cycle method
%
%  -au''+bu+cu = f,  x \in (x0,xL);
%     g_N = du(x0), g_D = u(xL);
%  or g_D = u(x0),  g_N = du(xL);
%
% where, g_N for Neumann boundary conditions and
%        g_D for Dirichlet boundary conditions,
%        a, b and c are constants.
%
% Copyright (C) Terence Yu.

N = size(node,1); Ndof = 2;

%% Assemble stiffness matrix
% All element matrices
para = pde.para;
a = para.a; b = para.b; c = para.c;
x1 = node(elem(:,1)); x2 = node(elem(:,2));
h = x2-x1;
k11 = a./h+b/2*(-1)+c.*h./6*2;
k12 = a./h*(-1)+b/2+c.*h./6;
k21 = a./h*(-1)+b/2*(-1)+c.*h./6;
k22 = a./h+b/2+c.*h./6*2;
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

%% Set solver
disp('Multigrid V-cycle Preconditioner with Gauss-Seidel Method');
fprintf('\n');
A = speye(N); A(freeDof,freeDof) = kk(freeDof,freeDof);
b = u; b(freeDof) = ff(freeDof);
u = mgVcycle(A,b,Pro,Res); % multigrid Vcycle