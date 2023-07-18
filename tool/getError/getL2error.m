function err = getL2error(node,elem,u,uh,quadOrder)
%getL2error gets L2 norm of the approximation error.
%
% Copyright (C) Long Chen, modified by Terence Yu.

%% Numbers
Nu = size(uh,1);  
% auxstructure
auxT = auxstructure(node,elem);
N = size(node,1);  NT = size(elem,1);  NE = size(auxT.edge,1);
NP1 = N;    NP2 = N + NE;   NP3 = N + 2*NE + NT;  NCR = NE;  

%% Default quadrature orders for different elements
if ~exist('quadOrder','var')
    switch Nu
        case NP1      % piecewise linear function P1 element
            quadOrder = 3;
        case NP2    % piecewise quadratic function
            quadOrder = 4;     
        case NP3    % P3 element
            quadOrder = 5;
        case NCR    % Crouzeix-Raviart linear element
            quadOrder = 3; 
    end
end

% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder); nQuad = length(weight);
% area of triangles
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
ve2 = z1-z3; ve3 = z2-z1;
area = 0.5*abs(-ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1));

%% P1-Lagrange
if Nu == NP1
    elem2dof = elem;
    phi = lambda; % basis functions
end

%% P2-Lagrange
if Nu == NP2
    elem2edge = auxT.elem2edge;
    elem2dof = [elem, elem2edge+N];
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = 4*lambda(:,2).*lambda(:,3);
    phi(:,5) = 4*lambda(:,1).*lambda(:,3);
    phi(:,6) = 4*lambda(:,2).*lambda(:,1);
end

%% P3-Lagrange
if Nu == NP3
    elem2edge = auxT.elem2edge;
    % sgnelem
    edge = auxT.edge; NE = size(edge,1);
    bdStruct = setboundary(node,elem);
    v1 = [2 3 1]; v2 = [3 1 2]; % L1: 2-3
    bdEdgeIdx = bdStruct.bdEdgeIdx; E = false(NE,1); E(bdEdgeIdx) = 1;
    sgnelem = sign(elem(:,v2)-elem(:,v1));
    sgnbd = E(elem2edge);    sgnelem(sgnbd) = 1;
    sgnelem(sgnelem==-1) = 0;
    elema = elem2edge + N*sgnelem + (N+NE)*(~sgnelem); % 1/3 point
    elemb = elem2edge + (N+NE)*sgnelem + N*(~sgnelem); % 2/3 point
    % local --> global
    elem2dof = [elem, elema, elemb, (1:NT)'+N+2*NE];
    phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
    phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
    phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
    phi(:,4) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1);
    phi(:,5) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);
    phi(:,6) = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,1)-1);
    phi(:,7) = 9/2*lambda(:,2).*lambda(:,3).*(3*lambda(:,3)-1);
    phi(:,8) = 9/2*lambda(:,3).*lambda(:,1).*(3*lambda(:,1)-1);
    phi(:,9) = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,2)-1);
    phi(:,10) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);
end

%% Crouzeix-Raviart linear element
if Nu == NCR
    elem2dof = auxT.elem2edge; 
    phi = 1-2*lambda; % basis functions
end

Ndof = size(elem2dof,2);

%% elementwise error
err = zeros(NT,1);
for p = 1:nQuad
    uhp = 0;
    for j = 1:Ndof
        uhp = uhp + uh(elem2dof(:,j))*phi(p,j);
    end
    % quadrature points in the x-y coordinate
    pz =  lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
    err = err + weight(p)*(u(pz) - uhp).^2;
end
err = area.*err;
% Modification
err(isnan(err)) = 0; % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(abs(sum(err)));