function err = getH1error(node,elem,Du,uh,quadOrder)
%getL2error gets H1 semi-norm of the approximation error.
%
% Copyright (C) Long Chen, modified by Terence Yu.

%% Numbers
Nu = size(uh,1);  
% auxstructure
auxT = auxstructure(node,elem);
N = size(node,1);  NT = size(elem,1);  NE = size(auxT.edge,1);
NP1 = N;   NP2 = N + NE;   NP3 = N + 2*NE + NT;  NCR = NE;

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
% auxstructure
auxT = auxstructure(node,elem);
% area
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
xi = [z2(:,1)-z3(:,1), z3(:,1)-z1(:,1), z1(:,1)-z2(:,1)];
eta = [z2(:,2)-z3(:,2), z3(:,2)-z1(:,2), z1(:,2)-z2(:,2)];
area = 0.5*(xi(:,1).*eta(:,2)-xi(:,2).*eta(:,1));
% gradbasis
Dlambdax = eta./repmat(2*area,1,3);
Dlambday = -xi./repmat(2*area,1,3);
Dlambda1 = [Dlambdax(:,1), Dlambday(:,1)];
Dlambda2 = [Dlambdax(:,2), Dlambday(:,2)];
Dlambda3 = [Dlambdax(:,3), Dlambday(:,3)];

%% P1-Lagrange
if Nu == NP1
    % elementwise d.o.f.s
    elem2dof = elem;
    % numerical gradient  (at p-th quadrature point)
    Duh =  repmat(uh(elem2dof(:,1)),1,2).*Dlambda1 ...
        +  repmat(uh(elem2dof(:,2)),1,2).*Dlambda2 ...
        +  repmat(uh(elem2dof(:,3)),1,2).*Dlambda3;
    % elementwise error
    err = zeros(NT,1);
    for p = 1:nQuad
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end
end

%% P2-Lagrange
if Nu == NP2
    % elementwise d.o.f.s
    elem2edge = auxT.elem2edge;
    elem2dof = [elem, elem2edge+N];
    % numerical gradient (at p-th quadrature point)
    err = zeros(NT,1);
    for p = 1:nQuad
        Dphip1 = (4*lambda(p,1)-1)*Dlambda1;
        Dphip2 = (4*lambda(p,2)-1)*Dlambda2;
        Dphip3 = (4*lambda(p,3)-1)*Dlambda3;
        Dphip4 = 4*(lambda(p,2)*Dlambda3+lambda(p,3)*Dlambda2);
        Dphip5 = 4*(lambda(p,3)*Dlambda1+lambda(p,1)*Dlambda3);
        Dphip6 = 4*(lambda(p,1)*Dlambda2+lambda(p,2)*Dlambda1);
        Duh =  repmat(uh(elem2dof(:,1)),1,2).*Dphip1 ...
            +  repmat(uh(elem2dof(:,2)),1,2).*Dphip2 ...
            +  repmat(uh(elem2dof(:,3)),1,2).*Dphip3 ...
            +  repmat(uh(elem2dof(:,4)),1,2).*Dphip4 ...
            +  repmat(uh(elem2dof(:,5)),1,2).*Dphip5 ...
            +  repmat(uh(elem2dof(:,6)),1,2).*Dphip6;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
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
    % elementwise d.o.f.s
    elem2dof = [elem, elema, elemb, (1:NT)'+N+2*NE];
    % numerical gradient (at p-th quadrature point)
    err = zeros(NT,1);
    for p = 1:nQuad
        Dphip1 = (27/2*lambda(p,1)*lambda(p,1)-9*lambda(p,1)+1).*Dlambda1;           
        Dphip2 = (27/2*lambda(p,2)*lambda(p,2)-9*lambda(p,2)+1).*Dlambda2; 
        Dphip3 = (27/2*lambda(p,3)*lambda(p,3)-9*lambda(p,3)+1).*Dlambda3;
        Dphip4 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda3+...
                lambda(p,3)*(6*lambda(p,2)-1).*Dlambda2);
        Dphip5 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda1+...
                 lambda(p,1)*(6*lambda(p,3)-1).*Dlambda3);        
        Dphip6 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda2+...
                 lambda(p,2)*(6*lambda(p,1)-1).*Dlambda1);
        Dphip7 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda2+...
                 lambda(p,2)*(6*lambda(p,3)-1).*Dlambda3);
        Dphip8 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda3+...
                 lambda(p,3)*(6*lambda(p,1)-1).*Dlambda1);  
        Dphip9 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda1+...
                 lambda(p,1)*(6*lambda(p,2)-1).*Dlambda2);  
        Dphip10 = 27*(lambda(p,1)*lambda(p,2)*Dlambda3+lambda(p,1)*lambda(p,3)*Dlambda2+...
                 lambda(p,3)*lambda(p,2)*Dlambda1);
        Duh =  repmat(uh(elem2dof(:,1)),1,2).*Dphip1 ...
            +  repmat(uh(elem2dof(:,2)),1,2).*Dphip2 ...
            +  repmat(uh(elem2dof(:,3)),1,2).*Dphip3 ...
            +  repmat(uh(elem2dof(:,4)),1,2).*Dphip4 ...
            +  repmat(uh(elem2dof(:,5)),1,2).*Dphip5 ...
            +  repmat(uh(elem2dof(:,6)),1,2).*Dphip6 ...
            +  repmat(uh(elem2dof(:,7)),1,2).*Dphip7 ...
            +  repmat(uh(elem2dof(:,8)),1,2).*Dphip8 ...
            +  repmat(uh(elem2dof(:,9)),1,2).*Dphip9 ...
            +  repmat(uh(elem2dof(:,10)),1,2).*Dphip10;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
end

%% Crouzeix-Raviart linear element
if Nu == NCR
    % elementwise d.o.f.s
    elem2dof = auxT.elem2edge;
    % numerical gradient  (at p-th quadrature point)
    Duh =  repmat(uh(elem2dof(:,1)),1,2).*(-2*Dlambda1) ...
        +  repmat(uh(elem2dof(:,2)),1,2).*(-2*Dlambda2) ...
        +  repmat(uh(elem2dof(:,3)),1,2).*(-2*Dlambda3);
    % elementwise error
    err = zeros(NT,1);
    for p = 1:nQuad
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2 + lambda(p,3)*z3;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end
end

%% Modification
err = area.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(abs(sum(err)));