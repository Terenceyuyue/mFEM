function err = getH1error1(node,elem,Du,uh,quadOrder)

% length
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); 
h = sqrt(sum((z2-z1).^2,2));
% derivatives of bases
Dlambda1 = -1./h; 
Dlambda2 = 1./h;

Nu = size(uh,1);    N = size(node,1);   nel = size(elem,1);
NP1 = N;    NP2 = N + nel;   NP3 = N + 2*nel;    

%% Default quadrature orders for different elements
if ~exist('quadOrder','var')
    switch Nu
        case NP1      % piecewise linear function P1 element
            quadOrder = 3;
        case NP2    % piecewise quadratic function
            quadOrder = 4;     
        case NP3    % P3 element
            quadOrder = 5;               
    end
end

% Gauss quadrature rule
[lambda,weight] = quadpts1(quadOrder); ng = length(weight);

%% P1-Lagrange
if Nu == NP1
    % elementwise d.o.f.s
    elem2dof = elem;
    % numerical gradient  (at p-th quadrature point)
    Duh =  uh(elem2dof(:,1)).*Dlambda1 ...
        +  uh(elem2dof(:,2)).*Dlambda2;
    % elementwise error
    err = zeros(nel,1);
    for p = 1:ng
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end
end

%% P2-Lagrange
if Nu == NP2
    % elementwise d.o.f.s
    elem2dof = [elem, (1:nel)'+N];
    % numerical gradient (at p-th quadrature point)
    err = zeros(nel,1);
    for p = 1:ng
        Dphip1 = (4*lambda(p,1)-1)*Dlambda1;
        Dphip2 = (4*lambda(p,2)-1)*Dlambda2;
        Dphip3 = 4*(lambda(p,1)*Dlambda2+lambda(p,2)*Dlambda1);
        Duh =  uh(elem2dof(:,1)).*Dphip1 ...
            +  uh(elem2dof(:,2)).*Dphip2 ...
            +  uh(elem2dof(:,3)).*Dphip3;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
end

%% P3-Lagrange
if Nu == NP3
    % elementwise d.o.f.s
    elem2dof = [elem, (1:nel)'+N, (1:nel)'+N+nel];
    % numerical gradient (at p-th quadrature point)
    err = zeros(nel,1);
    for p = 1:ng
        Dphip1 = (27/2*lambda(p,1)^2-9*lambda(p,1)+1).*Dlambda1;           
        Dphip2 = (27/2*lambda(p,2)^2-9*lambda(p,2)+1).*Dlambda2;        
        Dphip3 = 9/2*(lambda(p,1)*(3*lambda(p,1)-1).*Dlambda2 ...
                    + lambda(p,2)*(6*lambda(p,1)-1).*Dlambda1);       
        Dphip4 = 9/2*(lambda(p,2)*(3*lambda(p,2)-1).*Dlambda1 ...
                    + lambda(p,1)*(6*lambda(p,2)-1).*Dlambda2);        
        Duh =  uh(elem2dof(:,1)).*Dphip1 ...
            +  uh(elem2dof(:,2)).*Dphip2 ...
            +  uh(elem2dof(:,3)).*Dphip3 ...
            +  uh(elem2dof(:,4)).*Dphip4;
        % elementwise error
        pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
        err = err + weight(p)*sum((Du(pz)-Duh).^2,2);
    end    
end

%% Modification
err = h.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(abs(sum(err)));