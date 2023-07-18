function err = getL2error1(node,elem,u,uh,quadOrder)

% length
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); 
h = sqrt(sum((z2-z1).^2,2));

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
    elem2dof = elem;
    phi = lambda; % basis functions
end

%% P2-Lagrange
if Nu == NP2
    elem2dof = [elem, (1:nel)'+N];
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = 4*lambda(:,1).*lambda(:,2);
end

%% P3-Lagrange
if Nu == NP3
    elem2dof = [elem, (1:nel)'+N, (1:nel)'+N+nel];
    phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);
    phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2);
    phi(:,3) = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,1)-1);
    phi(:,4) = 9/2*lambda(:,2).*lambda(:,1).*(3*lambda(:,2)-1);
end

Ndof = size(elem2dof,2);

%% elementwise error
err = zeros(nel,1);
for p = 1:ng
    uhp = 0;
    for j = 1:Ndof
        uhp = uhp + uh(elem2dof(:,j))*phi(p,j);
    end
    % quadrature points in the x-y coordinate
    pz =  lambda(p,1)*z1 + lambda(p,2)*z2;
    err = err + weight(p)*(u(pz) - uhp).^2;
end
err = h.*err;
% Modification
err(isnan(err)) = 0; % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(abs(sum(err)));