function err = getL2error_Zienkiewicz(node,elem,u,uh,quadOrder)

if ~exist('quadOrder','var') || isempty(quadOrder)
    quadOrder = 3;
end

%% elem2dof
N = size(node,1); NT = size(elem,1); Ndof = 9;
elem2dof = [elem, elem+N, elem+2*N];

%% parameters used in computation
% area
z1 = node(elem(:,1),:); 
z2 = node(elem(:,2),:); 
z3 = node(elem(:,3),:);
e1 = z2-z3; e2 = z3-z1; e3 = z1-z2; % ei = [xi, etai]
area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
% xi, eta
xi = [e1(:,1), e2(:,1), e3(:,1)]; 
eta = [e1(:,2), e2(:,2), e3(:,2)]; 
% coefficients in the basis functions
c1 = xi(:,1);   d1 = eta(:,1);
c2 = -xi(:,2);  d2 = -eta(:,2);
c3 = -[c2,c1,zeros(NT,1)]; d3 = -[d2,d1,zeros(NT,1)];
c4 = [0.5*c1-c2,  0.5*c2-c1,  0.5*(c1+c2)];
d4 = [0.5*d1-d2,  0.5*d2-d1,  0.5*(d1+d2)];
% Gauss quadrature rule
[lambda,weight] = quadpts(quadOrder);
nQuad = length(weight); % number of integration points

%% compute L2 error
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pz = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof);
    for i = 1:3
        % \phi
        lam123 = lambda(p,1).*lambda(p,2).*lambda(p,3);
        base(:,i) = lambda(p,i).^2.*(3-2*lambda(p,i)) + 2*lam123;
        % \psi
        base(:,3+i) = lambda(p,i).^2.*(c1*lambda(p,2) + c2*lambda(p,1) + c3(:,i)) ...
            + c4(:,i)*lam123;
        % \zeta
        base(:,6+i) = lambda(p,i).^2.*(d1*lambda(p,2) + d2*lambda(p,1) + d3(:,i)) ...
            + d4(:,i)*lam123;
    end
    uhp = 0;
    for i = 1:Ndof
        uhp = uhp + uh(elem2dof(:,i)).*base(:,i);
    end
    err = err + weight(p)*sum((u(pz)-uhp).^2,2);
end
err = area.*err;

%% Modification
err(isnan(err)) = 0; % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(abs(sum(err)));