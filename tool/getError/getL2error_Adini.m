function err = getL2error_Adini(node,elem,u,uh,quadOrder)

if ~exist('quadOrder','var') || isempty(quadOrder)
    quadOrder = 3;
end

%% elem2dof
N = size(node,1); NT = size(elem,1); Ndof = 4*3;
elem2dof = [elem, elem+N, elem + 2*N];

%% parameters used in computation
% area
xiv = [-1 1 1 -1]; etav = [-1 -1 1 1]; % vertices of [-1,1]^2
x1 = node(elem(:,1),1); y1 = node(elem(:,1),2);
x2 = node(elem(:,2),1); y4 = node(elem(:,4),2);
a = x2-x1;  b = y4-y1;
area = a.*b;
% quadrature points on all elements [xa, xb]
[lambda,weight] = quadpts1(quadOrder);
xx = (2*lambda(:,1)-1)';  % each row corresponds to an element
ww = weight;
nr = length(ww); nI = nr^2;
xi = reshape(repmat(xx',nr,1),1,[]);
eta = reshape(ones(nr,1)*xx,1,[]);
weight = reshape(ww'*ww,1,[]);

%% compute L2 error
err = zeros(NT,1);
x0 = (x1+x2)./2; y0 = (y1+y4)./2;
for p = 1:nI
    % quadrature points in the x-y coordinate
    x = a*xi(p)/2 + x0; y = b*eta(p)/2 + y0;
    pz = [x,y];
     % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof);
    for i = 1:4
        tp = (1+xiv(i)*xi(p))*(1+etav(i)*eta(p));
        % \phi
        base(:,i) = 1/8*tp*(2+xiv(i)*xi(p)+etav(i)*eta(p)-xi(p)^2-eta(p)^2);
        % \psi
        base(:,4+i) = -1/16*a*xiv(i)*tp*(1-xi(p)^2);
        % \zeta
        base(:,8+i) = -1/16*b*etav(i)*tp*(1-eta(p)^2);
    end
    uhp = 0;
    for i = 1:Ndof
        uhp = uhp + uh(elem2dof(:,i)).*base(:,i);
    end
    err = err + weight(p)*sum((u(pz)-uhp).^2,2);
end
err = (area./4).*err;

%% Modification
err(isnan(err)) = 0; % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(abs(sum(err)));