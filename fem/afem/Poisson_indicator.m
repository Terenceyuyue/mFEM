function eta = Poisson_indicator(node,elem,uh,pde)
% This function returns the local error indicator of Poisson equation with
% homogeneous Dirichlet boundary condition in 2-D.

% Copyright (C) Terence Yu.

%% preparation for the computation
% auxiliary data
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge, ~, totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);  NE = size(edge,1);
elem2edge = reshape(totalJ,NT,3);
% area, diameter, he, ne
z1 = node(elem(:,1),:); z2 = node(elem(:,2),:); z3 = node(elem(:,3),:);
e1 = z2-z3; e2 = z3-z1; e3 = z1-z2; % ei = [xi, etai]
area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
e = node(edge(:,2),:)-node(edge(:,1),:);
he = sqrt(sum(e.^2,2));
diameter = max(he(elem2edge), [], 2);
ne = [-e(:,2),e(:,1)]; % stored in rows
ne = ne./repmat(he,1,2);
% gradbais
grad1 = [e1(:,2), -e1(:,1)]./repmat(2*area,1,2); % stored in rows
grad2 = [e2(:,2), -e2(:,1)]./repmat(2*area,1,2);
grad3 = -(grad1+grad2); 
% elem2dof
elem2dof = elem;
% sign of elementwise edges
sgnelem = sign([elem(:,3)-elem(:,2), ...
                elem(:,1)-elem(:,3), ...
                elem(:,2)-elem(:,1)]);

%% elementwise residuals
[lambda,weight] = quadpts(3);
elemRes = zeros(NT,1);
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*z1+lambda(p,2)*z2+lambda(p,3)*z3;
    fxy = pde.f(pxy);  
    elemRes = elemRes + weight(p)*fxy.^2;
end
elemRes = diameter.^2.*area.*elemRes;

%% connectivity list of jump integral
[~,weight1d] = quadpts1(4); ng = length(weight1d);
elem2dofP = zeros(NT,3*ng); % P: plus +
%elem2dofM = zeros(NT,3*ng); % M: minus -
for i = 1:3  % i-th side
    ei = elem2edge(:,i); sgni = sgnelem(:,i);
    idP = repmat(1:ng, NT,1);  
    idP(sgni<0,:) = repmat((ng:-1:1)+NE*ng, sum(sgni<0), 1);
    elem2dofP(:,(1:ng)+(i-1)*ng) = idP + (ei-1)*ng;
end
index = (elem2dofP>ng*NE);
elem2dofM = elem2dofP + (-ng*NE)*index + ng*NE*(~index);

%% elementwise interpolant at quadrature points
elemuhxP = zeros(NT,3*ng); 
elemuhyP = zeros(NT,3*ng);
basex = [grad1(:,1),grad2(:,1),grad3(:,1)]; % Dxphi1,Dxphi2,Dxphi3
basey = [grad1(:,2),grad2(:,2),grad3(:,2)];
Ndof = 3;
for p = 1:3*ng  
    % interpolation at the p-th quadrture point
    for i = 1:Ndof
        elemuhxP(:,p) = elemuhxP(:,p) + uh(elem2dof(:,i)).*basex(:,i);
        elemuhyP(:,p) = elemuhyP(:,p) + uh(elem2dof(:,i)).*basey(:,i);
    end
end
uhxI = zeros(2*NE*ng,1);      uhyI = zeros(2*NE*ng,1);
uhxI(elem2dofP) = elemuhxP;   uhyI(elem2dofP) = elemuhyP;
elemuhxM = uhxI(elem2dofM);   elemuhyM = uhyI(elem2dofM);

%% elementwise jump at quadrature points
elem2Jumpx = elemuhxP - elemuhxM;
elem2Jumpy = elemuhyP - elemuhyM;
elemJump = zeros(NT,1);
for i = 1:3
    hei = he(elem2edge(:,i));
    cei = hei;
    neix = ne(elem2edge(:,i),1); neiy = ne(elem2edge(:,i),2);
    Jumpnx = elem2Jumpx(:,(1:ng)+(i-1)*ng).*repmat(neix,1,ng);
    Jumpny = elem2Jumpy(:,(1:ng)+(i-1)*ng).*repmat(neiy,1,ng);
    Jumpn = (Jumpnx+Jumpny).^2;    
    elemJump = elemJump + cei.*hei.*(Jumpn*weight1d(:));
end

%% Local error indicator
eta = (abs(elemRes) + elemJump).^(1/2);