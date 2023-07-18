function eta = Poisson_indicator_old(node,elem,u,pde)
% This function returns the local error indicator of Poisson equation with
% homogeneous Dirichlet boundary condition in 2-D.

% Copyright (C) Long Chen (see estimateresidual.m), modified by Terence Yu.

%% auxiliary data
aux = auxgeometry(node,elem);
area = aux.area; diameter = aux.diameter;
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; 
edge = auxT.edge; 
edge2elem = auxT.edge2elem;

%% elementwise residuals
[lambda,weight] = quadpts(3);
NT = size(elem,1);  
elemRes = zeros(NT,1);
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    fxy = pde.f(pxy);  
    elemRes = elemRes + weight(p)*fxy.^2;
end
elemRes = diameter.^2.*area.*elemRes;


%% elementwise edge jumps
% edge vectors 
e = node(edge(:,2),:)-node(edge(:,1),:);
% scaled norm vectors he*ne
ne = [-e(:,2),e(:,1)]; % stored in rows
% information of left and right elements
k1 = edge2elem(:,1); k2 = edge2elem(:,2); 
% grad of uh 
Dlambda = gradbasis(node,elem);
gradLu = gradu(node,elem(k1,:),u,Dlambda(k1,:,:));
gradRu = gradu(node,elem(k2,:),u,Dlambda(k2,:,:));
% jump of gradu
Jumpu = gradLu-gradRu;
Jumpu(k1==k2,:) = gradLu(k1==k2,:);
% edgeJump
edgeJump = dot(Jumpu',ne').^2; edgeJump = edgeJump';
% elemJump
elemJump = sum(edgeJump(elem2edge),2);

%% Local error indicator
eta = (abs(elemRes) + elemJump).^(1/2);