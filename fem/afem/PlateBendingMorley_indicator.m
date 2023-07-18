function eta = PlateBendingMorley_indicator(node,elem,wh,pde)
% This function returns the local error indicator of Plate bending problems
% with Morley element applied.
%
% Copyright (C) Terence Yu.

%% preparation for the computation
% auxiliary data
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge, i1, totalJ] = unique(totalEdge,'rows');
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
Dlambdax = [grad1(:,1), grad2(:,1), grad3(:,1)];
Dlambday = [grad1(:,2), grad2(:,2), grad3(:,2)];
% elem2dof
N = size(node,1); Ndof = 6;
elem2dof = [elem, elem2edge+N];
% sgnelem, sgnbase
[~, i2] = unique(totalEdge(end:-1:1,:),'rows');  
i2 = size(totalEdge,1)+1-i2;            % last occurrence
bdEdgeIdx = (i1==i2);
sgnelem = sign([elem(:,3)-elem(:,2), ...
                elem(:,1)-elem(:,3), ...
                elem(:,2)-elem(:,1)]);
E = false(NE,1); E(bdEdgeIdx) = 1; sgnbd = E(elem2edge);
sgnelem(sgnbd) = 1;
sgnbase = ones(NT,Ndof); sgnbase(:,4:6) = sgnelem;

%% elementwise residuals --- the first term
% 2-D quadrature
[lambda,weight] = quadpts(4);
elemRes = zeros(NT,1);
for p = 1:length(weight)
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*z1+lambda(p,2)*z2+lambda(p,3)*z3;
    fxy = pde.f(pxy);  
    elemRes = elemRes + weight(p)*fxy.^2;
end
elemRes = diameter.^4.*area.*elemRes;

%% connectivity list of jump integral
[~,weight1d] = quadpts1(4); ng = length(weight1d);
elem2dofP = zeros(NT,3*ng); % P: plus +
%elem2dofM = zeros(NT,3*ng); % M: minus -
for i = 1:3  % i-th side
    ei = elem2edge(:,i); sgni = sgnelem(:,i);
    idP = repmat(1:ng, NT,1);  
    idP(sgni<0,:) = repmat((ng:-1:1)+NE*ng, sum(sgni<0), 1);
    elem2dofP(:,(1:ng)+(i-1)*ng) = idP + repmat((ei-1)*ng,1,ng);
end
index = (elem2dofP>ng*NE);
elem2dofM = elem2dofP + (-ng*NE)*index + ng*NE*(~index);

%% elementwise jump of wh --- the second term
% boundary quadrature points and weights on the reference triangle
[lambda1d,weight1d] = quadpts1(4);  ng = length(weight1d);
[~,id] = sort(lambda1d(:,1)); 
lambda1d = lambda1d(id,:); weight1d = weight1d(id);
lambdae1 = [zeros(ng,1), lambda1d(:,2), lambda1d(:,1)];
lambdae2 = [lambda1d(:,1), zeros(ng,1), lambda1d(:,2)];
lambdae3 = 1 - lambdae1 - lambdae2;
lambdaRef = [lambdae1; lambdae2; lambdae3];
% preparation for the computation of basis functions
% xi, eta
xi = [e1(:,1), e2(:,1), e3(:,1)]; 
eta = [e1(:,2), e2(:,2), e3(:,2)]; 
% sij
sb = zeros(NT,3,3);
for i = 1:3
    for j = 1:3
        sb(:,i,j) = -xi(:,i).*xi(:,j) - eta(:,i).*eta(:,j);
    end
end
% coefficients in the basis functions
L = he(elem2edge);
ind = [1 2 3; 2 3 1; 3 1 2];  % rotation index
it = ind(:,1); jt = ind(:,2); kt = ind(:,3);
c0 = zeros(NT,3); c1 = c0; c2 = c0;
for i = 1:3
    j = ind(i,2); k = ind(i,3);
    c0(:,i) = 1./(2*area.^2);
    c1(:,i) = sb(:,i,j)./L(:,j).^2;
    c2(:,i) = sb(:,k,i)./L(:,k).^2;
end
c3 = c1+c2;
% interpolant of wh at boundary quadrature points
elemwhP = zeros(NT,3*ng);
for p = 1:3*ng
    % basis functions at the p-th quadrture point
    base = zeros(NT,Ndof);
    % i = 1,2,3    
    base(:,it) = repmat(lambdaRef(p,it).^2,NT,1) ...
               + c3.*repmat(lambdaRef(p,jt).*lambdaRef(p,kt),NT,1) ...
               + c2.*repmat(lambdaRef(p,kt).*lambdaRef(p,it),NT,1) ...
               + c1.*repmat(lambdaRef(p,it).*lambdaRef(p,jt),NT,1);
    % i = 4,5,6
    ci = 2*repmat(area,1,3)./L(:,it);
    base(:,3+it) = ci.*repmat(lambdaRef(p,it).*(lambdaRef(p,it)-1),NT,1);
    % equipped with global sign  
    base = sgnbase.*base;     
    % interpolant at the p-th quadrture point
    for i = 1:Ndof
        elemwhP(:,p) = elemwhP(:,p) + wh(elem2dof(:,i)).*base(:,i);
    end
end
whI = zeros(2*NE*ng,1);   whI(elem2dofP) = elemwhP;
elemwhM = whI(elem2dofM);
% elementwise jump of wh
elem2Jumpval = (elemwhP - elemwhM).^2;
ce = 0.5*ones(NE,1); ce(bdEdgeIdx) = 1;
elemJumpval = zeros(NT,1);
for i = 1:3
    hei = he(elem2edge(:,i));  
    cei = ce(elem2edge(:,i)).*hei.^(-3);
    Jump = elem2Jumpval(:,(1:ng)+(i-1)*ng);
    elemJumpval = elemJumpval + cei.*hei.*(Jump*weight1d(:));
end

%% elementwise jump of whx, why --- the third term
% interpolant of whx, why at boundary quadrature points
elemwhxP = zeros(NT,3*ng); 
elemwhyP = zeros(NT,3*ng);
for p = 1:3*ng
    % derivatives of basis functions at the p-th quadrture point
    basex = zeros(NT,Ndof); basey = zeros(NT,Ndof);
    % i = 1,2,3
    lambdaRefpit = repmat(lambdaRef(p,it),NT,1);
    lambdaRefpjt = repmat(lambdaRef(p,jt),NT,1);
    lambdaRefpkt = repmat(lambdaRef(p,kt),NT,1);
    basex(:,it) = 2*lambdaRefpit.*Dlambdax(:,it) ...
        + c3.*(Dlambdax(:,jt).*lambdaRefpkt+lambdaRefpjt.*Dlambdax(:,kt)) ...
        + c2.*(Dlambdax(:,kt).*lambdaRefpit+lambdaRefpkt.*Dlambdax(:,it)) ...
        + c1.*(Dlambdax(:,it).*lambdaRefpjt+lambdaRefpit.*Dlambdax(:,jt));
    basey(:,it) = 2*lambdaRefpit.*Dlambday(:,it) ...
        + c3.*(Dlambday(:,jt).*lambdaRefpkt+lambdaRefpjt.*Dlambday(:,kt)) ...
        + c2.*(Dlambday(:,kt).*lambdaRefpit+lambdaRefpkt.*Dlambday(:,it)) ...
        + c1.*(Dlambday(:,it).*lambdaRefpjt+lambdaRefpit.*Dlambday(:,jt));
    % i = 4,5,6
    ci = 2*repmat(area,1,3)./L(:,it);
    basex(:,3+it) = ci.*(Dlambdax(:,it).*(lambdaRefpit-1)+lambdaRefpit.*Dlambdax(:,it));
    basey(:,3+it) = ci.*(Dlambday(:,it).*(lambdaRefpit-1)+lambdaRefpit.*Dlambday(:,it));
    % equipped with global sign 
    basex = sgnbase.*basex;  basey = sgnbase.*basey;
    % interpolant at the p-th quadrture point
    for i = 1:Ndof
        elemwhxP(:,p) = elemwhxP(:,p) + wh(elem2dof(:,i)).*basex(:,i);
        elemwhyP(:,p) = elemwhyP(:,p) + wh(elem2dof(:,i)).*basey(:,i);
    end
end
whxI = zeros(2*NE*ng,1);      whyI = zeros(2*NE*ng,1);
whxI(elem2dofP) = elemwhxP;   whyI(elem2dofP) = elemwhyP;
elemwhxM = whxI(elem2dofM);   elemwhyM = whyI(elem2dofM);
% elementwise jump of whx, why
elem2Jumpx = elemwhxP - elemwhxM;
elem2Jumpy = elemwhyP - elemwhyM;
elemJumpgrad = zeros(NT,1);
for i = 1:3
    hei = he(elem2edge(:,i));
    cei = ce(elem2edge(:,i)).*hei.^(-1);
    neix = ne(elem2edge(:,i),1); neiy = ne(elem2edge(:,i),2);
    Jumpnx = elem2Jumpx(:,(1:ng)+(i-1)*ng).*repmat(neix,1,ng);
    Jumpny = elem2Jumpy(:,(1:ng)+(i-1)*ng).*repmat(neiy,1,ng);
    Jumpn = (Jumpnx+Jumpny).^2;    
    elemJumpgrad = elemJumpgrad + cei.*hei.*(Jumpn*weight1d(:));
end

%% Local error indicator
elemJump = elemJumpval + elemJumpgrad;
eta = (abs(elemRes) + abs(elemJump)).^(1/2);