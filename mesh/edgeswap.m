function elem = edgeswap(node,elem,signSwap)
%EDGESWAP swaps the diagonal of a quadrilateral formed by two triangles
% signSwap = 1: swap edges with slope>0
% signSwap = -1: swap edges with slope<0
%
% Copyright (C) Long Chen, Modified by Terence Yu. 

if size(elem,2)~=3
    error('Error: not a triangulation');
end
if nargin==2, signSwap = 1; end % default: slope>0

%% Construct necessary data structure
auxT = auxstructure(node,elem);
edge = auxT.edge;
edge2elem = auxT.edge2elem;
clear auxT;

%% Find edges to be swapped
x1 = node(edge(:,1),1);  y1 = node(edge(:,1),2);  
x2 = node(edge(:,2),1);  y2 = node(edge(:,2),2);
% exclude horizontal and vertical edges
idEdgeSwap = false(size(edge,1),1);
eps = 1e-8;
isEdgeNeed = (abs(x1-x2)>eps) & (abs(y1-y2)>eps); % x1==x2, y1==y2
slope = (y2(isEdgeNeed)-y1(isEdgeNeed))./(x2(isEdgeNeed)-x1(isEdgeNeed));
idEdgeSwap(isEdgeNeed) = (sign(slope)==signSwap);

%% Swap edges
if any(idEdgeSwap)
    %   4 - 3            4 - 3
    %   | / |   ---->    | \ | 
    %   1 - 2            1 - 2
    % left and right triangles and local indices
    k1 = edge2elem(idEdgeSwap,1);  k2 = edge2elem(idEdgeSwap,2);
    i1 = edge2elem(idEdgeSwap,3);  i2 = edge2elem(idEdgeSwap,4);
    % the index of p1,p2,p3,p4 
    next = [2;3;1];  NT = size(elem,1);
    c4 = k1 + (i1-1)*NT;  % extract entries of A(:)
    c1 = k1 + (next(i1)-1)*NT; 
    c2 = k2 + (i2-1)*NT;
    c3 = k2 + (next(i2)-1)*NT;    
    p1 = elem(c1);  p2 = elem(c2); p3 = elem(c3); p4 = elem(c4);  
    % new triangles
    leftElem = [p1, p2, p4];  
    rightElem = [p3, p4, p2];
    elem(k1,:) = leftElem;
    elem(k2,:) = rightElem;
end