function [node,elem] = bisect(node,elem,elemMarked)
%bisect refines the elements using newest-node bisection
%
% Copyright (C) Long Chen, modified by Terence Yu.

%% Construct auxiliary data structure
aux = auxstructure(node,elem);
elem2edge = aux.elem2edge; edge = aux.edge; neighbor = aux.neighbor;
clear aux;
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% (peak) base set 
isCutEdge = false(NE,1); 
while sum(elemMarked)>0
    base = elem2edge(elemMarked,1);   isCutEdge(base) = true; 
    refineNeighbor = neighbor(elemMarked,1);
    baseNeighbor = elem2edge(refineNeighbor,1);
    elemMarked = refineNeighbor(~isCutEdge(baseNeighbor));
end

%% Add new nodes
Nbase = sum(isCutEdge);
edgeCutNumber = N + (1:Nbase);
edgebase = edge(isCutEdge,:);
node(edgeCutNumber,:) = (node(edgebase(:,1),:) + node(edgebase(:,2),:))/2;

%% Refine elements
edgeCutNumber = zeros(NE,1);
edgeCutNumber(isCutEdge) = (N+1:N+Nbase)';
for k = 1:2
    t = find(isCutEdge(elem2edge(:,1))==true);
    newNT = length(t);
    if newNT==0, break; end
    L = t;  R = NT+(1:newNT);
    p1 = elem(t,1); p2 = elem(t,2); p3 = elem(t,3);
    p4 = edgeCutNumber(elem2edge(t,1));
    elem(L,:) = [p4, p1, p2];  
    elem(R,:) = [p4, p3, p1];
    elem2edge(L,1) = elem2edge(t,3);
    elem2edge(R,1) = elem2edge(t,2);    
    NT = NT + newNT;
end
