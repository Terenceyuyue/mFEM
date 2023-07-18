function [uhI,nodeI,elemI] = interpolantCR(uh,node,elem)

% elem2edge
allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(allEdge,2);
[~, ~, totalJ] = unique(totalEdge,'rows');
NT = size(elem,1);
elem2edge = reshape(totalJ,NT,3);

% uhI,nodeI,elemI
elemI = reshape(1:3*NT,NT,3);
nodeI = node(elem(:),:);
uhI = zeros(3*NT,1);
uhI(1:NT) = -uh(elem2edge(:,1)) + uh(elem2edge(:,2)) + uh(elem2edge(:,3));
uhI(NT+(1:NT)) = uh(elem2edge(:,1)) - uh(elem2edge(:,2)) + uh(elem2edge(:,3));
uhI(2*NT+(1:NT)) = uh(elem2edge(:,1)) + uh(elem2edge(:,2)) - uh(elem2edge(:,3));
