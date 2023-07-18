function findelem(node,elem,varargin)
%Findelem highlights some elements
%
% Copyright (C) Terence Yu

hold on

NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

range = 1:NT;
if nargin==3, range = unique(varargin{1}); end
range = range(:);

center = zeros(length(range),2);
s = 1;
for iel = range(:)' % only valid for row vector
    index = elem{iel};
    V = node(index, :); 
    center(s,:) = polycentroid(V);
    s = s+1;
end
plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
    'MarkerFaceColor','y','MarkerSize',18);
text(center(:,1)-0.01,center(:,2),int2str(range),'FontSize',12,...
    'FontWeight','bold','Color','k');

hold off

end

function centroid = polycentroid(V)
% polycentroid returns the x,y coordinates of centroid of polygon in 2-D
% The vertices must be in cyclic order
%
%    V = [x,y]
%
% Copyright (C) Terence Yu.

V1 = circshift(V,-1);
area_components = V(:,1).*V1(:,2) - V1(:,1).*V(:,2);
ar = 0.5*(sum(area_components)); 
centroid = sum((V+V1).*repmat(area_components,1,2))/(6*ar);
end