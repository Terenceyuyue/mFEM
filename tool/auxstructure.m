function aux = auxstructure(node,elem)
%auxstructure: auxiliary data structure in 2-D
%
% Copyright (C) Long Chen, modified by Terence Yu.

dim = size(elem,2); NT = size(elem,1);

if dim == 4 % rectangle element
    elem = mat2cell(elem,ones(size(elem,1),1), 4); 
end

%% elem2edge
if ~iscell(elem) % Triangulation    
    totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
    [~, i1, totalJ] = unique(totalEdge,'rows');    
    elem2edge = reshape(totalJ,NT,3);
else   % rectangle
    % totalEdge
    shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
    T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
    v0 = horzcat(elem{:})'; % the starting points of edges
    v1 = horzcat(T1{:})'; % the ending points of edges
    totalEdge = sort([v0,v1],2);    
    % elem2edge
    [~, i1, totalJ] = unique(totalEdge,'rows');
    elemLen = cellfun('length',elem); % length of each elem
    elem2edge = mat2cell(totalJ',1,elemLen)';
end

%% edge, bdEdge
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
bdEdge = edge(s==1,:); % not counterclockwise

%% edge2elem
if ~iscell(elem)
    totalJelem = repmat((1:NT)',3,1);
else
    Num = num2cell((1:NT)');    Len = num2cell(elemLen);
    totalJelem = cellfun(@(n1,n2) n1*ones(n2,1), Num, Len, 'UniformOutput', false);
    totalJelem = vertcat(totalJelem{:});   
end
[~, i2] = unique(totalJ(end:-1:1),'rows');
i2 = length(totalEdge)+1-i2;
edge2elem = totalJelem([i1,i2]);
if ~iscell(elem) % triangulation 
    edge2elem = [edge2elem, ceil(i1/NT), ceil(i2/NT)];
end

%% neighbor
neighbor = zeros(NT,dim);
for iel = 1:NT
    indexE = elem2edge(iel,:);  
    ia = edge2elem(indexE,1); ib = edge2elem(indexE,2);
    ia(ia==iel) = ib(ia==iel);
    neighbor(iel,:) = ia';
end


aux.node = node; aux.elem = elem;
aux.elem2edge = elem2edge;
aux.edge = edge; aux.bdEdge = bdEdge;
aux.edge2elem = edge2elem;
aux.neighbor = neighbor;