function bdStruct= setboundary1(node,elem1,varargin)
% varargin: string for Neumann boundary

% edge number
bdEdge = find(accumarray(elem1(:),1) == 1);
nE = size(bdEdge,1); 

% initial as Dirichlet (true for Dirichlet, false for Neumann)
bdNodeIdxD = true(nE,1);  nodeEdge = node(bdEdge,:);
if size(node,2)==1, x = nodeEdge; end
if size(node,2)==2, x = nodeEdge(:,1); y = nodeEdge(:,2); end
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1})) 
    for i = 1:nvar
        bdNeumann = varargin{1};
        id = eval(bdNeumann);
        bdNodeIdxD(id) = false;
    end
end
bdStruct.Dirichlet = bdEdge(bdNodeIdxD); % Dirichlet boundary nodes
bdStruct.Neumann = bdEdge(~bdNodeIdxD); % Dirichlet boundary nodes