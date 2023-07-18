function t = my_delaunay(p,varargin)
% signSwap = varargin{1}
% signSwap = 1: swap edges with slope>0
% signSwap = -1: swap edges with slope<0 (default)
% signSwap ~= -1,1: don't swap edges

% determine the supertriangle
xmin = min(p(:,1)); xmax = max(p(:,1));
ymin = min(p(:,2)); ymax = max(p(:,2));
hx = xmax-xmin; hy = ymax-ymin;

eps = 0.05*hy;
pts = [
    xmin+hx/2,       ymax+hy;
    xmin-hx/2-2*eps, ymin-eps;
    xmax+hx/2+2*eps, ymin-eps
    ];
node = [pts; p];   elem = [1 2 3];

for k = 1:size(p,1)
    % find effected triangles assocaited with k-th point
    z1 = node(elem(:,1),:);
    z2 = node(elem(:,2),:);
    z3 = node(elem(:,3),:);
    id_tri = inCircle(p(k,:),z1,z2,z3);
    tri = elem(id_tri,:); % effected triangles
    edge_buffer = [tri(:,[2,3]); tri(:,[3,1]); tri(:,[1,2])];
    
    % obtain edges of the enclosing polygon
    aa = full(sparse(edge_buffer(:,1),edge_buffer(:,2),1));
    aa(aa+aa'>=2) = 0;
    [i,j] = find(aa);
    edge_poly = [i,j];
    
    % update elem
    N = 3+k; n_poly = size(edge_poly,1);
    tri_poly = [edge_poly, N*ones(n_poly,1)]; % new triangles
    elem = [elem(~id_tri,:); tri_poly];
end

% remove any triangles that use the supertriangle vertices
id = (elem(:,1)>3) & (elem(:,2)>3) & (elem(:,3)>3);
t = elem(id,:)-3;

% swap the diagonal of a quadrilateral formed by two triangles
signSwap = -1; % swap edges with slope<0
if nargin==2, signSwap = varargin{1}; end
t = edgeswap(p,t,signSwap);
end % end of my_delaunay

function id = inCircle(p0,p1,p2,p3)

% detect if point p0 lies in the circumcircle of triangle with vertices
% p1,p2,p3 stored in a counterclockwise order
% see Wikipedia: https://en.wikipedia.org/wiki/Delaunay_triangulation
 
if size(p0,1)==1, p0 = repmat(p0,size(p1,1),1); end

x0 = p0(:,1); y0 = p0(:,2); x1 = p1(:,1); y1 = p1(:,2);
x2 = p2(:,1); y2 = p2(:,2); x3 = p3(:,1); y3 = p3(:,2);
detA = - x0.^2.*x1.*y2 + x0.^2.*x1.*y3 + x0.^2.*x2.*y1 - x0.^2.*x2.*y3 - x0.^2.*x3.*y1 + ...
    x0.^2.*x3.*y2 + x0.*x1.^2.*y2 - x0.*x1.^2.*y3 - x0.*x2.^2.*y1 + x0.*x2.^2.*y3 + ...
    x0.*x3.^2.*y1 - x0.*x3.^2.*y2 + x0.*y1.^2.*y2 - x0.*y1.^2.*y3 - x0.*y1.*y2.^2 + ...
    x0.*y1.*y3.^2 + x0.*y2.^2.*y3 - x0.*y2.*y3.^2 - x1.^2.*x2.*y0 + x1.^2.*x2.*y3 + ...
    x1.^2.*x3.*y0 - x1.^2.*x3.*y2 + x1.*x2.^2.*y0 - x1.*x2.^2.*y3 - x1.*x3.^2.*y0 + ...
    x1.*x3.^2.*y2 - x1.*y0.^2.*y2 + x1.*y0.^2.*y3 + x1.*y0.*y2.^2 - x1.*y0.*y3.^2 - ...
    x1.*y2.^2.*y3 + x1.*y2.*y3.^2 - x2.^2.*x3.*y0 + x2.^2.*x3.*y1 + x2.*x3.^2.*y0 - ...
    x2.*x3.^2.*y1 + x2.*y0.^2.*y1 - x2.*y0.^2.*y3 - x2.*y0.*y1.^2 + x2.*y0.*y3.^2 + ...
    x2.*y1.^2.*y3 - x2.*y1.*y3.^2 - x3.*y0.^2.*y1 + x3.*y0.^2.*y2 + x3.*y0.*y1.^2 - ...
    x3.*y0.*y2.^2 - x3.*y1.^2.*y2 + x3.*y1.*y2.^2;

itol = 1e-6;
id = detA>itol;
end

