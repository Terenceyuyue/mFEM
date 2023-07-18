function showmesh(node,elem,options)
%Showmesh displays a mesh in 2D and 3D.
%
% Copyright (C) Terence Yu.

if nargin==2, options.FaceAlpha = 0.4; end

dim = size(node,2);

% ------------------- Triangulation --------------------------
% 2D
if ~iscell(elem) && dim==2
    h = patch('Faces', elem, 'Vertices', node);
end

% 3D 
isexpr = isfield(options,'expr'); % part of the 3-D mesh

if ~iscell(elem) && dim==3 && ~isexpr
    elemcut = elem;
    totalFace = [elemcut(:,[1,2,3]); elemcut(:,[1,2,4]); elemcut(:,[1,3,4]); elemcut(:,[2,3,4])];
    totalFace = sort(totalFace,2);
    face = unique(totalFace,'rows');  % triangles
    h = patch('Faces', face, 'Vertices', node);
end

% Part of the 3-D mesh, e.g. cavity (cylinder with hole): only need to plot surface
if ~iscell(elem) && dim==3 && isexpr      
    expr = options.expr;
    x = node(:,1);  y = node(:,2);  z = node(:,3); %#ok<NASGU>
    incl = find(eval(expr));
    elemcut = elem(any(ismember(elem,incl),2),:); % cutted cavity 
    
    surface = boundarymesh3D(elem);  % outer and inner surfaces of cavity
    surface0 = surface(any(ismember(surface,incl),2),:); % cutted surface
    surface1 = boundarymesh3D(elemcut); % surface of cutted cavity (consisting of cutted surface and the section)
    surface_section = setdiff(surface1, surface0, 'rows' ); % surface of the section
    
    h = patch('Faces', surface_section, 'Vertices', node);
    set(h,'facecolor',[0.9 0.9 0.9],'edgecolor','k');
    h = patch('Faces', surface0, 'Vertices', node);
end


% ------------------- Polygonal/polyhedron mesh -----------------------
if iscell(elem)
    if iscell(elem{1}), elem = vertcat(elem{:}); end
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    face = vertcat(tpad{:});  % polygon
    h = patch('Faces', face, 'Vertices', node);
end

if dim==3
    view(3); set(h,'FaceAlpha',options.FaceAlpha); % transparency
end

facecolor = [0.5 0.9 0.45];
if isfield(options,'facecolor')
    facecolor = options.facecolor;
end
set(h,'facecolor',facecolor,'edgecolor','k');
axis equal; axis tight;


