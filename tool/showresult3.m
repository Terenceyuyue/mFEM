function showresult3(node,elem,u,uh,options)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

if nargin==4, options.FaceAlpha = 1; end

clf;  % clear figures
set(gcf,'Units','normal');
set(gcf,'Position',[0.25,0.25,0.7,0.25]);

%% Find faces at z = zmin
allFace = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
totalFace = sort(allFace,2);
face = unique(totalFace,'rows'); NF = size(face,1);
zmin = min(node(:,3)); 
nodeFace = (node(face(:,1),:) + node(face(:,2),:) + node(face(:,3),:))/3;
z = nodeFace(:,3); %#ok<NASGU>
Idx = false(NF,1);
expr = sprintf('z==%f',zmin);
for i = 1:NF
    id = eval(expr);
    Idx(id) = true;
end
facemin = face(Idx,:);

%% Plot mesh
subplot(1,3,1);
showmesh(node,elem,options);

%% Plot numerical solution
subplot(1,3,2);
showsolution(node(:,1:2),facemin,uh(1:size(node,1),1));
title('Display numerical solution at z = z_{min}');

%% Plot exact solution
subplot(1,3,3);
ue = u(node); ue = ue(:,1);
showsolution(node(:,1:2),facemin,ue(1:size(node,1)));
title('Display exact solution at z = z_{min}');