clc;clear;close all

Domain = @Rectangle_Domain;
h0 = 0.2;
[fd,fh,BdBox,pfix] = feval(Domain);
[node,elem] = distmesh2d(fd,fh,h0,BdBox,pfix);

figure;
showmesh(node,elem);
plotcircle(0,0,1)

phi = @(p) dcircle(p,0,0,1);

deps = 1e-8;



%% 1. 去除所有顶点在外侧的单元

d = phi(node);
ip = (abs(d)<0.1*h0);
node = back2boundary(node,ip,phi);

d = phi(node);  Ix = (d>=0);
elemIf = Ix(elem);
idOut = (sum(elemIf,2)==3);
elem = elem(~idOut,:);

figure;
showmesh(node,elem);
plotcircle(0,0,1)

%% 2. 将很近和仅有一个顶点在外侧的三角形的该顶点拉回到边界
ip1 = 1;
while ~isempty(ip1)
    
    d = phi(node);
    ip = (abs(d)<0.1*h0);
    node = back2boundary(node,ip,phi);
    
    d = phi(node);  Ix = (d>=0);  % 去除恰好一个点在边界上的外部
    elemIf = Ix(elem);
    ix2 = (sum(elemIf,2)==3);
    elem = elem(~ix2,:); elemIf2 = elemIf(ix2,:);
    
    d = phi(node);  Ix = (d>0); % 将仅有一个外部点在外面的投影回边界
    elemIf = Ix(elem);
    ix1 = (sum(elemIf,2)==1);
    Elem1 = elem(ix1,:); elemIf1 = elemIf(ix1,:);
    ip1 = Elem1(elemIf1);
    node = back2boundary(node,ip1,phi);    
end

%% 2. Reorder the vertices
node = node(unique(elem(:)),:);
NT = size(elem,1);
if ~iscell(elem), elem = mat2cell(elem,ones(NT,1),length(elem(1,:))); end
[~,~,totalid] = unique(horzcat(elem{:})');
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';
elem = cell2mat(elem);

figure;
showmesh(node,elem);

[node,elem] = distmesh2dTemp(node,phi,fh,h0,BdBox,pfix);

figure;
showmesh(node,elem);
% findnode(node);
% findelem(node,elem);
% plotcircle(0,0,1)











