function tri = boundarymesh3D(elem)
% Triangulation
totalFace = [elem(:,[1,2,3]);elem(:,[1,2,4]);elem(:,[1,3,4]);elem(:,[2,3,4])];
totalFace = sort(totalFace,2);
[~,i1,totalJ] = unique(totalFace,'rows');
id = 1:max(totalJ);
id_rep = histcounts(totalJ,id);  %  num --> id
s = i1(id_rep==1);
tri = totalFace(s,:);

