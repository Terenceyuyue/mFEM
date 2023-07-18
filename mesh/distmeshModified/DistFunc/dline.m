
function d = dline(p,x1,y1,x2,y2)
% By convention, a point located at the left hand side of the line
% is inside the region and it is assigned a negative distance value.
a = [x2-x1,y2-y1]; a = a/norm(a);
b = [p(:,1)-x1,p(:,2)-y1];
d = b(:,1)*a(2) - b(:,2)*a(1);

