function d=dcylinder(p,xc,yc,r,z1,z2)

% (xc,yc,r): curved surface
% z1,z2: bottom and top

x = p(:,1); y = p(:,2); z = p(:,3);

d1 = sqrt((x-xc).^2+(y-yc).^2) - r;
d2 = -z-z1; % bottom
d3 = z-z2; % top

d = max(max(d1,d2),d3);
