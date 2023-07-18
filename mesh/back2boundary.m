function node = back2boundary(node,id,fd)

d = fd(node);
eps = 1e-8;
n1 = (fd([node(id,1)+eps,node(id,2)])-fd(node(id,:)))/eps;
n2 = (fd([node(id,1),node(id,2)+eps])-fd(node(id,:)))/eps;
nn = (n1.^2 + n2.^2); % better performance without sqrt
nn(nn<eps) = 1;
n1 = n1./nn; n2 = n2./nn;
node(id,:) = node(id,:)-[d(id).*n1,d(id).*n2];   % Project back to the boundary