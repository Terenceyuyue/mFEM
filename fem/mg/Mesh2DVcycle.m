function [node,elem,Pro,Res] = Mesh2DVcycle(node,elem,J)

if J<=1
    Pro = []; Res = []; return;
end

Pro = cell(J-1,1); Res = cell(J-1,1);
for j = 2:J    
    % auxiliary mesh data
    aux = auxstructure(node,elem);
    edge = aux.edge; elem2edge = aux.elem2edge;
    N = size(node,1); NT = size(elem,1); 
    NE = size(edge,1); Ndof = 3;
    
    % HB in the current level 
    HB = zeros(NE,3); 
    HB(:,1) = (1:NE)'+N; HB(:,2:3) = edge;
    
    % Prolongation matrix
    Coarse2Fine = (1:N)'; CoarseId = (1:N)';
    nCoarse = N; nf = NE; nFine = nCoarse+nf;
    ii = [Coarse2Fine; HB(:,1); HB(:,1)];
    jj = [CoarseId; HB(:,2); HB(:,3)];
    ss = [ones(nCoarse,1); 0.5*ones(nf,1); 0.5*ones(nf,1)];
    Pro{j-1} = sparse(ii,jj,ss,nFine,nCoarse);
    Res{j-1} = Pro{j-1}';
    
    % add new nodes: middle points of all edges
    node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;
    
    % add new elements: refine each triangle into four triangles as follows
    % 3
    % | \
    % 5- 4
    % |\ |\
    % 1- 6- 2
    t = 1:NT; p = zeros(NT,2*Ndof);
    p(:,1:3) = elem;
    p(:,4:6) = elem2edge + N;
    elem(t,:) = [p(t,1), p(t,6), p(t,5)];
    elem(NT+1:2*NT,:) = [p(t,6), p(t,2), p(t,4)];
    elem(2*NT+1:3*NT,:) = [p(t,5), p(t,4), p(t,3)];
    elem(3*NT+1:4*NT,:) = [p(t,4), p(t,5), p(t,6)];
end


