function [node,elem,Pro,Res] = Mesh1DVcycle(node,elem,J)

if J<=1
    Pro = []; Res = []; return;
end

Pro = cell(J-1,1); Res = cell(J-1,1);
for j = 2:J
    
    N = size(node,1); nel = size(elem,1);
    % HB in the current level    
    HB = zeros(nel,3);
    HB(:,1) = (1:nel)'+N;  % number of the new nodes
    HB(:,2:3) = elem;     % number of the left and right nodes
    
    % Prolongation matrix
    Coarse2Fine = (1:N)'; CoarseId = (1:N)';
    nCoarse = N; nf = nel; nFine = nCoarse+nf;
    ii = [Coarse2Fine; HB(:,1); HB(:,1)];
    jj = [CoarseId; HB(:,2); HB(:,3)];
    ss = [ones(nCoarse,1); 0.5*ones(nf,1); 0.5*ones(nf,1)];
    Pro{j-1} = sparse(ii,jj,ss,nFine,nCoarse);
    Res{j-1} = Pro{j-1}';
    
    % Refine the mesh
    [node,elem] = uniformrefine1(node,elem);
end