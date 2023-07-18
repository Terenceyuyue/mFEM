function [Pro,Res] = uniformtransferoperator1(elem,J)
%uniformtransferoperator1 returns the prolongation matrices and restriction
% matrices of meshes obtained by using uniformrefine1.m
% elem: the final mesh
% J: level number
% 
% Copyright (C) Terence Yu.

Pro = cell(J-1,1); Res = cell(J-1,1);

% 1----3----2

for j = J:-1:2
    
    %% HB in the level J-1 
    N = max(elem(:));  nel = size(elem,1);
    nel = nel/2;  % number of coarse mesh
    t1 = 1:nel; t2 = t1+nel;
    HB = zeros(nel,3);
    HB(:,1) = elem(t1, 2);  % number of the new nodes
    elem = [elem(t1, 1),  elem(t2, 2)]; % coarse mesh
    HB(:,2:3) = elem;
    %% Transfer matrices in the level J-1
    Coarse2Fine = unique(elem(:)); CoarseId = Coarse2Fine;
    nCoarse = length(CoarseId);  nFine = N;   nf = size(HB,1);
    ii = [Coarse2Fine; HB(:,1); HB(:,1)];
    jj = [CoarseId; HB(:,2); HB(:,3)];    
    ss = [ones(nCoarse,1); 0.5*ones(nf,1); 0.5*ones(nf,1)];
    Pro{j-1} = sparse(ii,jj,ss,nFine,nCoarse);
    Res{j-1} = Pro{j-1}';
    
end

