function [node,elem] = uniformrefine1(node,elem)
%uniformrefine1 uniformly refines a 1-D line segments.
%
% Copyright (C) Terence Yu.

N = size(node,1); nel = size(elem,1); np = 2;
%% Add new nodes: middle points of all segments
node(N+1:N+nel) = (node(elem(:,1))+node(elem(:,2)))/2;

%% Refine each line segment into two segments as follows
% 1 --- 3 --- 2
t = 1:nel; p = zeros(nel,2*np);
p(:,1:2) = elem;  p(:,3) = (1:nel)' + N;
elem(t,:) = [p(t,1), p(t,3)];  % 1 -- 3
elem(nel+1:2*nel,:) = [p(t,3), p(t,2)]; % 3 -- 2