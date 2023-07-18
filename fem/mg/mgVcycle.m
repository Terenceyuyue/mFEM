% ---------------------------------------------------------------------
%                         MG solver: mgVcycle.m
% ---------------------------------------------------------------------

% Copyright (C) Long Chen, modified by Terence Yu.

function u = mgVcycle(A,b,Pro,Res)

J = length(Pro)+1;
Ai = cell(J,1);  Ai{J} = A; % matrices in subspaces
for j = J-1:-1:1 % note that J-1:-1:1 is empty when J<=1 (solve it directly)
    Ai{j} = Res{j}*Ai{j+1}*Pro{j};
end

tol = 1e-9; 
Err = 10; u = zeros(size(b));
iter = 0; MaxIt = 100;
while Err>tol && iter<=MaxIt
    r = b-A*u;
    e = Vcycle(A,r,Ai,Pro,Res);
    u = u+e;
    Err = norm(e);  iter = iter + 1;
end

end % end of mgVcycle

% ---------------------------------------------------------------------
%                         Vcycle.m 
% ---------------------------------------------------------------------
function e = Vcycle(A,r,Ai,Pro,Res)
% Solve the residual equation Ae = r by multigrid V-cycle method

J = length(Ai); % level length

% If the problem is small enough, solve it directly
if J<=1
    e = A\r;    return;
end

ri = cell(J,1);            % residual in each level
ei = cell(J,1);            % correction in each level
ri{J} = r;

% ---------- Correction in each level ----------
option = 'forward';
for j = J:-1:2
    % % pre-smoothing: one step
    % ei{j} = R{j}*ri{j};
    ei{j} = smoother(Ai{j},ri{j},option);
    % update and restrict residual
    ri{j-1} = Res{j-1}*(ri{j}-Ai{j}*ei{j});
end
ei{1} = Ai{1}\ri{1}; % exact solver in the coarsest level

% -------------- prolongation and correction ------------
option = 'backward';
for j = 2:J
    % prolongation and correction
    ei{j} = ei{j}+Pro{j-1}*ei{j-1};
    % % post-smoothing: one step
    % ei{j} = ei{j} + R{j}'*(ri{j}-Ai{j}*ei{j});
    rij = ri{j}-Ai{j}*ei{j};
    ei{j} = ei{j} + smoother(Ai{j},rij,option);
end
e = ei{J};
end % end of Vcycle1D

% ---------------------------------------------------------------------
%                         smoother.m 
% ---------------------------------------------------------------------
function ei = smoother(Ai,ri,option)
switch option
    case 'forward'
        Ri = tril(Ai);  % Forward Gauss-Seidel   R = D+L
    case 'backward'
        Ri = triu(Ai);  % Backward Gauss-Seidel  R = D+U
end
ei = Ri\ri;
end % end of smoother;

