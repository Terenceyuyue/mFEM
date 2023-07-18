clc; clear; close all;
%% Parameters
maxIt = 5;
h = zeros(maxIt,1); NNdof = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an initial mesh
g = [2  2  2  2  2  2   % decomposed geometry matrix
    0  1  1 -1 -1  0
    1  1 -1 -1  0  0
    0  0  1  1 -1 -1
    0  1  1 -1 -1  0
    1  1  1  1  1  1
    0  0  0  0  0  0];
[p,e,t] = initmesh(g,'hmax',0.5); % initial mesh
node = p';  elem = t(1:3,:)';
bdNeumann = []; % only Dirichlet condition for elasticityNavier

%% Get the PDE data
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata1(para);

%% Finite element method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % set boundary
    bdStruct = setboundary(node,elem,bdNeumann);
    % set up solver type
    option.solver = 'mg';
    option.J = k+1;
    % solve the equation
    uh = elasticityNavier(node,elem,pde,bdStruct,option);
    uh = reshape(uh,[],2);
    % record and plot
    NNdof(k) = length(uh);
    h(k) = 1/(sqrt(size(node,1))-1);
    if NNdof(k)<2e3
        figure(1);
        showresult(node,elem,pde.uexact,uh(:,1));
        pause(1);
    end
    % compute error
    tru = eye(2); trDu = eye(4);
    errL2 = zeros(1,2);  errH1 = zeros(1,2); % square
    for id = 1:2
        uid = uh(:,id);
        u = @(pz) pde.uexact(pz)*tru(:, id);
        Du = @(pz) pde.Du(pz)*trDu(:, 2*id-1:2*id);
        errL2(id) = getL2error(node,elem,u,uid);
        errH1(id) = getH1error(node,elem,Du,uid);
    end
    ErrL2(k) = norm(errL2);
    ErrH1(k) = norm(errH1);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,NNdof,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order), L2-norm
% (2nd order) is observed.