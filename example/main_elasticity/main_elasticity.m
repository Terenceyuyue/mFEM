clc; clear; close all;
%% Parameters
maxIt = 5;
h = zeros(maxIt,1); NNdof = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an initial mesh
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 4; Ny = 4; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);
bdNeumann = 'abs(y-0)<1e-4 | abs(x-1)<1e-4'; % string for Neumann

%% Get the PDE data
lambda = 1; mu = 1;
para.lambda = lambda; para.mu = mu;
pde = elasticitydata(para);

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
    uh = elasticity(node,elem,pde,bdStruct,option);
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
