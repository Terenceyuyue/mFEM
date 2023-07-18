clc; clear; close all;

%% Parameters
maxIt = 5;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an intitial mesh
[node,elem] = squaremesh([0 1 0 1],0.25);
bdNeumann = 'abs(y-0)<1e-4 | abs(x-1)<1e-4';

%% Get the data of the pde
pde = Poissondata;

%% Finite Element Method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine(node,elem);
    % set boundary
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve the equation
    [uh,info] = PoissonCR(node,elem,pde,bdStruct);
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1/(sqrt(size(node,1))-1);
    if N(k) < 2e3  % show mesh and solution for small size
        figure(1);
        [uhI,nodeI,elemI] = interpolantCR(uh,node,elem);
        showresult(nodeI,elemI,pde.uexact,uhI);
        pause(1);
    end
    % compute error
    ErrL2(k) = getL2error(node,elem,pde.uexact,uh);
    ErrH1(k) = getH1error(node,elem,pde.Du,uh);
end

%% Plot convergence rates and display error table
figure(2);
showrateh(h,ErrH1,ErrL2);
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed.