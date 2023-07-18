clc; clear; close all; 
%% Parameters
maxIt = 5;
h = zeros(maxIt,1); NNdof = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Get the PDE data
pde = Poissondata();

%% Finite element method
for k = 1:maxIt  % level
    % generate an initial mesh
    [node,elem] = squaremesh([0 1 0 1],0.25,0.25);
    % Final mesh and transfer matrices    
    [node,elem,Pro,Res] = Mesh2DVcycle(node,elem,k);
    % set boundary
    bdNeumann = 'abs(x-1)<1e-4';
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve the equation
    uh = PoissonVcycle(node,elem,pde,bdStruct,Pro,Res);
    % record
    NNdof(k) = length(uh);
    h(k) = 1/(sqrt(size(node,1))-1);
    if size(node,1)<2e3
        figure(1); 
        showresult(node,elem,pde.uexact,uh);
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
disptable(colname,NNdof,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order), L2-norm
% (2nd order) is observed. 