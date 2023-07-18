clc; clear; close all;
%% Parameters
maxIt = 5;
h = zeros(maxIt,1); NNdof = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Get the PDE data
a = 1;  b = 1;  c = 0;
para = struct('a',a, 'b',b, 'c',c);
pde = pde1D(para);

%% Finite element method
for k = 1:maxIt % level
    % generate an initial mesh
    a = 0; b = 1; N0 = 5;  
    node = linspace(a,b,N0)'; % uniform
    elem = [(1:N0-1)', (2:N0)']; % initial mesh
    % Final mesh and prolongation matrices
    [node,elem,Pro,Res] = Mesh1DVcycle(node,elem,k);
    % set boundary
    bdNeumann = 'abs(x-1)<1e-4';
    bdStruct = setboundary1(node,elem,bdNeumann);
    % solve the equation
    uh = FEM1DVcycle(node,elem,pde,bdStruct,Pro,Res);
    % record
    NNdof(k) = length(uh);
    h(k) = 1/size(elem,1);
    if NNdof(k) < 80  % show solution for small size
        figure(1); 
        [node1,id] = sort(node);
        plot(node1,uh(id),'r-',node1,pde.uexact(node1), ...
            'k*','linewidth',2);
        pause(1);
    end
    % compute error
    ErrL2(k) = getL2error1(node,elem,pde.uexact,uh);
    ErrH1(k) = getH1error1(node,elem,pde.Du,uh);
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

