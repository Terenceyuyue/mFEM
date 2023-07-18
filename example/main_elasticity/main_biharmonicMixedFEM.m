clc; clear; close all;
%% Parameters
maxIt = 5; 
NNdof = zeros(maxIt,1);  h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);  erruH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);  errwH1 = zeros(maxIt,1);

%% Generate an initial mesh
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
Nx = 4; Ny = 4; h1 = (b1-a1)/Nx; h2 = (b2-a2)/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

%% Get the PDE data
pde = biharmonicdata;

%% Finite element method
for k = 1:maxIt
    % refine mesh
   [node,elem] = uniformrefine(node,elem); 
   % set boundary
   bdStruct = setboundary(node,elem);
   % solve the equation
   [u,w] = biharmonicMixedFEM(node,elem,pde,bdStruct);
   % record and plot
   NNdof(k) = length(u); 
   NT = size(elem,1); h(k) = 1/sqrt(NT);
   if NT < 2e3
       figure(1);
       showresult(node,elem,pde.uexact,u);
       pause(1);
   end
   % compute error
   erruL2(k) = getL2error(node,elem,pde.uexact,u);
   erruH1(k) = getH1error(node,elem,pde.Du,u);
   errwL2(k) = getL2error(node,elem,pde.wexact,w);
   errwH1(k) = getH1error(node,elem,pde.Dw,w);
end

%% Plot convergence rates and display error table
figure(2);
subplot(1,2,1);
showrateh(h, erruL2, erruH1, '||u - u_h||', '|u - u_h|_1');
subplot(1,2,2)
showrateh(h, errwL2, errwH1, '||w - w_h||', '|w - w_h|_1');    

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','||w-w_h||','|w-w_h|_1'};
disptable(colname,NNdof,[],h,'%0.3e',erruL2,'%0.5e',erruH1,'%0.5e',errwL2,'%0.5e',errwH1,'%0.5e');

%% Conclusion
%
% In general, the rate of convergence for u is optimal but for w is sub-optimal. 
% For linear element, optimal order for w is also observed.