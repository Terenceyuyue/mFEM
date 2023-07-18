clc;clear;close all;

%% Parameters
maxN = 1e4;     theta = 0.4;    maxIt = 100;
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
ErrL2 = zeros(1,maxIt);

%% Generate an initial mesh
Nx = 5; Ny = 5; h1 = 1/Nx; h2 = 1/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

%% Get the PDE data
pde = Poissondata_afem();

%% Adaptive Finite Element Method 
for k = 1:maxIt
    % Step 1: SOLVE
    bdStruct = setboundary(node,elem);
    option.solver = 'direct';
    u = Poisson(node,elem,pde,bdStruct,option);
    figure(1); 
    showresult(node,elem,pde.uexact,u);
    pause(0.05);
    ErrL2(k) = getL2error(node,elem,pde.uexact,u);
    
    % Step 2: ESTIMATE
    eta = Poisson_indicator(node,elem,u,pde);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = bisect(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt)
        bdStruct = setboundary(node,elem);
        u = Poisson(node,elem,pde,bdStruct,option);
        setp = k
        break;
    end    
end
figure,
plot(1:k,ErrL2(1:k),'k.-','linewidth',1);
