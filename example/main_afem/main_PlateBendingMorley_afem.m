clc;clear;close all;

%% Parameters
maxN = 1e4;     theta = 0.99;    maxIt = 100;
a1 = 0; b1 = 1; a2 = 0; b2 = 1;
ErrL2 = zeros(1,maxIt);

%% Generate an initial mesh
Nx = 4; Ny = 4; h1 = 1/Nx; h2 = 1/Ny;
[node,elem] = squaremesh([a1 b1 a2 b2],h1,h2);

%% Get the PDE data
pde = PlateBendingData4;

%% Adaptive Finite Element Method 
format shorte
for k = 1:maxIt
    % Step 1: SOLVE
    bdStruct = setboundary(node,elem);
    [wh,he] = PlateBendingMorley(node,elem,pde,bdStruct);
    figure(1); 
    showresult(node,elem,pde.uexact,wh);
    pause(0.05);
    ErrL2(k) = getL2error_Morley(node,elem,pde.uexact,wh);
    
    % Step 2: ESTIMATE
    eta = PlateBendingMorley_indicator(node,elem,wh,pde);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = bisect(node,elem,elemMarked);    
    
    if (size(node,1)>maxN) || (k==maxIt) || min(he)<1e-8
        bdStruct = setboundary(node,elem);
        wwh = PlateBendingMorley(node,elem,pde,bdStruct);
        step = k
        break;
    end    
end
figure,
plot(1:k,ErrL2(1:k),'k.-','linewidth',1);
