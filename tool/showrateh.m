function showrateh(h,varargin)
%showrateh displays convergence rates of an err sequence
%
%  e.g. showrateh(h,ErrL2,ErrH1,ErrH2,ErrI)
%       showrateh(h,ErrL2,'||u-uh||');
% Copyright (C) Terence Yu.

%% Determine number of Err
n = 0;
for i = 2:nargin
    varName = lower(inputname(i));
    if mycontains(varName, 'err'),  n = n + 1; end
end

%% Define the line styles
lab1 = {'r-*', 'b-s', 'b-o', 'k-d'};
lab2 = {'k.-', 'k--', 'r.-', 'k:.'};

%% Plot convergence rate according to self-defined strings
if nargin-1 > n
    str = cell(1,n);
    for i = 1:n
        r = showrate(h,varargin{i},lab1{i},lab2{i});
        str{2*i-1} = varargin{n+i};
        str{2*i} = ['O (h^{' num2str(r,2) '})'];
    end
    h_legend = legend(str,'location','best');
    set(h_legend,'FontSize',10);
    return;
end

%% Plot convergence rates according to the inputnames
str = cell(1,2*n);  
for i = 1:n
    var = upper(inputname(i+1));
    if mycontains(var, 'L2') 
        stri = '||u-u_h||';  
        r = showrate(h,varargin{i},lab1{1},lab2{1});
    end
    if mycontains(var, 'H1')
        stri = '|u-u_h|_1';  
        r = showrate(h,varargin{i},lab1{2},lab2{2});
    end
    if mycontains(var, 'H2')
        stri = '|u-u_h|_2'; 
        r = showrate(h,varargin{i},lab1{3},lab2{3});
    end
    if mycontains(var, 'I')
        stri = '||u_I-u_h||_E'; 
        r = showrate(h,varargin{i},lab1{4},lab2{4});        
    end
    str{2*i-1} = stri;
    str{2*i} = ['O (h^{' num2str(r,2) '})'];
end

h_legend = legend(str,'location','best');
set(h_legend,'FontSize',10);