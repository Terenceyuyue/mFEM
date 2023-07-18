function pde = pde1D(para)

syms x;
c1 = 0.5/exp(1); c2 = -0.5*(1+1/exp(1));
u = c1*exp(2*x)+c2*exp(x)+1/2;
% exact solution
uexact = eval(['@(x)',vectorize(u)]);  % transform to anonymous function

% rhs
f = -para.a*diff(u,2)+para.b*diff(u,1)+para.c*u;
f = eval(['@(x)',vectorize(f)]);

% boundary conditions
du = diff(u);  du = eval(['@(x)',vectorize(du)]);
Du = @(x) du(x); g_D = @(x) uexact(x);

pde = struct('f', f, 'uexact', uexact, 'g_D', g_D, 'Du', Du, 'para', para);