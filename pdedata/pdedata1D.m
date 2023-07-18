function pde = pdedata1D

% --------- given by the symbolic computation ------
[af,bf,cf,u,ux,rhs] = compute_rhs;

% a
    function val = a(p)
        x = p(:,1);
        val = af(x) + 0*x;
    end
% b
    function val = b(p)
        x = p(:,1);
        val = bf(x) + 0*x;
    end
% c
    function val = c(p)
        x = p(:,1);
        val = cf(x) + 0*x;
    end
% exact solution
    function val = uexact(p)
        x = p(:,1);
        val = u(x);
    end

% right side hand function
    function val = f(p)
        x = p(:,1);
        val = rhs(x);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% for Neumann boundary conditions
    function val = Du(p)
        x = p(:,1);
        val = ux(x);
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du',@Du, 'a',@a, 'b',@b, 'c',@c);

end

function [a,b,c,u,ux,f] = compute_rhs

    syms x;
    c1 = 0.5/exp(1); c2 = -0.5*(1+1/exp(1));
    u = c1*exp(2*x)+c2*exp(x)+1/2;

    ux = diff(u);  

    a = 1 + 0*x;  b = 1 + 0*x;  c = 1 + 0*x;
    f = -a*diff(u,2)+b*diff(u,1)+c*u; % note: -a

    % convert to anonymous functions
    a = matlabFunction(a,'Vars',{x});
    b = matlabFunction(b,'Vars',{x});
    c = matlabFunction(c,'Vars',{x}); 
    u = matlabFunction(u,'Vars',{x});    
    f = matlabFunction(f,'Vars',{x});    
    ux = matlabFunction(ux,'Vars',{x});  

end