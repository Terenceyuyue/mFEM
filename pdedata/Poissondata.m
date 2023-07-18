function pde = Poissondata()

pde = struct('uexact',@uexact, 'f',@f, 'Du',@Du,  'g_D',@g_D);

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = y.^2.*sin(pi*x);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = (pi^2*y.^2-2).*sin(pi*x);
    end
% Neumann boundary conditions (given by Du)
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [pi*y.^2.*cos(pi*x), 2*y.*sin(pi*x)];
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
end
