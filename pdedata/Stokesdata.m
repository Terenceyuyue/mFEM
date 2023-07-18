function pde = Stokesdata()


% --------- given by the symbolic computation ------
[u1,u1x,u1y,u2,u2x,u2y,pe,px,py,f1,f2] = compute_rhs();

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y),u2(x,y)];
    end
    function val = pexact(p)
        x = p(:,1); y = p(:,2);
        val = pe(x,y);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [f1(x,y),f2(x,y)];
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% Neumann boundary conditions ( right side hand )
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [u1x(x,y), u1y(x,y),u2x(x,y), u2y(x,y)];
    end

   function val = Dpexact(p)
        x = p(:,1); y = p(:,2);
        val = [px(x,y), py(x,y)];
    end

pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du,'pexact', @pexact, 'Dpexact',@Dpexact);
end

function [u1,u1x,u1y,u2,u2x,u2y,p,px,py,f1,f2] = compute_rhs()
    syms x y;
    % exact solution
    u1 = -(2^8)*(x^2-2*x^3+x^4)*(2*y-6*y^2+4*y^3);
    u2 = 2^8*(2*x-6*x^2+4*x^3)*(y^2-2*y^3+y^4);
    p = -(2^8)*(2-12*x+12*x^2)*(y^2-2*y^3+y^4);
    % derivative
    u1x = diff(u1,x);      u1y = diff(u1,y);
    u2x = diff(u2,x);      u2y = diff(u2,y);
    px = diff(p,x);        py = diff(p,y);
    % Lap(u)
    Lapu1 = diff(u1,x,2)+diff(u1,y,2);
    Lapu2 = diff(u2,x,2)+diff(u2,y,2);
    % f
    f1 = -Lapu1 + px;  f2 = -Lapu2 + py;
    
    % convert to anonymous functions
    u1 = matlabFunction(u1,'Vars',{x,y}); 
    u2 = matlabFunction(u2,'Vars',{x,y});  
    p = matlabFunction(p,'Vars',{x,y});
    
    u1x = matlabFunction(u1x,'Vars',{x,y}); u1y = matlabFunction(u1y,'Vars',{x,y}); 
    u2x = matlabFunction(u2x,'Vars',{x,y}); u2y = matlabFunction(u2y,'Vars',{x,y});
    px = matlabFunction(px,'Vars',{x,y}); py = matlabFunction(py,'Vars',{x,y});
    
    f1 = matlabFunction(f1,'Vars',{x,y});    
    f2 = matlabFunction(f2,'Vars',{x,y}); 
end
