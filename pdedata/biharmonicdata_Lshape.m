function pde = biharmonicdata_Lshape()


% --------- given by the symbolic computation ------
[u,ux,uy,rhs] = compute_rhs();

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = u(x,y);
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = rhs(x,y);
    end
% Dirichlet boundary conditions
    function val = g_D(p)
        val = uexact(p);
    end
% derivative
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val = [ux(x,y), uy(x,y)];
    end


pde = struct('uexact',@uexact, 'f',@f, 'g_D',@g_D, 'Du', @Du);
end

function [u,ux,uy,f] = compute_rhs()
    syms x y;

%     a = 0.5444837;
%     w = 3*pi/2;
% 
%     % polar
%     r = sqrt(x^2+y^2);
%     t = atan2(y,x);
% 
%     % g
%     g = (sin((a-1)*w)/(a-1) - sin((a+1)*w)/(a+1))*(cos((a-1)*t) - cos((a+1)*t)) ...
%         - (sin((a-1)*t)/(a-1) - sin((a+1)*t)/(a-1))*(cos((a-1)*w) - cos((a+1)*w));
% 
%     % exact solution
%     u = (1-r^2*x^2)^2*(1-r^2*y^2)^2*r^(1+a)*g;

    u = x^2*y^2*(1-x)^2*(1-y)^2*sin(pi*x);
    % derivative
    ux = diff(u,x);      uy = diff(u,y);    
    u1111 = diff(u,x,4);
    u1212 = diff(diff(u,x,2),y,2);
    u2121 = u1212;
    u2222 = diff(u,y,4);
    % f
    f = u1111+u1212+u2121+u2222;
    
    % convert to anonymous functions
    u = matlabFunction(u,'Vars',{x,y});     
    ux = matlabFunction(ux,'Vars',{x,y}); uy = matlabFunction(uy,'Vars',{x,y}); 
    
    f = matlabFunction(f,'Vars',{x,y});    
end
