function pde = pde1D_variable(para)

if isa(para.acoef,'double'), para.acoef = @(x) para.acoef+0*x; end
if isa(para.bcoef,'double'), para.bcoef = @(x) para.bcoef+0*x; end
if isa(para.ccoef,'double'), para.ccoef = @(x) para.ccoef+0*x; end

pde = struct('f', @f, 'uexact', @u, 'g_D', @g_D, 'g_N', @g_N);

    function val = f(x)
        val = para.acoef(x).*du2(x) + para.bcoef(x).*du(x) ...
            + para.ccoef(x).*u(x);
    end
    function val = u(x)
        val = sin(pi*x);
    end
    function val = du(x)
        val = pi*cos(pi*x);
    end
    function val = du2(x)
        val = -pi^2*sin(pi*x);
    end
    function val = g_D(x)
        val = u(x);
    end
    function val = g_N(x)
        val = du(x);
    end

end

