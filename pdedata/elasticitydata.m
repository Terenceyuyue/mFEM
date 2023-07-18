function pde = elasticitydata(para)

% ------ Lame constants ------
if nargin == 0
    lambda = 1;   mu = 1;
else
    lambda = para.lambda; mu = para.mu;
end

pde = struct('lambda',lambda,'mu',mu, 'f', @f, 'uexact',@uexact, 'g_D',@g_D, 'g_N',@g_N, 'Du',@Du);

    function val = f(p)
        val = mu*2*pi^2.*uexact(p);
    end
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [cos(pi*x).*cos(pi*y), sin(pi*x).*sin(pi*y)];
    end
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        val(:,1) = -pi*sin(pi*x).*cos(pi*y); % u1x
        val(:,2) = -pi*cos(pi*x).*sin(pi*y); % u1y
        val(:,3) = pi*cos(pi*x).*sin(pi*y);  % u2x
        val(:,4) = pi*sin(pi*x).*cos(pi*y);  % u2y
    end
    function val = g_D(p)
        val = uexact(p);
    end
    function val = g_N(p)
        x = p(:,1); y = p(:,2);
        E11 = -pi*sin(pi*x).*cos(pi*y); E22 = pi*sin(pi*x).*cos(pi*y);
        E12 = 0.5*(-pi*cos(pi*x).*sin(pi*y) + pi*cos(pi*x).*sin(pi*y));
        sig11 = (lambda+2*mu)*E11 + lambda*E22;
        sig22 = (lambda+2*mu)*E22 + lambda*E11;
        sig12 = 2*mu*E12;
        val = [sig11,sig22,sig12];
    end
end