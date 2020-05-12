function [xk, iter] = lineBGFS( f, x0, tol, maxiter )
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
%
% Same parameters and results as lineDFP,
% with the exception that the (iBGFS) update formula is used.
%In : f ... function to minimize
% x0 ... initial point
% tol ... tolarance
% maxiter ... upper bound for iterations
%
% Out: xf ... final approximation of x*
% iter ... number of iterations used
%
% Initial values: H = I
% Criterios de paro: || gk || <= tol or ||s|| <= 1e-7


    n = length(x0);
    iter = maxiter;
    xk = x0;
    g = apGrad(f, xk);
    H = speye(n);
    
    for k = 1:maxiter
        if norm(g,'inf') <= tol 
            iter = k-1;
            break
        end
        dk = -H*g;%Paso 1
        assert( dot(dk, g) < 0 )
        [alpha, gnew] = lineSearch( f, xk, dk, g );%paso 2
        
        s = alpha*dk;
        xk = xk + s;
        gamma = gnew - g;
        irho = dot(gamma,s); % inversa de rho
        
        assert(irho >0)
        Hgam = H*gamma/irho;
        H = H - (s*Hgam' + Hgam*s') + ((dot(gamma, Hgam)+1)/irho*s)*s';
        g = gnew;
    end
    
    
end