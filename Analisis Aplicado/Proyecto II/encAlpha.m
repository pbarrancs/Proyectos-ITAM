function [alpha, gnew] = encAlpha( f, xk, dk, gk )
% Purpose: find the first alpha of the form (1/2)^k that satisfies the Wolfe condition (W1).
% If (W2) is not satisfied, throw an error.
%
% In : f ... function to minimize
% xk ... current point
% dk ... descend direction
% gk ... gradient of f in xk
%
% Out: alpha ... parameter satisfying the two Wolfe conditions.
% gnew ... gradient of f in xk + alpha*dk
%
% parameters: 

    alpha = 1;
    c1 = 1e-4;
    c2 = 0.99;
    fk = f(xk);
    
    slope0 =  dot(gk,dk);
        while f(xk + alpha*dk) > fk + alpha*c1*slope0
            alpha = 0.5* alpha;
        end
        gnew = apGrad(f, xk + alpha*dk);
        assert( dot(gnew, dk) >= c2*slope0, 'W2  failed');
end