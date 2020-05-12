function [pC] = pCauchy( B, g, delta )
% In : B ... (symmetric matrix) approximates the hessian of f in xk
% g ... (vector) gradient of f in xk
% delta ... trust region radius
%
% Out: pC ... The Cauchy point

    
    pk = g/ norm(g);
    %Calculamos apha 
    valor = dot(g,B*g);
    if valor > 0
        alphaStar = norm(g)/(delta*valor);
    else
        alphaStar = 1;
    end
    pC = -delta*min(1,alphaStar)*pk;
    
end

