function [xk, iter] = mRCSR1(f, x0,tol,maxiter, radioMax)
% Trust region method decribed on the previous page
%
% In : f ... (handle) function to be optimized
% x0 ... (vector) initial point
% itmax ... (natural number) upper bound for number of iterations
%
% Out: x ... (vector) last approximation of a stationary point
% iter... (natural number) iterations used
% msg ... (string) message that says whether (or not) a minimum was found


    eta = 0.1;
    r = 1e-6;
    n = length(x0);
    radio = radioMax;
    iter = 0;
    xk = x0;
    g = apGrad(f, xk);
    H = speye(n);
    B = H;
    
    while norm(g, 'inf') > tol
        %P1
        s = -H*g;
        if dot(s,g) < 0
            if norm(s) > radio
                s = radio*s/norm(s);
            end
        else
            s = pCauchy( B, g, radio);
        end
         
        %P2
        coc = -(f(xk) - f(xk+s))/(dot(g,s) + 0.5*dot(s,B*s));
        gnew = apGrad(f, xk + s);
        gamma = gnew -g;
        
        %P3
        if(coc > eta)
            xk = xk + s;
            g = gnew;
            iter = iter + 1;
            if iter == maxiter
                break
            end
        end
        
        if coc > 0.75 %P4
            if norm(s) > 0.8*radio
                radio = min(2*radio, radioMax);
            end
        elseif coc < 0.1 %P5
            radio = 0.5*radio;
        end
            
        v = gamma -B*s;
        %P6
        if abs(dot(v,s)) >= r*norm(v)*norm(s)
            B = B + v/dot(v,s)*v';
            u = s-H*gamma;
            H = H + u/dot(u,gamma)*u';
        end
    end
end 
