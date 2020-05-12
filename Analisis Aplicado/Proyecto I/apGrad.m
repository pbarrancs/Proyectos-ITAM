function [gradf] = apGrad( f, x )
    n = length(x);
    H = spdiags( nthroot(eps, 3)*(abs(x)+1), 0, n,n);
    
    gradf = zeros(n,1);
    for i = 1:n
        gradf(i) = ( f(x+H(:,i)) - f(x-H(:,i)) )*0.5/H(i,i);
    end
end