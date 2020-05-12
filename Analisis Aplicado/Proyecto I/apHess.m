function [D2f] = apHess(f, x)
    n = length(x);
    H = spdiags( nthroot(eps, 4)*(abs(x)+1), 0, n,n);
    
    D2f = zeros(n,n);
    % parte triangular inferior
    for j = 1:n-1
        for i = j+1:n
            D2f(i,j) = 0.25*( f(x+H(:,i)+H(:,j)) -f(x-H(:,i)+H(:,j)) ...
                              -f(x+H(:,i)-H(:,j)) + f(x-H(:,i)-H(:,j)) )/(H(i,i)*H(j,j));
        end
    end
    % parte superior
    D2f = D2f+D2f';
    
    % diagonal
    for i = 1:n
       D2f(i,i) = 0.25*( f(x+2*H(:,i)) -2*f(x) + f(x-2*H(:,i)) )/(H(i,i)*H(i,i));
    end
end