function [A,b,c] = generaKleeMinty(m)
% purpose: Version del Simplex en la Fase II
% minimizar c^T x
% sujeto a Ax <= b , x >= 0 , b >= 0
% In : m ... indica la dimension de la matriz A
%            del vecor b, y del vector c.
%
% Out: A ... matriz de restricciones
%      b ... vector de restricciones tal que Ax <= b
%      c ... vector de costo

% Generamos la Matriz A
A = ones(m); 
A = 2*tril(A)- eye(m);

% Generamos el vector b
b = zeros(m,1);
    for p = 1:m
        b(p) = 2^(p) -1;
    end
    
%Generamos el vector c
c = ones(m,1);
end

