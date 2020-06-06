function [g] = rest(x)
% Descripcion


    %% Inicializamos
    n = length(x)/2;
    r = x(1:n);
    theta = x((n + 1):2*n);
    g = [];

    %% Primer Restriccion
    %Combinaciones (i,j)
    for i = 1:(n-1)
        for j = i+1:n
            g = [g; 1-(r(i)^2 + r(j)^2 - 2*r(i)*r(j)*cos(theta(i)-theta(j)))]; %Concatemanos renglones
        end
    end
    
    %% Segunda Restricción
    for i = 1:n
        g = [g; r(i)];
    end
    %Opcional para fijar el polinomio entre 0 y 1
    for i = 1:n
        g = [g;1 - r(i)];
    end
    
    %% Tercera Restricción
    %Parte 1
    for i = 1:n
        g = [g; pi - theta(i)];
    end
    %Parte 2
    for i = 1:n
        g = [g; theta(i)];
    end

    %% Cuarta Restricción
    for i = 1:n-1
        g = [g; theta(i+1) - theta(i)];
    end
    
    %% Fijo el origen (opcional)
    g = [g; - theta(1)];
    g = [g; -r(1)];
end