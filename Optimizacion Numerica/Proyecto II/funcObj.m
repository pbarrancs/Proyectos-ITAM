function [f] = funcObj(x)
% Se divide el vector en la parte superior tiene las r's y en la parte
% inferior las thetas. Otra version podria ser pensandolo las entradas con
% una Matriz de 2xn
    
    %% Inicializamos 
    np = length(x)/2; 
    r = x(1:np);
    theta = x(np+1:2*np);
    f = 0;
    
    %% Suma
    for i = 1:(np-1)
        f = f + r(i)*r(i+1)*sin(theta(i+1)-theta(i));
    end
    
    %% factorizamos el (1/2) de la suma y como estamos maximizando f*(-1)
    f = (-0.5)*f; 
end