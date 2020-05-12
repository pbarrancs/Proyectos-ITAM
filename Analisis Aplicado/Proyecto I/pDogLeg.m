function [p] = pDogLeg(B, g, delta)
% In    : B ... matriz s.p.d que aproxima la hessiana de f en xk
%         g ... (vector) gradiente de f en xk
%       delta ... radio de región de confianza
%
% Out: p ... Punto Dogleg


alfa = ((g)'*(g))/(g'*B*g);
normag = norm(g);
%% Primero calculamos Pu
Pu = -alfa*g;

%% Checamos condicion para Pu y Pb
    if alfa >= delta/normag
        p = -((delta)/normag) * g; %Dogleg es Cauchy
    else
        Pb = -inv(B)*g; %Dirección de Newton
        if norm(Pb) <= delta
            p = Pb; %Dogleg es la dirección de Newton
        else
            %Coeficientes del polinomio cuadrático respecto a alfa
            alfaroots = [dot((Pb-Pu),(Pb-Pu));
                          2*dot(Pu,(Pb-Pu));
                          (dot(Pu,Pu) -(delta*delta))];
            %buscamos alfa entre 0 y 1 como mínimo del polinomio anterior
            r = roots(alfaroots);
            if r(1) > 0 && r(1) < 1
                p = Pu + r(1)*(Pb - Pu);
            else
                p = Pu + r(2)*(Pb - Pu);
            end
        end
    end
    
end

