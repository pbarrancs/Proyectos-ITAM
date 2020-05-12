function [alphaOptimo, gnew] = lineSearch(f, xk, dk, gk)
%% Parametros
alphaOld = 0;
alphaNew = 1;
alphaMax = 500;
c1 = 1e-4;
c2 = 0.99;

%%
%Definimos handle functions para realizar el algoritmo
slope = dot(gk,dk);
Phi = @(x) f(xk + x*dk);
Leg = @(y) f(xk) + c1*y*slope; 
PhiP = @(z) dot(apGrad(f, xk + z*dk),dk);

%%
%Algoritmo para encontrar un alpha de manera sofisticada
while alphaNew > 0 && alphaNew < alphaMax
    if Phi(alphaNew) > Leg(alphaNew) || Phi(alphaNew) >= Phi(alphaOld)
        alphaAux = zoom(alphaOld, alphaNew, f, xk, gk, dk);
        break
        elseif abs(PhiP(alphaNew)) <= -c2*slope
        alphaAux = alphaNew;
        break
        elseif PhiP(alphaNew) >= 0
            alphaAux = zoom(alphaNew, alphaOld, f, xk, gk, dk);
        break
    else
        alphaOld = alphaNew;
        alphaNew = 2*alphaNew;
    end
end
    alphaOptimo = alphaAux;
    gnew = apGrad(f, xk + alphaOptimo*dk);
end