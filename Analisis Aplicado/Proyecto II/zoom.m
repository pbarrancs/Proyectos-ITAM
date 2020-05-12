function [alphaOptimo] = zoom(alphaLo, alphaHi, f, xk, gk, dk)

%Define los parametros
c1 = 1e-4;
c2 = 0.99;

slope0 = dot(gk,dk);

%Definimos nuestras funciones 
Phi = @(x) f(xk + x*dk);
PhiP = @(y) dot(apGrad(f,xk + y*dk),dk);
Line = @(z) f(xk) + c1*z*slope0;

%Algoritmo para zoom
while 1
    alphaB = (alphaLo + alphaHi)/2;
    if Phi(alphaB) > Line(alphaB) || Phi(alphaB) >= Phi(alphaLo)
        alphaHi = alphaB;
    elseif abs(PhiP(alphaB)) <= -c2*slope0
        alphaAux = alphaB;
        break
    elseif PhiP(alphaB)*(alphaHi - alphaLo) >= 0
        alphaHi = alphaLo;
        alphaLo = alphaB;
    else
        alphaLo = alphaB;
    end
end

alphaOptimo = alphaAux;

end
