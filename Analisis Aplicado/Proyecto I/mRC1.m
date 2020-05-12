function [x, msg, k] = mRC1(f, x0, itmax)
% Método de región de confianza usando el punto de Cauchy.
%
% In :  f ... (handle) funcion por optimizar
%      x0 ... (vector) punto inicial
%   itmax ... (natural number) cota superior para el número de iteraciones.
%
% Out:  x ... (vector) última aproximación de un punto estacionario
%     msg ... (string) Mensaje que dice si se encontró on mínimo o no


xk = x0;
eta = 0.1;
tol = 10^-5;
deltaMax = 1.5; 
deltak = 1;
k = 0;

Bk = apHess(f, xk);
gk = apGrad(f, xk);

    while k < itmax && norm(gk) > tol 
        
        pk = pCauchy(Bk, gk, deltak);
        rok = (f(xk) - f(xk + pk)) / (f(xk) - (f(xk)+ gk'*pk +0.5*(pk'*Bk*pk)));
        norm_pk = norm(pk);
        
        if rok < 0.25
            deltak = 0.25*(norm_pk);
        else
            if rok > 0.75 && norm_pk == deltak
                deltak = min(2*deltak, deltaMax);
            end
        end
        
        if rok > eta
            xk = xk + pk;
            Bk = apHess(f, xk);
            gk = apGrad(f, xk);
        end
        k = k + 1;
    end
    %%AHORA CREAMOS EL MENSAJE%%
    
    x = xk;
    %Checamos si el método convirgió antes de alcanzar el número máximo de
    %iteraciones:
    if norm(gk) < tol 
        %Checamos si Bk es spd pues esto nos dirá si se trata de un mínimo
        eigMin = eigs(Bk,1, 'smallestabs');
        if eigMin >= 0
            msg = "Se encontró el mínimo: " + x + " en " + k + "iteraciones";
        end
    else
        msg = "El metodo no convergio a un mínimo y se detuvo tras " ...
            + k + "iteraciones";
    end
end

