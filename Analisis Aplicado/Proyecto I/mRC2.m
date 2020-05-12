function [x, msg, k] = mRC2(f, x0, itmax)
% M�todo de regi�n de confianza usando punto Dogleg
% In : f ... (handle) funcion a optimizar
% x0 ... (vector columna) punto inicial
% itmax ... (numero natural) cota superior para el n�mero de iteraciones
%
% Out: x ... (vector) ultima aproximacion del punto estacionario
% msg ... (string) Mensaje que informa si se encontr� o no el m�nimo
xk = x0;
eta = 0.1;
tol = 10^-5;
deltaMax = 1.5; 
deltak = 1;
k = 0;              % contador para iteraciones

gk = apGrad(f, xk);
Bk = apHess(f, xk);
n = length(x0); %Necesitamos la dimension del vector en caso de que Bk necesite shift
%Checamos si Bk es spd
eigMin = eigs(Bk,1, 'smallestabs');
if eigMin <= 0
    %Bk necesita un shift 
    Bk = Bk + eye(n)*(10^(-12) - 9.0*eigMin/8.0);
end

    while k < itmax && norm(gk) > tol 

        %Ahora encontramos el punto pk y el valor de rho (rok)
        pk = pDogLeg(Bk, gk, deltak);
        rok = (f(xk) - f(xk + pk)) / (f(xk) - (f(xk)+ gk'*pk +0.5*(pk'*Bk*pk)));
        norm_pk = norm(pk);
        if rok < 0.25
                deltak = 0.25 * deltak;
        else
            if rok > 0.75 && norm_pk == deltak
                deltak = min(2*deltak, deltaMax); %hacemos m�s grande la RC
            end
        end   
        if rok > eta
            %Tenemos que actualizar xk, luego tambien actualizamos gk y Bk pues
            %son dependientes de xk. Por �ltimo, checamos que Bk sea spd
            xk = xk + pk;
            gk = apGrad(f, xk);
            Bk = apHess(f, xk);
            eigMin = eigs(Bk,1, 'smallestabs');
            if eigMin <= 0
                %Bk necesita un shift 
                Bk = Bk + eye(n)*(10^(-12) - 9.0*eigMin/8.0);
            end
        end
        k = k + 1;
        
    end
    x = xk;
    
    %Checamos si el m�todo convirgi� antes de alcanzar el n�mero m�ximo de
    %iteraciones:
    if norm(gk) < tol 
        %Checamos si Bk es spd pues esto nos dir� si se trata de un m�nimo
        eigMin = eigs(Bk,1, 'smallestabs');
        if eigMin >= 0
            msg = "Se encontr� el m�nimo: " + x + " en " + k + "iteraciones";
        end
    else
        msg = "El metodo no convergio a un m�nimo y se detuvo tras " ...
            + k + "iteraciones";
    end
end
