function [x, fx, iter, kkt] = metpuntint(funobj, fung, x)
% Aproxima una solución del problema
%
% Min f(x)
% Sujeto a g(x) >= 0 .
%
% donde las funciones, f : Rn -> R y g : Rn -> Rp
% por medio del método de puntos interiores y
% y actualización cuasi-Newton a la matriz hessiana del lagrangiano.
%
% In
% funobj: cadena de caracteres con el nombre del código de Matlab
%         con la expresión algebraica de f(x).
% fung: cadena de caracteres con el nombre del código de Matlab
%       con la expresión algebraica de g(x).
% x0: es un punto inicial para el proceso de optimización que
%     satisface, g(x0) > 0.
% Out
% x: aproximación al minimo del problema.
% fx: valor de la función objetivo en el optimo x.
% iter: número de iteraciones en el proceso.
%--------------------------------------------------------------------------
% Parámetros y valores dentro de la función:
% maxiter = 200, número máximo de iteraciones permitidas
% tol = 10^-5, criterio de paro a las condiciones de Karush-Kuhn-Tucker
% mu, vector de multiplicadores de Lagrange inicial igual al vector de unos
% z, vector de holgura inicial igual al vector de unos
% In matriz inicial para la Hessiana del Lagrangeano en xo
%-------------------------------------------------------------------------
% Optimización Numérica
%         ITAM
% Pablo Barranco Soto 151528
% Jose Miguel Carmona Jiménez 157600
% Diego Villegas Aguilar 143456
%--------------------------------------------------------------------------
% Comentarios:
%   Si el sistema lineal se malcondiciona por alguna entrada de la inversa de
%   Z, usamos el sistema completo para obtener mejores resutltados. Para
%   mas informacion, consulte nuestro pdf.
%--------------------------------------------------------------------------
%% Inicilizamos nuestras variables
% Hiperparametros
tol = 1e-05;
maxiter = 200;
Ftol = 1e-06;
bandera = true;

%Valores Iniciales
iter = 0;
n = length(x);
fx = feval(funobj,x);
gx = feval(fung,x); %g(x)
p = length(gx);
mu = ones(p,1); % Multiplicador de Lagrange
z = ones(p,1);% Variable de Holgur
gamma = 1.0; % Variable de Barrera Logarítmica
             
grad = gradiente(funobj,x); %d/dx f(x)
A = jacobiana1(fung,x); % A(x)
Z = diag(z);
U = diag(mu);
e = ones(p,1);
B = eye(n);

%Condiciones Necesarias de Primer Orden con hariable de holgura
v1 = grad - A'*mu;
v2 =  -gx + z;
v3 = U*(Z*e);
cnpo = [v1; v2;v3];
kkt = norm(cnpo);

%para graficar 
iters = iter;
cnpos = kkt; 
fxs = -fx;
%% Iteraciones del algoritmo
while(norm(cnpo,'inf') > tol && iter  < maxiter && bandera)   
%Sistema Lineal resolvemos Método de Newton
    C = Z\U;
    M = B + A'*(C*A);
    if( rcond(M) < 1e-5) %Usamos el sistema Completo si esta mal condicionado
         M = [B           -A'       zeros(n,p);...
        -A         zeros(p)       eye(p);...
        zeros(p,n)      Z             U];
    
        ld = -[v1; v2; v3 - gamma*e];
        dw = M\ld;
        dx = dw(1:n);
        dm = dw(n+1:n+p);
        dz = dw(n+p+1:n+p+p);
    else % De lo contrario usamos el reducido
        ld = -v1 - (A'*C)*gx + gamma*A'*(Z\e);
        dx = M\ld;
        dm = gamma*(Z\e) - C*(A*dx + gx);
        dz = A*dx + gx - z;
    end
    
    %% Encontramos alpha (paso de recorte)
    vu = ones(p,1);
    vz = ones(p,1);

    for i = 1:p
         if(dm(i) < 0)
             vu(i) = -(mu(i)/dm(i));
         end

         if(dz(i) < 0)
             vz(i) = -(z(i)/dz(i));
         end
    end
    alfau = min(vu);
    alfaz = min(vz);
    alfa = (0.95)*min([alfau alfaz 1]);
    
    %% Nuevo Punto
    x0 = x;                 %guardamos el viejo punto
    x = x + alfa*dx;
    mu = mu + alfa*dm;
    z = z + alfa*dz;

    %% Actualizamos
    grad0 = grad;
    grad = gradiente(funobj,x);
    A0 = A;
    
    A = jacobiana1(fung,x);
    Z = diag(z);
    U = diag(mu);
    fx = feval(funobj,x);
    gx = feval(fung,x);
    gamma = gamma/10;
    
    v1 = grad - A'*mu;
    v2 =  -gx + z;
    v3 = U*(Z*e);
    cnpo = [v1; v2;v3];
    
    
    %% BGFS
    sk = alfa*dx; %sk = x -x0;
    yk = (grad - A'*mu) - (grad0 - A0'*mu);
    w =0.2;
    Betha = sk'*(B*sk);
    
    if sk'*yk > 0
        B = B + ((yk*yk')/(sk'*yk)) - ((B*sk)*(sk'*B)/(Betha));
    else
        if sk'*yk >= w*Betha
            tk = 1;
        else
            tk = ((1-w)*Betha)/(Betha - sk'*yk);
        end
        rk = B*sk + tk*(yk - B*sk);
        B = B  + ((rk*rk')/(sk'*rk)) - ((B*sk*(sk'*B))/(Betha)) ;
    end
    
    if rcond(B) < 1e-5
        B = eye(n);
    end
    
    %% Actualizacion de otros parametros
    iter = iter + 1;
    kkt = norm(cnpo);
    
    %para graficar 
    iters =[iters; iter];
    cnpos =[cnpos;kkt]; 
    fxs = [fxs;-fx];
    
    %% Condicion de paro f(x)
    if abs(fx- feval(funobj, x0)) < Ftol
        bandera = false;
    end
end

%% Graficación de resultados
figure(2)
bar(iters,cnpos)
title('Condiciones Necesarias de Primer orden') 
legend({'cnpo'},'Location','northeast')

figure(3)
stem(iters,fxs)
title('Valor de f(x)') 
legend({'f(x)'},'Location','northwest')

end
