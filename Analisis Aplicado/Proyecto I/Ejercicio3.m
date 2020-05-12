clear;   close all;   clc;

%% Define f with two arguments  and  so that
% it can be evaluated for matrices of values X, Y
f = @(x,y) (x + 2*y - 7).^2 + (2*x + y - 5).^2;
f2 = fbooth();


%% define point and trust region radius
x0    = [5; 5.5];
gk = apGrad(f2,x0);
Bk = apHess(f2,x0);



%% Level sets for Cauchy point
stepsize =  0.1;  % mas chico hace los conjuntos de nivel mas detallados
[X,Y] = meshgrid(-1:stepsize:7);
z = f(X,Y);
niveles = [1, 1:13];
figure(1)
contour(X,Y,z, niveles)
axis equal


%% Iteration History for Cauchy 
hold on
deltaMax = 1.5; 
deltak = 1;
k = 0;
xk =  x0;
eta = 0.1;
B = zeros([1 13]);
C = zeros([1 13]);
for i = 1:14
    pk = pCauchy(Bk, gk, deltak);
        rok = (f2(xk) - f2(xk + pk)) / (f2(xk) - (f2(xk)+ gk'*pk +0.5*(pk'*Bk*pk)));
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
            Bk = apHess(f2, xk);
            gk = apGrad(f2, xk);
        end
        B(i) = pk(1) + xk(1);
        C(i) = pk(2) + xk(2);
        k = k + 1;
end

plot(B,C,'--d')
title('Cauchy Point')
hold off



%% Level sets for DogLeg point
stepsize =  0.1;  % mas chico hace los conjuntos de nivel mas detallados
[X,Y] = meshgrid(-1:stepsize:7);
z = f(X,Y);
niveles = [1, 1:13];
figure(2)
contour(X,Y,z, niveles)

axis equal



%% Iteration History for Dogleg
% Reiniciamos nuestras variables
E = zeros([1 4]);
F = zeros([1 4]);
xk = x0;
delta = 1.5;
gk = apGrad(f2,x0);
Bk = apHess(f2,x0);
for i = 1:5

        %Ahora encontramos el punto pk y el valor de rho (rok)
        pk = pDogLeg(Bk, gk, deltak);
        rok = (f2(xk) - f2(xk + pk)) / (f2(xk) - (f2(xk)+ gk'*pk +0.5*(pk'*Bk*pk)));
        norm_pk = norm(pk);
        if rok < 0.25
                deltak = 0.25 * deltak;
        else
            if rok > 0.75 && norm_pk == deltak
                deltak = min(2*deltak, deltaMax); %hacemos más grande la RC
            end
        end   
        if rok > eta
            %Tenemos que actualizar xk, luego tambien actualizamos gk y Bk pues
            %son dependientes de xk. Por último, checamos que Bk sea spd
            xk = xk + pk;
            gk = apGrad(f2, xk);
            Bk = apHess(f2, xk);
            eigMin = eigs(Bk,1, 'smallestabs');
            if eigMin <= 0
                %Bk necesita un shift 
                Bk = Bk + eye(n)*(10^(-12) - 9.0*eigMin/8.0);
            end
        end
        E(i) = pk(1) +xk(1);
        F(i) = pk(2) + xk(2);
        k = k + 1;
end
hold on
plot(E,F,'r--d')
title('Dogleg')
hold off