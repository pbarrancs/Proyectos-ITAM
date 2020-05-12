%La función booth tiene un único mínimo global en (1,3)
%La Hessiana es constante con eigenvalores 18 y 2
close all; clc;
format long

x0 = [0.5;2.1];
delta = 1;
f =  fbooth(); 
%Podria ser más eficiente dar la definicion de f directamente, sin embargo,
%dado el alcance de los ejercicios, no es imperativo y hace más cómoda la
%escritura del código.

%Definimos el modelo cuadrático: necesitamos la hessiana Bk y el gradiente
%g

gk = apGrad(f,x0);
Bk = apHess(f,x0);

mk = @(pk) f(pk) + gk' * pk + 0.5 * pk'*Bk*pk;

%Calculamos la dirección de Newton, Cauchy y Dogleg (checar esta parte)
pN = -inv(Bk)*gk;
pC = pCauchy(Bk,gk,delta);
pDog = pDogLeg(Bk,gk,delta);

%Graficamos el modelo en la región de confianza junto con las direcciones
%de Newton(verde), Cauchy(rojo) y dogleg(Azul).
%usamos coordenadas polares: r en [0, delta], theta en [0,2*pi]

fsurf(@(r,theta) x0(1)+r*cos(theta), @(r,theta)  x0(2)+r*sin(theta), @(r,theta) ...
    mk([r*cos(theta);r*sin(theta)]), [0,delta,0,2*pi])

hold on
dirN = quiver3(x0(1), x0(2), f(x0),pN(1), pN(2), 0, 0,'Color', 'green',   'LineWidth', 2);
dirC = quiver3(x0(1), x0(2), f(x0), pC(1), pC(2), 0, 0,'Color', 'red',   'LineWidth', 2);
dirDog = quiver3(x0(1), x0(2), f(x0), pDog(1), pDog(2), 0, 0,'Color', 'blue', 'LineWidth', 2);

regionConfianza = viscircles(x0', delta,'LineStyle', '--', 'Color', 'm', 'LineWidth',1);

legend([dirN, dirC, dirDog, regionConfianza], {'Dirección Newton', ... 
    'Dirección Cauchy', 'Dirección Dogleg', 'Región de confianza'});

hold off
grid on