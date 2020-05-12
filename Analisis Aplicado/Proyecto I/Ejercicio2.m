%Nuestra Función es globalmente convexa
%La función booth tiene un único mínimo global en (1,3)
%La Hessiana es constante con eigenvalores 18 y 2
close all; clc;
format long



%deltaMax = 1.5;
%Definimos x0 de tal manera que la distancia con xOpt = (1,3) mayor q 3*deltaMax
max = [1;3];
x0 = [5; 5.5];
f =  fbooth();
%Podria ser más eficiente dar la definicion de f directamente, sin embargo,
%dado el alcance de los ejercicios, no es imperativo y hace más cómoda la
%escritura del código.


%Definimos el modelo cuadrático: necesitamos la hessiana Bk y el gradiente
%g

gk = apGrad(f,x0);
Bk = apHess(f,x0);


% Usaremos itMax muy grande para que no influya en ver que metodo se
% aproxima mas rapido.

itmax = 1000;

[x1, msg1, iter1] = mRC1(f, x0, itmax);
[x2, msg2, iter2] = mRC2(f, x0, itmax);

%Errores
error1 = norm(x1-max,'inf');
error2 = norm(x2-max,'inf');

%Errores relativos 
errorREL1 = norm(x1-max./max,'inf');
errorREL2 = norm(x2-max./max,'inf');


%Genera la Tabla
fprintf('\n\t%s \t\t%s \t\t%s\n' ,'Método', 'Iteraciones','Error');
fprintf('\t------------------------------------------');

fprintf('\n\t%s \t\t%d \t\t\t%.4d', 'pCauchy', iter1, error1);
fprintf('\n\t%s \t\t%d \t\t\t%.4d', 'pDogleg', iter2, error2);
fprintf('\n');
fprintf('\t------------------------------------------');
fprintf('\n');