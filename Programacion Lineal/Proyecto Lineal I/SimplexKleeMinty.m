clear; close all; clc;

% Genera la tabla
fprintf('\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
tiempos = zeros(10,1);
for m = 3:10
    [A,b,c] = generaKleeMinty(m);
    t1 = tic;
    [xo, zo, ban, iter] = mSimplexFaseII(A,b,c);
    t2 = toc(t1);
    tiempos(m)= t2;
    fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');