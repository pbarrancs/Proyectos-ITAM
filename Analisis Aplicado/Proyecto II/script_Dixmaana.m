%Script para funcion Dixmaana
clear all; close all; clc;

%Parametros
tol = 1e-5;
maxIter = 100;
deltaMax = 1.5;

%Variables
% n = 10;
% x0 = 2 * ones(n,1);
f = @Dixmaana;
% 
% [x1, iter1] = lineBGFS(f, x0, tol, maxIter);
% [x2, iter2] = mRCSR1(f, x0, tol, maxIter, deltaMax);
% [x3, iter3] = lineBGFSLM(f, x0, tol, maxIter, 1);


% Generamos la tabla para lineBGFS, n = 200;
fprintf('\t\t\t\t\t\t%s' ,'Funcion Dixmaana')
fprintf('\n\n\t\t%s' ,'lineBGFS')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'n', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
    n = 200;
    x0 = 2 * ones(n,1);
    t1 = tic;
    [x, iter] = lineBGFS(f, x0, tol, maxIter);
    t2 = toc(t1);
    fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', n, iter,t2);
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

% Generamos la tabla para mRCSR1
fprintf('\n\t\t%s' ,'mRCSR1')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'n', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
    n=200;
    x0 = 2 * ones(n,1);
    t1 = tic;
    [x, iter] = mRCSR1(f, x0, tol, maxIter,deltaMax);
    t2 = toc(t1);
    fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', n, iter,t2);
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');


% Genera la tabla lineBFGSLM para n = 200
n = 200;
fprintf('\n\t\t%s' ,'lineBGFSLM con n = 200')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29    
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 29)
        x0 = 2 * ones(n,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');


% Genera la tabla lineBFGSLM para n = 1000
n = 1000;
fprintf('\n\t\t%s' ,'lineBGFSLM con n = 1000')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29    
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 29)
        x0 = 2 * ones(n,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

