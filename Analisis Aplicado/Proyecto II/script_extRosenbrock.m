%Script para funcion Rosenbrock
clear all; close all; clc;

%Variables
% n = 10;
% x0 = [-1.2;1];
% x0 = repmat(x0,n/2,1);
f = @extendedRosenbrock;

%Parametros
tol = 1e-5;
maxIter = 1000;
deltaMax = 1.5;

% [x1, iter1] = lineBGFS(@extendedRosenbrock, x0, tol, maxIter);
% [x2, iter2] = mRCSR1(@extendedRosenbrock, x0, tol, maxIter,deltaMax);
% [x3, iter3] = lineBGFSLM(@extendedRosenbrock, x0, tol, maxIter, 3);

% Genera la tabla para lineBGFS
fprintf('\t\t\t\t\t\t%s' ,'Funcion Rosenbrock')
fprintf('\n\n\t\t%s' ,'lineBGFS')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'n', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for n = 2:200
    if(n == 2||n == 10 )
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFS(f, x0, tol, maxIter);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t\t%.2i \t\t\t\t%.16f', n, iter,t2);
    end
    
    if(n == 100 || n == 200 )
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFS(f, x0, tol, maxIter);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t%.16f', n, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

% Generamos la tabla para mRCSR1
fprintf('\n\t\t%s' ,'mRCSR1')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'n', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for n = 2:200
    
    if(n == 2 )
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = mRCSR1(f, x0, tol, maxIter,deltaMax);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t\t%.2i \t\t\t\t%.16f', n, iter,t2);
    end
    
    if(n == 10 )
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = mRCSR1(f, x0, tol, maxIter,deltaMax);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t\t%.2i \t\t\t%.16f', n, iter,t2);
    end
    
    if(n == 100 || n == 200 )
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = mRCSR1(f, x0, tol, maxIter,deltaMax);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t%.16f', n, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

% Generamos la tabla para lineBGFSLM , n = 2
n = 2;
fprintf('\n\t\t%s' ,'lineBGFSLM, n = 2 ')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 27)
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');


% Generamos la tabla para lineBGFSLM , n = 10
n = 10;
fprintf('\n\t\t%s' ,'lineBGFSLM, n = 10 ')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 27)
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');


% Generamos la tabla para lineBGFSLM , n = 100
n = 100;
fprintf('\n\t\t%s' ,'lineBGFSLM, n = 100 ')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 27)
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

% Generamos la tabla para lineBGFSLM , n = 200
%obs con m = 17 entra en un loop
n = 200;
fprintf('\n\t\t%s' ,'lineBGFSLM, n = 200 ')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29
    if(m == 1 || m == 3 || m == 5 || m == 27)
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');

% Generamos la tabla para lineBGFSLM , n = 1000
n = 1000;
fprintf('\n\t\t%s' ,'lineBGFSLM, n = 1000 ')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s\n', 'm', 'Numero de iteraciones', 'Tiempo');
fprintf('\t--------------------------------------------------------------------------');
for m = 1:29
    if(m == 1 || m == 3 || m == 5 || m == 17 || m == 27)
        x0 = [-1.2;1];
        x0 = repmat(x0,n/2,1);
        t1 = tic;
        [x, iter] = lineBGFSLM(f, x0, tol, maxIter, m);
        t2 = toc(t1);
        fprintf('\n\t%.2i \t\t\t\t\t\t%.2i \t\t\t\t%.16f', m, iter,t2);
    end
end
fprintf('\n');
fprintf('\t--------------------------------------------------------------------------');
fprintf('\n');