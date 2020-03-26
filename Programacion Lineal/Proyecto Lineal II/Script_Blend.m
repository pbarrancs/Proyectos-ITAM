close all; clear all; clc;
load('blend.mat');
A = full(A);
[x0, lam0, z0, iter] = mPI_Newton(c, A, b);
fprintf('\nScript Problema Blend');
%Valor optimo es
z = (c')*x0;
fprintf('\n\tValor optimo: %1.6f', z);
%Numero de Iteraciones
fprintf('\n\tNumero de iteraciones: %i' , iter);
%El valor de ||diag(x)z||
d = norm(diag(x0)*z0,inf);
fprintf('\n\tEl valor de la norma infinito de ||diag(x)z|| es: %1.16f\n', d);