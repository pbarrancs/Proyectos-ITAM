function [x0, lam0, z0, iter] = mPI_Newton(c, A, b)
% Resuelve el problema primal y dual por el metodo descrito arriba.
% Haz a lo mas 200 iteraciones.
% En cada iteración muestre el número de iteración y la norma de F(wk)
%
% In : A ... mxn matrix
%      b ... vector columna con tantas filas como A
%      c ... vector columna con tantas columnas como A
%
% Out: xo ... solución  óptima del problema primal
%       lamo, zo ... solución  óptima del problema dual
%       iter ... el numeró de iteraciones que hizo el método
%
% Parametros:
%        tol = 1e-9; x0 = ones(n,1); lam0 = zeros(m,1); z0 = ones(n,1);

%Declaramos sus parametros e inicializamos las variables
[m,n] = size(A);
tol = 1e-9;
x0 = ones(n,1);
lam0 = zeros(m,1);
z0 = ones(n,1);

iter=0;
sigma = 0.1;
e = ones(n,1);

F = [A*x0 - b; (A')*lam0 + z0 - c; diag(x0)*diag(z0)*e];

while norm(F,inf) > tol && iter<200
        mu = (1/n)*(x0'*z0);
        J = [A,zeros(m,m),zeros(m,n);zeros(n,n),A',eye(n,n);diag(z0),zeros(n,m),diag(x0)];
        delta = J\(-F+[zeros(n,1);zeros(m,1);sigma*mu*e]);
        
        
        dx = delta(1:n,1);
        dlam = delta(n+1:n+m,1);
        dz = delta(m+n+1:m+2*n);

        ax = recorte(x0,dx);
        az = recorte(z0,dz);

        x0 = x0 + (999/1000)*(ax*dx);
        lam0 = lam0 + (999/1000)*(az*dlam);
        z0 = z0 + (999/1000)*(az*dz);

        F = [A*x0 - b; (A')*lam0 + z0 - c; diag(x0)*diag(z0)*e];

            if ax*az > 0.8
                sigma = max(10^-4,sigma/10);
            end
        iter = iter+1;
    end
end

