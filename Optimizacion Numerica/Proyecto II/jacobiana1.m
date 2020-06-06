function [Jx] = jacobiana1(fname,x)
% Esta funcion aproxima por diferencia finitas la
% matriz jacobiana de fname en x, donde
% fnama:R^n --> R^m.

% In
% fname.- cadena con el nombre de la funcion.
% x    .- vector de dimension n.
% Out
% Jx.- matriz  de orden nxm.

h = 1.e-06;

n = length(x);
Fx = feval(fname,x);
m = length( Fx);
Jx = zeros(m,n);


for k = 1:n
    xh = x; 
    xh(k) = xh(k) + h;
    Fxh = feval(fname,xh);
    Jx(:,k) = (Fxh - Fx)/ h;
end
