function  [gfx] = gradiente(fname,x)
% Aproximacion del gradiente por diferencias hacia adelante 
% de una funcion  de R^n a R.

%In
% fname .- cadena con el nombre de la funcion.
% x .- vector columna de dimension n.
% Out
% gfx - vector columna de dimension n, es la aproximacion al
%       gradiente en x.
%--------------------------------------------------------------

h = 1.e-05;                    % tama�o de paso
fx = feval(fname,x);           % se eval�a la funci�n en el punto
n = length(x);                 % dimensi�n del vector
gfx = zeros(n,1);              % vector gradiente de salida
xt = x;

for i = 1:n
    xt(i) = xt(i) +h ;          % punto en la direcci�n i
    fxh = feval(fname,xt);      % evaluaci�n en la direcci�n i
    gfx(i) = ( fxh - fx)/ h;    % diferencia hacia adelante
    xt(i) = x(i);               % regreso a la copia de x
end



