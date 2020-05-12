function [f] = fbooth()
%Funcion Booth
%Entradas: x vector de dos coordenadas

f =  @(x) (x(1) + 2*x(2) - 7)^2 + (2*x(1) + x(2) - 5)^2;
% Df, Hf como posibles par�metros extra
% Df = @(x) [10*x(1) + 8*x(2) -34; ...
%             8*x(1) + 10*x(2) -38];
% 
% Hf = [10 8; 8 10];
%
% Como el c�digo en general calula el gradiente y la hessiana no s� si es
% buena idea calcularlos aqu�

end


