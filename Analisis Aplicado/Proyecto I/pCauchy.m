function [pC] = pCauchy( B, g, delta )
% In :  B ... (Matriz simétrica) Aproximación de la hessiana de f en xk
%       g ... (vector) gradiente de f en xk
%   delta ... radio de región de confianza
%
% Out : pC ... punto de Cauchy

%Por razones de eficiencia, guardamos el término cuadrático y la norma de
%g pues no queremos hacer llamadas innecesarias a funciones
termCuadr = g'*B*g;
normg = norm(g);

%% Calculamos alfak 
%Si el termino cuadrático es mayor a 0
%obtenemos el mínimo entre (normg.^3)/(delta*termCuadr),
%en caso contrari alfak es 1.
alfak = 1;
    if termCuadr > 0
        alfa = (normg.^3)/(delta*termCuadr);
        alfak = min(alfa,1);
    end

%% Calculamos pk
pk = -(delta/normg) * g;

%% Obtenemos el punto de Cauchy
pC = 0.99*alfak*pk;

end

