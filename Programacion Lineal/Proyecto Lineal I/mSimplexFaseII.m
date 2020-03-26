function [xo, zo, ban, iter] = mSimplexFaseII(A,b,c)
% purpose: Version del Simplex en la Fase II
% minimizar c^T x
% sujeto a Ax <= b , x >= 0 , b >= 0
% In : A ... mxn matrix
%      b ... column vector with as many rows as A
%      c ... column vector with as many entries as one row of A
%
% Out: xo ... SFB optima del problema
%      zo ... valor optimo del problema
%      ban ... indica casos:
%              -1 ... si el conjunto factible es vacio
%               0 ... si se encontro una solucion optima
%               1 ... si la funcion objectivo no es acotada.
%   iter ... es el numero de iteraciones (cambios de variables basicas)
%            que hizo el metodo

    %Dimension del problema
    [m,n] = size(A);%donde m es el numero de rengoles y n el de columnas

    %Inicializamos nuestra variable que contara las iteraciones
    iter = 0;

    % Estandarizamos el problema añadiendo variables de holgura
    A = [A,eye(m)];
    c = [c;zeros(m,1)]; %vector de costos con variables de holgura
    
    N = 1:n;% Define el conjunto de los indices de las variables no basicas 
    B = n+1:n+m;%Define el conjunto de los indices de las variables basicas
    
    % Fase II revisada 
    if any(b<0) % condicion de que las b>0 para que el problema tenga conjunto factible
        %fprintf('El conjunto factible es vacio\n\n');
        xo = 'NA';
        zo = 'NA';
        ban = -1; %la funcion objectivo no es acotada.
    else
        while 1
        Ab = A(:,B); %Submatriz basica
        An = A(:,N); %Submatriz no basica
        cb = c(B); %Vector de costos basico
        cn = c(N); %Vector de costos no basico
        
        
        % Paso 1
        lambda = (inv(Ab))'*cb; %
        rnt = (lambda')*An - (cn'); %vector de ahorro 
        if any(rnt > 0) %De lo contrario ya termino
            [~,e] = max(rnt); % Mayor descenso. Busca el indice mayor de r
            
            
            % Paso 2
            h = inv(Ab)*b; %Pivoteamos
            He = inv(Ab)*An(:,e); %Renglon de entrada
            if any(He > 0) %De lo contrario el sistema es no acotado
                indicesP = find(He > 0); %encuentra los indices de las entradas positivas de He
                [~,im] = min(h(indicesP)./He(indicesP)); %encuentra el indice que minimiza h/He
                s = indicesP(im); %Variable de salida
                
                
                % Paso 3
                aux = B(s); %acutualizamos los indices de las basicas y no basucas
                B(s) = N(e);
                N(e) = aux;
            else
                %fprintf('El problema es no acotado\n\n');
                xo = 'NA';
                zo = 'NA';
                ban = 1; %el conjunto factible es vacio
                break
            end
        else
            xo = [zeros(m+n,1)];
            x1 = inv(Ab)*b;
            xo(B)= x1; %SBF del problema
            zo = (lambda')*b; %valor optimo del problema
            ban = 0; %se encontro una solucion optima
            break
        end
        iter = iter + 1; %Contabilizamos las iteraciones
        end
    end
return