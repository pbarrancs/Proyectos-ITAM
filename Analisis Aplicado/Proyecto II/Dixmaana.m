function [f] = Dixmaana(x)

n = length(x);
m = floor(n/3);

%% Definimos los parametros de la funcion
alpha = 1;
beta = 0;
gamma = 0.125;
delta = 0.125;
k1 = 0;
k2 = 0;
k3 = 0;
k4 = 0;

%% Partimos la funcion en cuatro sumas
term1 = 0;
term2 = 0;
term3 = 0;
term4 = 0;

%%
for i = 1:n
    term1 = term1 + alpha*(x(i)^2)*(i/n)^k1;
    if i <= n-1
        term2 = term2 + beta*(x(i)^2)*((x(i+1) + x(i+1)^2)^2)*(i/n)^k2;
    end
    if i <= 2*m
        term3 = term3 + gamma*(x(i)^2)*(x(i+m)^4)*(i/n)^k3;
    end
    if i <= m
        term4 = term4 + delta*x(i)*x(i+2*m)*(i/n)^k4;
    end     
end

f = 1 + term1 + term2 + term3 + term4;

end
