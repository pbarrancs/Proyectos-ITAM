

% % generar dimensiones del problema
% m = round(10*exp(log(20)*rand()));
% n = round(10*exp(log(20)*rand()));
% % generar A, b, c
% sigma = 100;
% A = round(sigma*randn(m,n));
% b = round(sigma*abs(randn(m,1)));
% c = round(sigma*randn(n,1));

n = zeros(1,50);
m = zeros(1,50);
x = zeros(1,50);
ban = zeros(1,50);
iter = zeros(1, 50);
for i =1:50
        % generar dimensiones del problema
    m(i) = round(10*exp(log(20)*rand()));
    n(i) = round(10*exp(log(20)*rand()));
    
        % generar A, b, c
    sigma = 100;
    A = round(sigma*randn(m(i),n(i)));
    b = round(sigma*abs(randn(m(i),1)));
    c = round(sigma*randn(n(i),1));
    
        %Llamamos a nuestra funcion
    [xo, zo, ban(i), iter(i)] = mSimplexFaseII(A,b,c);
    x(i) = min(m(i),n(i));
end


%Hacemos la Grafica 
acotado = find(ban==0) % vector que tiene en sus entradas los indices acotados
noAcotado = find(ban ==1) % vector que tiene en sus entradas los indices no acotados
scatter( x(acotado), iter(acotado), 'b', 'filled')
hold on
scatter( x(noAcotado), iter(noAcotado), 'r', 's', 'filled')
hold off


xlabel('min(m,n)', 'fontsize', 14);
ylabel('#it', 'fontsize', 14);
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
grid on

