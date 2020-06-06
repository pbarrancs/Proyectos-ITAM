clear all ;clc;close all;
warning('off','all'); % Obtenemos pocas matrices malcondicionadas

%% Generamos nuestros punto inicial
prompt = 'Introduzca numero de vertices (np) y teclee Enter: ';
np = input(prompt);
thetas = linspace(1,np,np)*(pi/np);
rs = ones(1,np)*0.2;
rs(1) = 0;
thetas(np) = 0;      
thetas = sort(thetas);
x0 = [rs,thetas]';
xi = rs.*cos(thetas);
yi = rs.*sin(thetas);
xi = [xi, 0];
yi = [yi, 0];


%% Nuestra Solucion
funobj = 'funcObj';
fung = 'rest';
etime = cputime;
[x, fx, iter, cnpo] = metpuntint(funobj, fung, x0);
timeN = cputime - etime;

r_sol = [x(1:np);x(1)];
theta_sol = [x(np+1:2*np);x(np+1)];
x_sol = r_sol.*cos(theta_sol);
y_sol = r_sol.*sin(theta_sol);



%% Solucion Matlab
etime = cputime;
[xm,fxm,exitflag,output] = fmincon(funobj,x0,[],[],[],[],[],[],'restmatlab');
timeM = cputime - etime;
rm_sol = [xm(1:np);xm(1)];
thetam_sol = [xm(np+1:2*np);xm(np+1)];
xm_sol = rm_sol.*cos(thetam_sol);
ym_sol = rm_sol.*sin(thetam_sol);


%% Ploteo de las Soluciones
figure(1)
plot(xi,yi, ':') %puntos iniciales
hold on
plot(x_sol, y_sol) %metpuntint
hold on
plot(xm_sol,ym_sol,'--') %fmincon
hold off
title('Resultados') 
legend({'Cond inicial','metpuntint','fmincon'},'Location','northeast')


%% Comparaciones
fprintf('\n\t%s \t\t%s \t\t%s \t\t\t%s \t\t\t%s\n' ,'Método', 'Iteraciones','Tiempo','Cond KKT','f(x)');
fprintf('\t-----------------------------------------------------------------------------------');
fprintf('\n\t%s \t\t%d \t\t%.4d  \t\t%.4d  \t\t%.4d', 'metpuntint', iter, timeN, cnpo,-fx );
fprintf('\n\t%s \t\t%d \t\t%.4d  \t\t%.4d  \t\t%.4d', 'fmincon',  output.iterations, timeM,output.firstorderopt,-fxm);
fprintf('\n');
fprintf('\t-----------------------------------------------------------------------------------');
fprintf('\n');