%Script: AFIRO

clear all

load sctap1
A = full(A);
n = length(c);
Q = eye(n); F = eye(n); d = zeros(n,1);

disp('Opción 1')

figure(1)
[x1,y1,u1,z1,iter1,h1,bag1,ftime1] = qpintpointpc(Q, A, F, c, b, d);

fprintf('El valor óptimo es : %f\n', h1);
fprintf('El número de iteraciones fue : %f\n', iter1);
fprintf('El tiempo transcurrido fue : %f\n', ftime1);
fprintf('La norma de las condiciones suficientes de primer orden en el punto final es : %f\n', bag1(iter1));

disp('Opción 2')

figure(2)
[x2,y2,u2,z2,iter2,h2,bag2,ftime2] = qpintpointpc2(Q, A, F, c, b, d);

fprintf('El valor óptimo es : %f\n', h2);
fprintf('El número de iteraciones fue : %f\n', iter2);
fprintf('El tiempo transcurrido fue : %f\n', ftime2);
fprintf('La norma de las condiciones suficientes de primer orden en el punto final es : %f\n', bag2(iter2));

disp('Opción 3')

figure(3)
[x3,y3,u3,z3,iter3,h3,bag3,ftime3] = qpintpointpc3(Q, A, F, c, b, d);

fprintf('El valor óptimo es : %f\n', h3);
fprintf('El número de iteraciones fue : %f\n', iter3);
fprintf('El tiempo transcurrido fue : %f\n', ftime3);
fprintf('La norma de las condiciones suficientes de primer orden en el punto final es : %f\n', bag3(iter3));

disp('Opción 4')

figure(4)
[x4,y4,u4,z4,iter4,h4,bag4,ftime4] = qpintpointpc4(Q, A, F, c, b, d);

fprintf('El valor óptimo es : %f\n', h4);
fprintf('El número de iteraciones fue : %f\n', iter4);
fprintf('El tiempo transcurrido fue : %f\n', ftime4);
fprintf('La norma de las condiciones suficientes de primer orden en el punto final es : %f\n', bag4(iter4));

fprintf('\n\t\t%s' ,'Resultados SCTAP1')
fprintf('\n\n\t%s \t\t\t%s \t\t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n', 'Opcion', 'Numero de iteraciones', 'Tiempo', 'f(x*)', 'Norm(V)');
fprintf('\n\t--------------------------------------------------------------------------------------------------');
    fprintf('\n\t%s \t\t%.2i \t\t\t\t\t\t\t\t%.2i \t\t%.4f \t\t\t%f' ,'Opcion 1' , iter1,ftime1,h1,bag1(iter1));
    fprintf('\n\t%s \t\t%.2i \t\t\t\t\t\t\t\t%.2i \t\t%.4f \t\t\t%f ','Opcion 2' , iter2,ftime2,h2,bag2(iter2));
    fprintf('\n\t%s \t\t%.2i \t\t\t\t\t\t\t\t%.2i \t\t%.4f \t\t\t%f','Opcion 3' , iter3,ftime3,h3,bag3(iter3));
    fprintf('\n\t%s \t\t%.2i \t\t\t\t\t\t\t\t%.2i \t\t%.4f \t\t\t%f','Opcion 4' , iter4,ftime4,h4,bag4(iter4));
fprintf('\n');
fprintf('\t----------------------------------------------------------------------------------------------------');
fprintf('\n');