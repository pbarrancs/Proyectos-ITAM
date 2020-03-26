function [x,y,u,z,iter,h,bag,ftime] = qpintpointpc3(Q, A, F, c, b, d)
% Resuelve el problema cuadr´atico por el m´etodo de puntos interiores.
% Número maximo de iteraciones permitidas es, maxiter = 100.
% La variable de salida,iter, es el n´umero de iteraciones que usa el metodo.
% Se inicia con iter = 0 y se incrementa en uno en cada itercai´on.
% Se detiene el metodo cuando ?Fo(wk)?2 ? 10?5 o k = maxiter.
%
%Método de puntos interiores para
%   Min (1/2)*Q*x + c'*x
%   s.a. A*x = b
%   F*x >= d
%
%------------------------------------------------
n = length(c);
m = length(b);
p = length(d);

maxiter = 100;
tol = 1.e-05;
iter = 0;

x = A\b;
y = zeros(m,1);
u = ones(p,1);
z = ones(p,1);
gamasize = (1/2)*(u'*z)/p;

v1 = Q*x + c + A'*y - F'*u;
v2 = A*x - b;
v3 = -F*x + z + d;
v4 = diag(u)*diag(z)*ones(p,1);

v = [v1;v2;v3;v4];

etime = cputime;
while(norm(v) > tol && iter < maxiter)
    %Resolvemos el primer sistema
    K = [Q              A'         -F'         zeros(n,p);
         A           zeros(m)    zeros(m,p)    zeros(m,p);
         -F          zeros(p,m)  zeros(p)      eye(p);
         zeros(p,n)  zeros(p,m)  diag(z)       diag(u)];
     ld = -[v1; v2; v3; v4];
     
     Dw = K\ld;
%      Dx = Dw(1:n);
%      Dy = Dw(n+1:n+m);
     Du = Dw(n+m+1:n+m+p);
     Dz = Dw(n+m+p+1:n+m+p+p);
     %----------------------------------------------------------
     %Recortar el paso
     vu = ones(p,1);
     vz = ones(p,1);
     
     for i = 1:p
         if(Du(i) < 0)
             vu(i) = -(u(i)/Du(i));
         end
         
         if(Dz(i) < 0)
             vz(i) = -(z(i)/Dz(i));
         end
     end
     %sale el valor alfa
     alfau = min(vu);
     alfaz = min(vz);
     alfa = (0.995)*min([alfau alfaz 1]);
     
     %Actualizar mu* y z*
%      
%      x = x + alfa*Dx;
%      y = y + alfa*Dy;
     u0 = u + alfa*Du;
     z0 = z + alfa*Dz;
     
%      u0 = u + alfa*Du;%u0
%      z0 = z + alfa*Dz;%z0
     
     %Actualizamos gama
     gama = (1/2)*(z0'*u0)/p;%u0,z0
     gama = (gama/gamasize)^3;
     gama = gama*gamasize;
     
     
     %Duda Tenemos que reactualizar w y K
     %Resolvemos el segundo Sistema lineal 
     
%      v1 = Q*x + c + A'*y - F'*u;
%      v2 = A*x - b;
%      v3 = -F*x + z + d;
%      v4 = diag(u)*diag(z)*ones(p,1);

     
%      ld2 = -[v1; v2; v3; v4 + diag(Dz)*diag(Du)*ones(p,1) - gama*ones(p,1)]; %Duda es Dz o z*
%      ld2 = -[v1; v2; v3; v4 + Dz*Du'*ones(p,1) - gama*ones(p,1)];
     
%     ld2 = ld + [zeros(n,1);zeros(m,1);zeros(p,1); -diag(Dz)*diag(Du)*ones(p,1) + gama*ones(p,1)];
    ld2 = ld + [zeros(n,1);zeros(m,1);zeros(p,1); -diag(z0)*diag(u0)*ones(p,1) + gama*ones(p,1)];
    
     Dw2 = K\ld2;
     Dx2 = Dw2(1:n);
     Dy2 = Dw2(n+1:n+m);
     Du2 = Dw2(n+m+1:n+m+p);
     Dz2 = Dw2(n+m+p+1:n+m+p+p);
     
     %--------------------------------------------------------------
     % Hacemos el segundo recorte
     vu = ones(p,1);
     vz = ones(p,1);
     
     for i = 1:p
         if(Du2(i) < 0)
             vu(i) = -(u(i)/Du2(i));
         end
         
         if(Dz2(i) < 0)
             vz(i) = -(z(i)/Dz2(i));
         end
     end
     alfau = min(vu);
     alfaz = min(vz);
     alfa = (0.995)*min([alfau alfaz 1]);
     %--------------------------------------------------------------
     
     
     %Actualizar
     x = x + alfa*Dx2;
     y = y + alfa*Dy2;
     u = u + alfa*Du2;
     z = z + alfa*Dz2;
     %Condiciones KKT
     v1 = Q*x + c + A'*y - F'*u;
     v2 = A*x - b;
     v3 = -F*x + z + d;
     v4 = diag(u)*diag(z)*ones(p,1);
     
     v = [v1;v2;v3;v4];
     iter = iter + 1;
     bag(iter) = norm(v);
     %disp(sprintf('%2.0f %2.10f', iter, norm(v)))
end

ftime = cputime - etime;
nbag = length(bag);
iter;
bag;
semilogy(1:nbag, bag, 'b', 1:nbag, bag, '--r', 'Linewidth', 3)
h = (1/2)*x'*Q*x + c'*x;
end

