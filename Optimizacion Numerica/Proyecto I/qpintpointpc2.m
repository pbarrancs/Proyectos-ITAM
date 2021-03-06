function [x,y,u,z,iter,h,bag,ftime] = qpintpointpc2(Q, A, F, c, b, d)
% Resuelve el problema cuadr�atico por el m�etodo de puntos interiores.
% N�mero maximo de iteraciones permitidas es, maxiter = 100.
% La variable de salida,iter, es el n�umero de iteraciones que usa el metodo.
% Se inicia con iter = 0 y se incrementa en uno en cada itercai�on.
% Se detiene el metodo cuando Fo(wk)2 ? 10?5 o k = maxiter.
%
%M�todo de puntos interiores para
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
gama = 1.0;

v1 = Q*x + c + A'*y - F'*u;
v2 = A*x - b;
v3 = -F*x + z + d;
v4 = diag(u)*diag(z)*ones(p,1);

v = [v1;v2;v3;v4];
r1 = -v1; r2 = -v2; r3 = -v3; r4 = -v4;
ry = 0; rx = 0;

etime = cputime;
while(norm(v) > tol && iter < maxiter)
%     QQ = Q + F'*(diag(z)\diag(u))*F; 
%     r = r1 +F'\diag(z)*((r4 + gama*ones(p,1)) -(diag(u)*r3));
    QQ = Q + F'*(diag(z)\diag(u))*F; 
%     r = r1 +F'*diag(z)\((r4 + gama*ones(p,1)) -(diag(u)*r3));
    r = r1 +F'*diag(z)\(r4  -(diag(u)*r3));
    
    K = [QQ              A'    ;
         A           zeros(m)  ];
     ld = [r; r2];
%      ld = -[r;v3];
     
     if(cond(K) <= 10^5)
         ry = ry + 1 ;
         Dw = K\ld;  %%Cambio Wk
         Dx = Dw(1:n);
         Dy = Dw(n+1:n+m);
         Du = diag(z)\r4 - diag(z)\diag(u)*r3 - diag(z)\diag(u)*F*Dx;
         Dz = r3 + F*Dx;
     else
         rx = rx + 1;
         K = [Q              A'         -F'         zeros(n,p);
             A           zeros(m)    zeros(m,p)    zeros(m,p);
             -F          zeros(p,m)  zeros(p)      eye(p);
             zeros(p,n)  zeros(p,m)  diag(z)       diag(u)];
         ld = -[v1; v2; v3; v4 - gama*ones(p,1)];
         
         Dw = K\ld;
         Dx = Dw(1:n);
         Dy = Dw(n+1:n+m);
         Du = Dw(n+m+1:n+m+p);
         Dz = Dw(n+m+p+1:n+m+p+p);
     end

    %----------------------------------------------------------
    % Recortar el paso
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
     alfau = min(vu);
     alfaz = min(vz);
     alfa = (0.995)*min([alfau alfaz 1]);

     
     %sale el valor alfa
     %----------------------------------------------------------
     %Actualizar
     x = x + alfa*Dx;
     y = y + alfa*Dy;
     u = u + alfa*Du;
     z = z + alfa*Dz;
     
     %Condiciones KKT
     v1 = Q*x + c + A'*y - F'*u;
     v2 = A*x - b;
     v3 = -F*x + z + d;
     v4 = diag(u)*diag(z)*ones(p,1);
     r1 = -v1; r2 = -v2; r3 = -v3; r4 = -v4;
     v = [v1;v2;v3;v4];
     gama = (1/2)*(u'*z)/p;
     iter = iter + 1;
     bag(iter) = norm(v);
     %disp(sprintf('%2.0f %2.10f', iter, norm(v)))
end

ftime = cputime - etime;
iter;
ry;
rx;
bag;
nbag = length(bag);

semilogy([1:nbag], bag, 'b', [1:nbag], bag, '--r', 'Linewidth', 3)
h = (1/2)*x'*Q*x + c'*x;
end

