clear all

%Script: afiro_qp.m

load afiro
A = full(A);
n = length(c);
Q = eye(n); F = eye(n); d = zeros(n,1);

%[x, y, u, z, iter,h] = qpintpointpc3(Q, A, F, c, b,d);
%[x,y,u,z,iter,h] = qpintpointpc(Q, A, F, c, b, d);
%[x,y,u,z,iter,h] = qpintpointpc2(Q, A, F, c, b, d);
%[x,y,u,z,iter,h] = qpintpointpc4(Q, A, F, c, b, d)