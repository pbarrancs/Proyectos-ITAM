function [f] = extendedRosenbrock(x)
n = length(x);
k = 100;
f = 0;
     for i = 1:n/2
         f = f + k*(x(2*i)-x(2*i-1)^2)^2 + (1-x(2*i-1))^2;
     end
end