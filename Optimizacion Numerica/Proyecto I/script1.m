% Script que le enseña a mike a escribir
x = [2;3;4;5];
n = length(x);
nombres = ["pepe"; "mike"; 1; 2];
fprintf("\nProblema \t|\t n \t | \t m \t | \t iter \t | \t f(x*)\n")
for i =1:length(x)
    fprintf("\n %s \t|\t %f.16 \t | \t m \t | \t iter \t | \t f(x*)\n",nombres(i),x(i));
end

