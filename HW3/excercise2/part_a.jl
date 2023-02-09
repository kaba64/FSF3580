using LinearAlgebra
using Plots
using LaTeXStrings

n = 5;
x = randn(n);
y = randn(n);
y = y;
i=argmax(abs.(y));
rho = -sign(y[i]);
alpha = (rho*norm(x))/norm(y);
z=x-alpha*y;
u = z/norm(z);
p = Matrix(I,n,n)-2*u*u';
display(norm(p*x-alpha*y));
