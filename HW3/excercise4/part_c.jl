using LinearAlgebra
using Random
using BenchmarkTools
using Plots

function Jordan(eps)
    A = [π 1;0 π+eps];
    λeA, veA = eigen(A);
    Fj = veA * Diagonal(exp.([λeA[1],λeA[2]])) * veA^(-1);
    return Fj;
end

function analytical(eps)
    A = [π 1;0 π+eps];
    λeA = eigvals(A);
    α = ((π+eps)*exp(π)-π*exp(π+eps))/eps;
    β = (exp(π+eps)-exp(π))/eps;
    Fa = α*I + β*A;
    return Fa;
end
n = 120;
eps_ar = zeros(n);
error_r  = zeros(n);
alpha = 0.4^(1/5);
for i in 1:n
    eps_ar[i] = alpha^i;
end

for i in 1:n
    Fj = Jordan(eps_ar[i]);
    Fa = analytical(eps_ar[i]);
    error_r[i] = norm(Fj-Fa)
end

t = plot(eps_ar,error_r,xaxis=:log,yaxis=:log,legend=false);
xlabel!("ϵ")
ylabel!("||f(A)-F||");
savefig(t,"excercise_4c.png");
