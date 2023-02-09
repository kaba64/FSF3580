using LinearAlgebra
using Random
using SparseArrays
using BenchmarkTools
using LaTeXStrings

include("Benchmark_GMRES.jl")
gr();
n=500;
alpha = 100;
Random.seed!(0);
b = rand(n);
#m_iter = [5,10,20,50,100];
m = 100;
method = "B";
if(method=="G")
    A = sprand(n,n,0.5);
    A = A + alpha*sparse(I, n, n);
    A=A/norm(A,1);
    b_g = @benchmarkable Benchmark_GMRES(A,b,m);
    tune!(b_g);
    g_t = run(b_g);
    println("GMRES timing with m=$m and n=$n");
    display(g_t);
    x = Benchmark_GMRES(A,b,m);
    println("GMRES residual with m=$m and n=$n");
    display(norm(A*x - b));
else
# ********* Backslash ******************
    b_b = @benchmarkable A\b;
    tune!(b_b);
    b_t = run(b_b);
    println("Backslash timing with n=$n");
    display(b_t);
    xbs = A\b;
    println("Backslash residual with n=$n");
    display(norm(A*xbs - b));
end
